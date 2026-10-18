#!/usr/bin/env python3
"""Select the unique protein chains of every mmCIF file in a folder, rebuild their gaps,
renumber and centre them, and write each one with its domain and family annotation.

For every *.cif.gz (or *.cif) file in mmcif_files_location:
  1. Read the structure, keep the first model and drop hydrogens. Alternate locations:
     at each position keep the residue type and the conformer with the highest
     occupancy; of overlapping whole-residue alternates (two different residues, each on
     its own altloc) keep the better-occupied one.
  2. Chains are author chains (auth_asym_id). Protein chains (polypeptide polymer) are
     processed; chains with no ATOM records (only HETATM entries) and nucleic-acid chains
     are only listed in the summary.
  3. Signature of a protein chain:
       ATOM-record sequence - polymer residues with coordinates, in file order, residue
           names as written (modified residues such as SEP or MSE keep their names)
       HETATM groups of the chain whose code is in hetatm_list, counted by code
       waters of the chain with a structural/functional role, counted by role
     Chains with identical signatures form one unique chain. Its representative is the
     member with the most observed heavy atoms, then the lowest mean B, then the first.
  4. For each representative:
       - rebuild internal gaps of up to MAX_GAP_LENGTH residues (from SEQRES) and missing
         heavy atoms of standard residues with PDBFixer/OpenMM, one small local fragment
         at a time; added atoms get occupancy 0.00; unmodelled terminal tails stay out
       - renumber: polymer residues 1..N, then the HETATM groups, then the waters
       - move the centre of mass (all written atoms, atomic masses) to the origin
       - write <pdb>_<chain>/<pdb>_<chain>.cif, .pdb and _residue_map.csv (every SEQRES
         position with its author number, new number and status observed/rebuilt/missing,
         then the HETATM groups and waters with their old and new numbers)
  5. One row per chain in unique_chain_summary.csv.

Domain mapping is a separate, quick step (--map-domains) that can run any time after the
processing: for every chain folder it reads the residue map, writes <pdb>_<chain>_domains.csv
with the Pfam, CATH, SCOP2 and UniProt ranges in the new numbering, and fills the domains
column of the summary. Rerun it after the domain database is updated; the structures are
not processed again.

Water roles (each water counted once, under the first role it meets; heavy atoms of the
asymmetric unit; partners are polymer atoms of any chain and atoms of listed HETATM groups):
  metal       first shell of a metal atom of a listed group (free ion or cofactor metal,
              e.g. Fe of HEM): deposited metalc link or distance <= METAL_CUTOFF
  ligand      H-bond (<= HBOND_CUTOFF) to a listed group's N/O/S and to the polymer or
              to a second listed group
  buried      at most MAX_EXPOSED of its solvent-accessible surface open (1.4 A probe,
              other waters removed) and at least one H-bond to the polymer or a listed group
  interface   H-bonds to the polymers of two or more chains
  structural  >= MIN_HBONDS H-bonds to polymer N/O atoms of at least two residues

A rebuilt loop is kept only if every peptide bond across it is <= MAX_PEPTIDE_BOND; it gets
three tries (different random seeds), otherwise the gap stays open and is noted.

Domain and family data come from bulk files, downloaded once with --prepare-domain-data
into domain_data_location and turned into one SQLite database; --map-domains only reads
that database (no network):
  UniProt   SIFTS pdb_chain_uniprot (ranges) + UniProt query of all entries with a 3D
            structure (entry and protein names)
  Pfam      SIFTS pdb_pfam_mapping (ranges, names) + Pfam-A.clans (descriptions)
  CATH      CATH-B daily boundaries (author numbering) + cath-names
  SCOP2     SIFTS pdb_chain_scop2b_sf_uniprot (superfamily domains with ranges) and
            pdb_chain_scop2_uniprot (family per chain) + SCOP2 cla/des files (names)
The InterPro files (protein2ipr.dat.gz, ~13 GB, and entry.list) are downloaded into raw/ and
kept for reference, but InterPro is not used in the mapping.
Rerun --prepare-domain-data after deleting domain_data_location to refresh the data.

Safety: every file runs in its own process, and each PDBFixer fragment in a short-lived
child process with a time limit (GAP_TIMEOUT). A file that exceeds FILE_TIMEOUT or
WORKER_MAX_MEMORY_GB, or whose process dies, is retried once without rebuilding. New files
start only while MIN_FREE_MEMORY_GB is available; below LOW_MEMORY_GB the largest worker
is stopped and its file rerun alone. Ctrl-C stops cleanly; --resume continues a run.

Usage: 1. unique chains:  python py_unique_chain_selector.py [mmcif_files_location]
              [-o unique_chain_location] [--workers N] [--resume | --overwrite] [--limit N] [--no-rebuild]
       2. domain data (once, any time before step 3):
              python py_unique_chain_selector.py --prepare-domain-data [--without-interpro]
              (~14 GB download; --without-interpro skips the 13 GB InterPro file, which is kept
              for reference only)
       3. domain mapping: python py_unique_chain_selector.py --map-domains [-o unique_chain_location]
Requires Python >= 3.12 with gemmi, numpy, scipy, tqdm, pdbfixer and openmm (see
environment.yml). Runs on Linux, on the CPU; no GPU needed.
"""
import os

# One thread per library: the script runs several processes instead, and a process that
# holds no thread pools can be forked safely.
for _var in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENMM_CPU_THREADS'):
    os.environ.setdefault(_var, '1')

import warnings

warnings.filterwarnings('ignore')   # library warnings (OpenMM, numpy) would clutter the progress bar

import argparse
import csv
import gzip
import hashlib
import io
import multiprocessing as mp
import multiprocessing.connection
import re
import shutil
import signal
import sqlite3
import string
import sys
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict, deque

import gemmi
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm

import pdbfixer  # noqa: F401  imported once here so the fragment processes share it
from openmm import unit

# ---- inputs ---------------------------------------------------------------------------
hetatm_list = "MG, CA, ZN, NA, K, MN, FE, FE2, CO, NI, CU, CL, PO4, SF4, FES, HEM, HEC, FAD, FMN, NAD, NAP, NDP, COA, SAM, SAH, PLP, ATP, ADP, AMP, GTP, GDP, NAG, BMA, MAN, FUC, GAL, CLR, SEP, TPO, PTR, KCX, CSO, PCA, PSU, 5MC, OMG, OMC, 5MU, OMU, MA6, UR3, 2MG"  # 52 het_codes
mmcif_files_location = "mmcif_files/"                 # input folder (relative paths start at the
unique_chain_location = "unique_chains/"              #   current directory; change them here or
domain_data_location = "unique_chain_references/"     #   on the command line)

# ---- settings -------------------------------------------------------------------------
HBOND_CUTOFF = 3.5            # A, water O to an N/O(/S) partner
MIN_HBONDS = 3                # polymer H-bonds that make a water 'structural'
MAX_EXPOSED = 0.05            # largest open share of a 'buried' water's surface
METAL_CUTOFF = {'MG': 2.5, 'ZN': 2.6, 'FE': 2.6, 'MN': 2.6, 'CO': 2.6, 'NI': 2.6, 'CU': 2.6,
                'CA': 2.9, 'NA': 2.9, 'K': 3.3}
METAL_CUTOFF_DEFAULT = 3.0
MAX_GAP_LENGTH = 30           # longer internal gaps are not rebuilt
MAX_PEPTIDE_BOND = 1.8        # A; a rebuilt loop with a longer C-N bond is rejected
GAP_TIMEOUT = 180             # s for one PDBFixer fragment; after that the gap stays open
FILE_TIMEOUT = 3600           # s for one file; after that it is retried without rebuilding
WORKER_MAX_MEMORY_GB = 6.0    # a file above this is stopped and retried without rebuilding
MIN_FREE_MEMORY_GB = 2.5      # a new file starts only with this much memory available
LOW_MEMORY_GB = 1.0           # below this the largest worker is stopped and rerun alone
DOWNLOAD_TIMEOUT = 60         # s without data before a download is retried

# ---- constants ------------------------------------------------------------------------
HET_CODES = {code.strip().upper() for code in hetatm_list.split(',') if code.strip()}
PEPTIDE = {gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD}
HET_ENTITY = {gemmi.EntityType.NonPolymer, gemmi.EntityType.Branched}
WATER_ROLES = ('metal', 'ligand', 'buried', 'interface', 'structural')
VDW_RADIUS = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'P': 1.8, 'SE': 1.9}
VDW_DEFAULT, PROBE = 1.8, 1.4
_n = np.arange(96) + 0.5                       # 96 points spread evenly on a unit sphere
_phi, _theta = np.arccos(1 - 2 * _n / 96), np.pi * (1 + 5 ** 0.5) * _n
SPHERE = np.c_[np.cos(_theta) * np.sin(_phi), np.sin(_theta) * np.sin(_phi), np.cos(_phi)]
AA_HEAVY = {  # heavy atoms of the standard amino acids (OXT left out)
    'ALA': 'N CA C O CB', 'ARG': 'N CA C O CB CG CD NE CZ NH1 NH2', 'ASN': 'N CA C O CB CG OD1 ND2',
    'ASP': 'N CA C O CB CG OD1 OD2', 'CYS': 'N CA C O CB SG', 'GLN': 'N CA C O CB CG CD OE1 NE2',
    'GLU': 'N CA C O CB CG CD OE1 OE2', 'GLY': 'N CA C O', 'HIS': 'N CA C O CB CG ND1 CD2 CE1 NE2',
    'ILE': 'N CA C O CB CG1 CG2 CD1', 'LEU': 'N CA C O CB CG CD1 CD2', 'LYS': 'N CA C O CB CG CD CE NZ',
    'MET': 'N CA C O CB CG SD CE', 'PHE': 'N CA C O CB CG CD1 CD2 CE1 CE2 CZ', 'PRO': 'N CA C O CB CG CD',
    'SER': 'N CA C O CB OG', 'THR': 'N CA C O CB OG1 CG2',
    'TRP': 'N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2', 'TYR': 'N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH',
    'VAL': 'N CA C O CB CG1 CG2'}
AA_HEAVY = {k: set(v.split()) for k, v in AA_HEAVY.items()}
THREE = {gemmi.find_tabulated_residue(n).one_letter_code.upper(): n for n in AA_HEAVY}
FRAGMENT_CHAIN_IDS = string.ascii_uppercase + string.ascii_lowercase + string.digits
ANCHOR = 4                    # residues kept on each side of a gap in its fragment
EXTENSIONS = ('.cif.gz', '.cif', '.mmcif.gz', '.mmcif')
_SIFTS = 'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/'
BULK_FILES = {                # local name -> URL (downloaded into domain_data_location/raw)
    'pdb_chain_uniprot.tsv.gz': _SIFTS + 'pdb_chain_uniprot.tsv.gz',
    'pdb_pfam_mapping.tsv.gz': _SIFTS + 'pdb_pfam_mapping.tsv.gz',
    'pdb_chain_scop2_uniprot.tsv.gz': _SIFTS + 'pdb_chain_scop2_uniprot.tsv.gz',
    'pdb_chain_scop2b_sf_uniprot.tsv.gz': _SIFTS + 'pdb_chain_scop2b_sf_uniprot.tsv.gz',
    'cath-b-newest-all.gz': 'https://download.cathdb.info/cath/releases/daily-release/newest/cath-b-newest-all.gz',
    'cath-names.txt': 'https://download.cathdb.info/cath/releases/latest-release/cath-classification-data/'
                      'cath-names.txt',
    'Pfam-A.clans.tsv.gz': 'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz',
    'scop-cla-latest.txt': 'https://www.ebi.ac.uk/pdbe/scop/files/scop-cla-latest.txt',
    'scop-des-latest.txt': 'https://www.ebi.ac.uk/pdbe/scop/files/scop-des-latest.txt',
    'entry.list': 'https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list',
    'protein2ipr.dat.gz': 'https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz',
}
UNIPROT_FIELDS = 'accession,id,protein_name,organism_name'
UNIPROT_STREAM = ('https://rest.uniprot.org/uniprotkb/stream?query=structure_3d:true'
                  f'&fields={UNIPROT_FIELDS}&format=tsv&compressed=true')
UNIPROT_SEARCH = 'https://rest.uniprot.org/uniprotkb/search?query={}&fields=' + UNIPROT_FIELDS + '&format=tsv&size=500'
DOMAIN_DB = 'domain_annotations.sqlite'
DOMAIN_SOURCES = ('Pfam', 'CATH', 'SCOP2', 'UniProt')      # written in this order; InterPro is not used
SUMMARY_COLUMNS = ['pdbid_chainid', 'pdb_id', 'chain', 'chain_type', 'status', 'unique_chain', 'group_members',
                   'n_members', 'sequence_length', 'het_count', 'hets', 'water_metal', 'water_ligand',
                   'water_buried', 'water_interface', 'water_structural', 'water_significant', 'water_total',
                   'residues_rebuilt', 'atoms_added', 'domains', 'notes', 'source_file', 'sequence']
DOMAIN_COLUMNS = ['database', 'accession', 'name', 'description', 'residues_new_numbering', 'n_residues',
                  'seqres_range', 'author_range', 'details']
RESIDUE_MAP_COLUMNS = ['new_number', 'residue', 'seqres_number', 'author_number', 'status']
RANGE = re.compile(r'^(-?\d+[A-Za-z]?)-(-?\d+[A-Za-z]?)$')


# ---- reading and cleaning a structure ---------------------------------------------------

def entry_id(path):
    name = os.path.basename(path)
    for ext in EXTENSIONS:
        if name.lower().endswith(ext):
            return name[:-len(ext)].lower()
    return name.lower()


def read_entry(path):
    """First model without hydrogens, the entity descriptions, and per author chain the
    author residue numbers in SEQRES order as [(author number + ins. code, label_seq)]."""
    block = gemmi.cif.read(path)[0]
    st = gemmi.make_structure_from_block(block)
    if len(st) == 0 or not any(len(ch) for ch in st[0]):
        raise ValueError('no atom coordinates in the file')
    while len(st) > 1:
        del st[len(st) - 1]
    st.remove_hydrogens()
    st.setup_entities()
    st.merge_chain_parts()        # one Chain per author chain: polymer, then its ligands and waters
    ids = [gemmi.cif.as_string(v) for v in block.find_values('_entity.id')]
    names = [gemmi.cif.as_string(v) for v in block.find_values('_entity.pdbx_description')]
    scheme = defaultdict(list)
    for row in block.find('_pdbx_poly_seq_scheme.', ['pdb_strand_id', 'seq_id', 'pdb_seq_num', 'pdb_ins_code']):
        icode = row.str(3) if row.str(3) not in ('.', '?') else ''
        if row.str(1).isdigit():
            scheme[row.str(0)].append((f'{row.str(2)}{icode}', int(row.str(1))))
    return st, dict(zip(ids, names)), scheme


def resolve_altlocs(model):
    """Keep one conformer everywhere, chosen by occupancy (ties: first in the file)."""
    for chain in model:
        # Different residue types at one position (microheterogeneity)
        drop, i = [], 0
        while i < len(chain):
            j = i + 1
            while j < len(chain) and chain[j].seqid == chain[i].seqid:
                j += 1
            if j - i > 1:
                best = max(range(i, j), key=lambda k: sum(a.occ for a in chain[k]))
                drop += [k for k in range(i, j) if k != best]
            i = j
        for k in reversed(drop):
            del chain[k]
        # Whole-residue alternates: residues whose atoms all sit on one altloc and that
        # overlap a residue on another altloc (e.g. AGS 'A' / ADP 'B' at one site)
        alt_only = []
        for k, res in enumerate(chain):
            locs = {a.altloc for a in res}
            if len(locs) == 1 and '\x00' not in locs:
                alt_only.append((k, locs.pop(), float(np.mean([a.occ for a in res]))))
        if len(alt_only) > 1:
            xyz, owner = [], []
            for n, (k, _, _) in enumerate(alt_only):
                for a in chain[k]:
                    xyz.append(a.pos.tolist())
                    owner.append(n)
            parent = list(range(len(alt_only)))

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x
            for p, q in cKDTree(np.array(xyz)).query_pairs(1.5):
                a, b = owner[p], owner[q]
                if alt_only[a][1] != alt_only[b][1]:
                    parent[find(a)] = find(b)
            clusters = defaultdict(list)
            for n in range(len(alt_only)):
                clusters[find(n)].append(n)
            drop = []
            for members in clusters.values():
                occ = defaultdict(list)
                for n in members:
                    occ[alt_only[n][1]].append(alt_only[n][2])
                if len(occ) > 1:
                    best = max(occ, key=lambda loc: np.mean(occ[loc]))
                    drop += [alt_only[n][0] for n in members if alt_only[n][1] != best
                             and chain[alt_only[n][0]].entity_type != gemmi.EntityType.Polymer]
            for k in sorted(drop, reverse=True):
                del chain[k]
        # Alternate conformers inside a residue
        for res in chain:
            occ = {}
            for a in res:
                if a.altloc != '\x00':
                    occ.setdefault(a.altloc, []).append(a.occ)
            if not occ:
                continue
            best = max(occ, key=lambda loc: np.mean(occ[loc]))
            for k in reversed(range(len(res))):
                if res[k].altloc not in ('\x00', best):
                    del res[k]
            for a in res:
                a.altloc = '\x00'


def classify_chains(st, descriptions):
    entities = {e.name: e for e in st.entities}
    chains = []
    for ch in st[0]:
        polymer = [r for r in ch if r.entity_type == gemmi.EntityType.Polymer]
        entity = entities.get(polymer[0].entity_id) if polymer else None
        if not polymer:
            kind = 'none'
        elif entity is not None and entity.polymer_type in PEPTIDE:
            kind = 'protein'
        elif entity is not None and entity.polymer_type in (gemmi.PolymerType.Dna, gemmi.PolymerType.Rna,
                                                            gemmi.PolymerType.DnaRnaHybrid):
            kind = 'nucleic acid'
        else:
            kind = 'other polymer'
        chains.append({'order': len(chains), 'name': ch.name, 'kind': kind, 'polymer': polymer, 'entity': entity,
                       'description': descriptions.get(polymer[0].entity_id, '') if polymer else '',
                       'hets': [r for r in ch if r.entity_type in HET_ENTITY and r.name in HET_CODES],
                       'waters': [r for r in ch if r.entity_type == gemmi.EntityType.Water]})
    return chains


# ---- water roles ------------------------------------------------------------------------

def exposed_fraction(center, tree, heavy, radius):
    near = tree.query_ball_point(center, VDW_RADIUS['O'] + 2 * PROBE + max(VDW_RADIUS.values()))
    if not near:
        return 1.0
    points = center + (VDW_RADIUS['O'] + PROBE) * SPHERE
    dist = np.linalg.norm(points[:, None, :] - heavy[near][None, :, :], axis=2)
    return float((~(dist < radius[near] + PROBE).any(axis=1)).mean())


def assign_water_roles(st, chains):
    """Set chain['roles'] = one role (or None) per water of the chain."""
    water_pos, water_ref, water_index = [], [], {}
    part_pos, part_kind, part_res, part_chain, part_name, part_elem = [], [], [], [], [], []
    heavy, radius = [], []
    uid = 0
    for ci, (ch, info) in enumerate(zip(st[0], chains)):
        info['roles'] = [None] * len(info['waters'])
        wi = 0
        for res in ch:
            uid += 1
            if res.entity_type == gemmi.EntityType.Water:
                oxygen = next((a for a in res if a.element.name == 'O'), None)
                if oxygen is not None:
                    water_pos.append(oxygen.pos.tolist())
                    water_ref.append((ci, wi))
                water_index[(ch.name, str(res.seqid))] = (ci, wi)
                wi += 1
                continue
            polymer = res.entity_type == gemmi.EntityType.Polymer
            listed = res.entity_type in HET_ENTITY and res.name in HET_CODES
            for a in res:
                elem = a.element.name.upper()
                heavy.append(a.pos.tolist())
                radius.append(VDW_RADIUS.get(elem, VDW_DEFAULT))
                if polymer and elem in ('N', 'O'):
                    kind = 'polymer'
                elif listed and a.element.is_metal:
                    kind = 'metal'
                elif listed and elem in ('N', 'O', 'S'):
                    kind = 'het'
                else:
                    continue
                part_pos.append(a.pos.tolist())
                part_kind.append(kind)
                part_res.append(uid)
                part_chain.append(ch.name)
                part_name.append(a.name)
                part_elem.append(elem)
    if not water_pos:
        return
    linked = set()                               # deposited water-metal links, no symmetry
    for con in st.connections:
        if con.type != gemmi.ConnectionType.MetalC or con.asu != gemmi.Asu.Same:
            continue
        for w, m in ((con.partner1, con.partner2), (con.partner2, con.partner1)):
            key = (w.chain_name, str(w.res_id.seqid))
            if key in water_index and m.res_id.name in HET_CODES:
                linked.add(water_index[key])
    water_pos = np.array(water_pos)
    part_xyz = np.array(part_pos).reshape(-1, 3)
    part_tree = cKDTree(part_xyz) if part_pos else None
    heavy, radius = np.array(heavy).reshape(-1, 3), np.array(radius)
    heavy_tree = cKDTree(heavy) if len(heavy) else None
    reach = max(HBOND_CUTOFF, METAL_CUTOFF_DEFAULT, *METAL_CUTOFF.values())
    neighbours = part_tree.query_ball_point(water_pos, reach) if part_tree is not None else [[]] * len(water_pos)
    for (ci, wi), center, near in zip(water_ref, water_pos, neighbours):
        role = 'metal' if (ci, wi) in linked else None
        polar, residues, poly_chains, ligands = set(), set(), set(), set()
        near = list(near) if role is None else []
        for p, d in zip(near, np.linalg.norm(part_xyz[near] - center, axis=1) if near else []):
            if part_kind[p] == 'metal':
                if d <= METAL_CUTOFF.get(part_elem[p], METAL_CUTOFF_DEFAULT):
                    role = 'metal'
                    break
            elif d <= HBOND_CUTOFF:
                if part_kind[p] == 'polymer':
                    polar.add((part_res[p], part_name[p]))
                    residues.add(part_res[p])
                    poly_chains.add(part_chain[p])
                else:
                    ligands.add(part_res[p])
        if role is None:
            if ligands and (polar or len(ligands) > 1):
                role = 'ligand'
            elif (polar or ligands) and heavy_tree is not None \
                    and exposed_fraction(center, heavy_tree, heavy, radius) <= MAX_EXPOSED:
                role = 'buried'
            elif len(poly_chains) > 1:
                role = 'interface'
            elif len(polar) >= MIN_HBONDS and len(residues) > 1:
                role = 'structural'
        chains[ci]['roles'][wi] = role


# ---- rebuilding with PDBFixer -------------------------------------------------------------

def _isolated_child(func, args, conn):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        conn.send((True, func(*args)))
    except BaseException as exc:  # report every failure to the parent
        conn.send((False, f'{type(exc).__name__}: {exc}'))
    finally:
        conn.close()


def run_isolated(func, args, timeout):
    """Run func(*args) in a forked child; (True, result) or (False, reason)."""
    ctx = mp.get_context('fork')
    recv, send = ctx.Pipe(duplex=False)
    proc = ctx.Process(target=_isolated_child, args=(func, args, send), daemon=True)
    proc.start()
    send.close()
    try:
        if recv.poll(timeout):
            return recv.recv()
        return False, f'time limit {timeout} s'
    except (EOFError, OSError):
        return False, 'fragment process ended unexpectedly'
    finally:
        if proc.is_alive():
            proc.kill()
        proc.join()
        recv.close()


def fixer_job(pdb_text, insert, targets, seed=1):
    """PDBFixer on one fragment: insert = (position in chain 0, residue names) or None;
    missing atoms are added to the chains in targets only. Returns, per target chain, its
    residues as (name, [(atom, element, x, y, z), ...]) in Angstrom."""
    from pdbfixer import PDBFixer
    fixer = PDBFixer(pdbfile=io.StringIO(pdb_text))
    fixer.missingResidues = {(0, insert[0]): list(insert[1])} if insert else {}
    fixer.findMissingAtoms()
    fixer.missingAtoms = {r: a for r, a in fixer.missingAtoms.items() if r.chain.index in targets}
    fixer.missingTerminals = {}
    fixer.addMissingAtoms(seed=seed)
    pos = fixer.positions.value_in_unit(unit.angstrom)
    chains = list(fixer.topology.chains())
    return [[(res.name, [(a.name, a.element.symbol if a.element else '', *pos[a.index]) for a in res.atoms()])
             for res in chains[t].residues()] for t in targets]


def fragment_pdb(groups):
    """PDB text with each residue list in its own chain. Modified residues are written as
    their standard parent (SEP as SER, MSE as MET): PDBFixer has no template for them, and
    one next to a gap stops the rebuilt loop from joining that side. Residues are numbered
    1, 2, 3... in each chain: rebuilt residues have no number yet, and PDBFixer's reader
    treats consecutive residues with the same number as one."""
    st = gemmi.Structure()
    model = gemmi.Model(1)
    for cid, residues in zip(FRAGMENT_CHAIN_IDS, groups):
        ch = gemmi.Chain(cid)
        for k, r in enumerate(residues, 1):
            parent = parent_name(r.name) if r.name not in AA_HEAVY else None
            r = r.clone()
            r.seqid = gemmi.SeqId(k, ' ')
            if parent:
                r.name, r.het_flag = parent, 'A'
                for i in reversed(range(len(r))):
                    if r[i].name not in AA_HEAVY[parent]:
                        del r[i]
            ch.add_residue(r)
        model.add_chain(ch)
    st.add_model(model)
    return st.make_pdb_string(gemmi.PdbWriteOptions(minimal=True))


def contiguous(residues, i, j):
    a, b = residues[i].label_seq, residues[j].label_seq
    return j == i + 1 and (a is None or b is None or b == a + 1)


def environment(residues, exclude, distances, limit):
    """Contiguous runs of residues not in exclude, given each residue's distance to the
    site (np.inf = not near), nearest runs first, at most limit runs."""
    runs, cur = [], []
    for i in range(len(residues)):
        if i in exclude or not np.isfinite(distances[i]):
            continue
        if cur and not contiguous(residues, cur[-1], i):
            runs.append(cur)
            cur = []
        cur.append(i)
    if cur:
        runs.append(cur)
    runs.sort(key=lambda run: min(distances[i] for i in run))
    return runs[:limit]


def coords(res):
    return np.array([a.pos.tolist() for a in res]).reshape(-1, 3)


def new_atom(name, elem, x, y, z, b_iso):
    a = gemmi.Atom()
    a.name, a.element, a.pos, a.occ, a.b_iso = name, gemmi.Element(elem or name[:1]), gemmi.Position(x, y, z), 0.0, b_iso
    return a


def add_missing_atoms(res, built, b_iso):
    present = {a.name for a in res}
    added = 0
    for name, elem, x, y, z in built:
        if name not in present:
            res.add_atom(new_atom(name, elem, x, y, z, b_iso))
            added += 1
    return added


def longest_peptide_bond(residues, first, last):
    """Longest C(k)-N(k+1) distance for k = first..last-1 in a fragment result."""
    def atom(res, name):
        return next((np.array(a[2:]) for a in res[1] if a[0] == name), None)
    worst = 0.0
    for k in range(first, last):
        c, n = atom(residues[k], 'C'), atom(residues[k + 1], 'N')
        worst = max(worst, float(np.linalg.norm(c - n)) if c is not None and n is not None else np.inf)
    return worst


def parent_name(name):
    if name in AA_HEAVY:
        return name
    info = gemmi.find_tabulated_residue(name)
    return THREE.get(info.one_letter_code.upper()) if info is not None else None


def rebuild(residues, seqres, notes):
    """Fill internal gaps and missing heavy atoms of standard residues. Returns the new
    residue list and (residues built, atoms added)."""
    if not residues:
        return residues, 0, 0
    backbone_ok = [{'N', 'CA', 'C'} <= {a.name for a in r} for r in residues]
    if sum(backbone_ok) < 0.5 * len(residues):
        notes.append('mostly backbone-incomplete (e.g. CA-only): not rebuilt')
        return residues, 0, 0
    b_iso = max((a.b_iso for r in residues for a in r), default=0.0)
    built_residues = atoms_added = 0
    inserts = {}
    labels = [r.label_seq for r in residues]
    gap_at = [] if None in labels else [i + 1 for i in range(len(residues) - 1) if labels[i + 1] - labels[i] > 1]
    if None in labels:
        notes.append('no label_seq_id: gaps not rebuilt')
    for g in gap_at:
        first, last = labels[g - 1] + 1, labels[g] - 1
        names = [parent_name(n) for n in seqres[first - 1:last]] if last <= len(seqres) else [None]
        if len(names) > MAX_GAP_LENGTH or None in names:
            why = f'longer than {MAX_GAP_LENGTH}' if len(names) > MAX_GAP_LENGTH else 'unknown residue type'
            notes.append(f'gap {first}-{last} ({len(names)} res) not rebuilt: {why}')
            continue
        lo = max([x for x in gap_at if x < g] + [0, g - ANCHOR])
        hi = min([x for x in gap_at if x > g] + [len(residues), g + ANCHOR])
        anchor = set(range(lo, hi))
        ends = [next((a for a in residues[g - 1] if a.name == 'C'), residues[g - 1][len(residues[g - 1]) - 1]),
                next((a for a in residues[g] if a.name == 'N'), residues[g][0])]
        a_pos, b_pos = (np.array(e.pos.tolist()) for e in ends)
        reach = min(8.0 + 0.25 * len(names), 20.0)
        seg = b_pos - a_pos
        dist = np.full(len(residues), np.inf)
        for i, r in enumerate(residues):
            if i not in anchor:
                p = coords(r)
                t = np.clip(((p - a_pos) @ seg) / max(seg @ seg, 1e-9), 0, 1)
                d = np.linalg.norm(p - (a_pos + t[:, None] * seg), axis=1).min()
                dist[i] = d if d <= reach else np.inf
        runs = environment(residues, anchor, dist, len(FRAGMENT_CHAIN_IDS) - 1)
        groups = [[residues[i] for i in range(lo, hi)]] + [[residues[i] for i in run] for run in runs]
        fragment, worst = fragment_pdb(groups), np.inf
        for seed in (1, 2, 3):                   # a loop that does not close gets new seeds
            ok, result = run_isolated(fixer_job, (fragment, (g - lo, names), [0], seed), GAP_TIMEOUT)
            if not ok or len(result[0]) != (hi - lo) + len(names):
                break
            worst = longest_peptide_bond(result[0], g - lo - 1, g - lo + len(names))
            if worst <= MAX_PEPTIDE_BOND:
                break
        if not ok:
            notes.append(f'gap {first}-{last} ({len(names)} res) not rebuilt: {result}')
            continue
        chain0 = result[0]
        if len(chain0) != (hi - lo) + len(names):
            notes.append(f'gap {first}-{last} not rebuilt: unexpected fragment result')
            continue
        if worst > MAX_PEPTIDE_BOND:
            notes.append(f'gap {first}-{last} ({len(names)} res) not rebuilt: loop did not close '
                         f'(C-N {worst:.1f} A after 3 tries)')
            continue
        new = []
        for k, (name, atoms) in enumerate(chain0[g - lo:g - lo + len(names)]):
            r = gemmi.Residue()
            r.name, r.het_flag, r.entity_type, r.label_seq = name, 'A', gemmi.EntityType.Polymer, first + k
            r.seqid = gemmi.SeqId(0, ' ')
            for atom in atoms:
                r.add_atom(new_atom(*atom, b_iso))
                atoms_added += 1
            new.append(r)
        flanks = chain0[:g - lo] + chain0[g - lo + len(names):]
        for i, (_, atoms) in zip(range(lo, hi), flanks):
            if residues[i].name in AA_HEAVY:     # modified residues went in as stand-ins
                atoms_added += add_missing_atoms(residues[i], atoms, b_iso)
        inserts[g] = new
        built_residues += len(new)
    if inserts:
        merged = []
        for i, r in enumerate(residues):
            merged += inserts.get(i, [])
            merged.append(r)
        residues = merged
    # Missing heavy atoms of standard residues, in clusters of nearby residues
    incomplete = [i for i, r in enumerate(residues)
                  if r.name in AA_HEAVY and AA_HEAVY[r.name] - {a.name for a in r}
                  and {'N', 'CA', 'C'} <= {a.name for a in r}]
    clusters = []
    for i in incomplete:
        ca = residues[i].find_atom('CA', '*').pos
        if clusters and len(clusters[-1]) < 8 and residues[clusters[-1][0]].find_atom('CA', '*').pos.dist(ca) <= 15:
            clusters[-1].append(i)
        else:
            clusters.append([i])
    for cluster in clusters:
        target_runs = []
        for i in cluster:
            if target_runs and contiguous(residues, target_runs[-1][-1], i):
                target_runs[-1].append(i)
            else:
                target_runs.append([i])
        tree = cKDTree(np.vstack([coords(residues[i]) for i in cluster]))
        dist = np.array([tree.query(coords(r))[0].min() if len(r) else np.inf for r in residues])
        dist[dist > 8.0] = np.inf
        runs = environment(residues, set(cluster), dist, len(FRAGMENT_CHAIN_IDS) - len(target_runs))
        groups = [[residues[i] for i in run] for run in target_runs] + [[residues[i] for i in run] for run in runs]
        ok, result = run_isolated(fixer_job, (fragment_pdb(groups), None, list(range(len(target_runs)))),
                                  GAP_TIMEOUT)
        if not ok:
            notes.append(f'missing atoms of {len(cluster)} residues not added: {result}')
            continue
        for run, built in zip(target_runs, result):
            for i, (_, atoms) in zip(run, built):
                atoms_added += add_missing_atoms(residues[i], atoms, b_iso)
    return residues, built_residues, atoms_added


# ---- domain and family data: bulk download and database (once) ------------------------------

def download(url, path):
    """Download url to path unless it is there; resumes a partial download, retries."""
    if os.path.exists(path):
        return path
    part, error = path + '.part', None
    for attempt in range(6):
        done = os.path.getsize(part) if os.path.exists(part) else 0
        headers = {'User-Agent': 'py_unique_chain_selector'}
        if done:
            headers['Range'] = f'bytes={done}-'
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=DOWNLOAD_TIMEOUT) as resp:
                if done and resp.status != 206:          # the server ignored the range: start again
                    done = 0
                size = resp.headers.get('Content-Length')
                total = int(size) + done if size else None
                with open(part, 'ab' if done else 'wb') as fh, \
                        tqdm(total=total, initial=done, unit='B', unit_scale=True, leave=False,
                             desc=os.path.basename(path)) as bar:
                    while chunk := resp.read(1 << 20):
                        fh.write(chunk)
                        bar.update(len(chunk))
            if total is None or os.path.getsize(part) == total:
                os.replace(part, path)
                return path
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            error = exc
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f'download failed after 6 tries: {url} ({error})')


def check_md5(path, md5_url):
    """Compare a file with the MD5 sum published next to it; delete it if they differ."""
    with urllib.request.urlopen(urllib.request.Request(md5_url, headers={'User-Agent': 'py_unique_chain_selector'}),
                                timeout=DOWNLOAD_TIMEOUT) as resp:
        expected = resp.read().decode().split()[0]
    digest = hashlib.md5()
    with open(path, 'rb') as fh:
        while chunk := fh.read(1 << 24):
            digest.update(chunk)
    if digest.hexdigest() != expected:
        os.remove(path)
        raise RuntimeError(f'{os.path.basename(path)} is damaged (MD5 mismatch); rerun to download it again')


def tsv_rows(path):
    """Rows of a SIFTS table as dicts (the first line holds the release, the second the header)."""
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        line = fh.readline()
        if line.startswith('#'):
            line = fh.readline()
        columns = line.rstrip('\n').split('\t')
        for line in fh:
            yield dict(zip(columns, line.rstrip('\n').split('\t')))


def download_uniprot_names(path, accessions):
    """Entry name, protein name and organism of the UniProt entries used by SIFTS: one
    stream of all entries with a 3D structure, then small batched queries for the rest."""
    if os.path.exists(path):
        return
    stream = download(UNIPROT_STREAM, path + '.structures.gz')
    with gzip.open(stream, 'rt', encoding='utf-8') as fh:
        lines = fh.read().splitlines()
    found = {line.split('\t', 1)[0] for line in lines[1:]}
    missing = sorted(accessions - found)
    for i in tqdm(range(0, len(missing), 100), desc='UniProt names not in the stream', leave=False):
        query = urllib.parse.quote(' OR '.join(f'accession:{a}' for a in missing[i:i + 100]))
        request = urllib.request.Request(UNIPROT_SEARCH.format(query), headers={'User-Agent': 'py_unique_chain_selector'})
        try:
            with urllib.request.urlopen(request, timeout=DOWNLOAD_TIMEOUT) as resp:
                lines += resp.read().decode('utf-8').splitlines()[1:]
        except (urllib.error.URLError, OSError, TimeoutError):
            pass                                          # names stay empty for these few entries
        time.sleep(0.5)
    with gzip.open(path + '.part', 'wt', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')
    os.replace(path + '.part', path)
    os.remove(stream)


def build_database(raw, db_path):
    """One table, one row per chain and annotation; ranges in SEQRES numbering (label_seq)
    except CATH, whose boundaries use author numbering."""
    def p(name):
        return os.path.join(raw, name)
    tmp = db_path + '.building'
    if os.path.exists(tmp):
        os.remove(tmp)
    con = sqlite3.connect(tmp)
    con.execute('PRAGMA journal_mode = OFF')
    con.execute('PRAGMA synchronous = OFF')
    con.execute('CREATE TABLE domains (pdb_id TEXT, chain TEXT, db TEXT, accession TEXT, name TEXT, '
                'description TEXT, numbering TEXT, ranges TEXT, author_range TEXT, details TEXT)')
    con.execute('CREATE TABLE info (item TEXT, value TEXT)')
    insert = 'INSERT INTO domains VALUES (?,?,?,?,?,?,?,?,?,?)'
    rows, counts = [], Counter()

    def flush(force=False):
        if rows and (force or len(rows) >= 200000):
            con.executemany(insert, rows)
            counts.update(r[2] for r in rows)
            rows.clear()

    def as_int(text):
        try:
            return int(text)
        except (TypeError, ValueError):
            return None

    # UniProt
    names = {}
    with gzip.open(p('uniprot_names.tsv.gz'), 'rt', encoding='utf-8') as fh:
        for line in fh:
            f = (line.rstrip('\n').split('\t') + [''] * 4)[:4]
            if f[0] and f[0] != 'Entry':
                names[f[0]] = (f[1], f'{f[2]} ({f[3]})' if f[3] else f[2])
    for r in tqdm(tsv_rows(p('pdb_chain_uniprot.tsv.gz')), desc='UniProt', unit=' rows', leave=False):
        rb, re_, sb, se = (as_int(r.get(k)) for k in ('RES_BEG', 'RES_END', 'SP_BEG', 'SP_END'))
        if None in (rb, re_, sb, se):
            continue
        pdb, acc = r['PDB'].lower(), r['SP_PRIMARY']
        entry, desc = names.get(acc, ('', ''))
        auth = f'{r["PDB_BEG"]}-{r["PDB_END"]}' if r.get('PDB_BEG') and r.get('PDB_END') else ''
        rows.append((pdb, r['CHAIN'], 'UniProt', acc, entry, desc, 'seqres', f'{rb}-{re_}', auth,
                     f'UniProt residues {sb}-{se}'))
        flush()
    # Pfam
    clans = {}
    with gzip.open(p('Pfam-A.clans.tsv.gz'), 'rt', encoding='utf-8') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 5:
                clans[f[0]] = f[4] + (f' (clan {f[2]})' if f[2] else '')
    for r in tqdm(tsv_rows(p('pdb_pfam_mapping.tsv.gz')), desc='Pfam', unit=' rows', leave=False):
        acc = r['PFAM_ACCESSION'].split('.')[0]
        auth = (f'{r["AUTH_PDBRES_START"]}{r["AUTH_PDBRES_START_INS_CODE"].strip()}-'
                f'{r["AUTH_PDBRES_END"]}{r["AUTH_PDBRES_END_INS_CODE"].strip()}')
        rows.append((r['PDB'].lower(), r['CHAIN'], 'Pfam', acc, r['PFAM_NAME'], clans.get(acc, ''), 'seqres',
                     f'{r["PDB_START"]}-{r["PDB_END"]}', auth,
                     f'UniProt {r["UNIPROT_ACCESSION"]} residues {r["UNP_START"]}-{r["UNP_END"]}'))
        flush()
    # CATH (boundaries in author numbering, e.g. "0-153:A" or "10-60:A,70A-90:A")
    cath_names = {}
    with open(p('cath-names.txt'), encoding='utf-8', errors='replace') as fh:
        for line in fh:
            parts = line.split(None, 2)
            if not line.startswith('#') and len(parts) == 3:
                cath_names[parts[0]] = parts[2].strip().lstrip(':')
    segment = re.compile(r'^(-?\d+[A-Za-z]?)-(-?\d+[A-Za-z]?):(\S+)$')
    with gzip.open(p('cath-b-newest-all.gz'), 'rt', encoding='utf-8') as fh:
        for line in tqdm(fh, desc='CATH', unit=' domains', leave=False):
            f = line.split()
            if len(f) < 4:
                continue
            domain, code = f[0], f[2]
            per_chain = defaultdict(list)
            for seg in f[3].split(','):
                m = segment.match(seg)
                if m:
                    per_chain[m.group(3)].append(f'{m.group(1)}-{m.group(2)}')
            levels = code.split('.')
            lineage = ' / '.join(cath_names.get('.'.join(levels[:k]), '') for k in range(1, len(levels)))
            for chain, segs in per_chain.items():
                rows.append((domain[:4].lower(), chain, 'CATH', code, cath_names.get(code, ''), lineage, 'author',
                             ','.join(segs), ','.join(segs), f'domain {domain}'))
            flush()
    # SCOP2
    des = {}
    with open(p('scop-des-latest.txt'), encoding='utf-8') as fh:
        for line in fh:
            if not line.startswith('#') and line.strip():
                node, _, name = line.strip().partition(' ')
                des[node] = name
    sf_info, fa_info = {}, {}
    with open(p('scop-cla-latest.txt'), encoding='utf-8') as fh:
        for line in fh:
            f = line.split()
            if line.startswith('#') or len(f) < 11:
                continue
            lv = dict(x.split('=') for x in f[10].split(','))
            lineage = f'class: {des.get(lv.get("CL"), "")}; fold: {des.get(lv.get("CF"), "")}'
            sf_info[f[5]] = (des.get(lv.get('SF'), ''), f'superfamily {lv.get("SF")}; {lineage}')
            fa_info[f[0]] = (des.get(lv.get('FA'), ''),
                             f'family {lv.get("FA")}; superfamily: {des.get(lv.get("SF"), "")}; {lineage}')
    for r in tsv_rows(p('pdb_chain_scop2b_sf_uniprot.tsv.gz')):
        name, desc = sf_info.get(r['SF_DOMID'], ('', ''))
        rows.append((r['PDB'].lower(), r['CHAIN'], 'SCOP2', f'SF-DOMID:{r["SF_DOMID"]}', name, desc, 'seqres',
                     f'{r["RES_BEG"]}-{r["RES_END"]}', f'{r["PDB_BEG"]}-{r["PDB_END"]}', 'superfamily domain'))
        flush()
    for r in tsv_rows(p('pdb_chain_scop2_uniprot.tsv.gz')):
        name, desc = fa_info.get(r['FA_DOMID'], ('', ''))
        rows.append((r['PDB'].lower(), r['CHAIN'], 'SCOP2', f'FA-DOMID:{r["FA_DOMID"]}', name, desc, 'seqres', '',
                     '', 'family assigned to the chain (the bulk data give no residue range)'))
        flush()
    flush(force=True)
    con.execute('CREATE INDEX domains_pdb ON domains (pdb_id)')
    for name in ('pdb_chain_uniprot.tsv.gz', 'pdb_pfam_mapping.tsv.gz'):
        with gzip.open(p(name), 'rt', encoding='utf-8') as fh:
            con.execute('INSERT INTO info VALUES (?, ?)', (f'release ({name})', fh.readline().strip('# \n')))
    con.execute('INSERT INTO info VALUES (?, ?)', ('built', time.strftime('%Y-%m-%d %H:%M')))
    for db, n in counts.items():
        con.execute('INSERT INTO info VALUES (?, ?)', (f'rows {db}', str(n)))
    con.commit()
    con.close()
    os.replace(tmp, db_path)
    return counts


def prepare_domain_data(data_dir, with_interpro=True):
    """Download the bulk files (resumable) and build the domain database (Pfam, CATH, SCOP2,
    UniProt). The InterPro files are downloaded and kept in raw/ but not used; with
    with_interpro=False the 13 GB protein2ipr file is not downloaded."""
    raw = os.path.join(data_dir, 'raw')
    os.makedirs(raw, exist_ok=True)
    for name, url in BULK_FILES.items():
        if name == 'protein2ipr.dat.gz' and not with_interpro:
            continue
        print(f'download {name}', file=sys.stderr)
        fresh = not os.path.exists(os.path.join(raw, name))
        download(url, os.path.join(raw, name))
        if name == 'protein2ipr.dat.gz' and fresh:
            print('check MD5 of protein2ipr.dat.gz', file=sys.stderr)
            check_md5(os.path.join(raw, name), url + '.md5')
    accessions = {r['SP_PRIMARY'] for r in tsv_rows(os.path.join(raw, 'pdb_chain_uniprot.tsv.gz'))}
    print(f'UniProt names for {len(accessions)} entries', file=sys.stderr)
    download_uniprot_names(os.path.join(raw, 'uniprot_names.tsv.gz'), accessions)
    print('build the domain database', file=sys.stderr)
    counts = build_database(raw, os.path.join(data_dir, DOMAIN_DB))
    print('domain database ready: ' + ', '.join(f'{db} {n:,}' for db, n in sorted(counts.items())), file=sys.stderr)


# ---- domain and family data: lookup for one entry ---------------------------------------

def seqres_numbers(text):
    nums = set()
    for part in filter(None, text.split(',')):
        a, _, b = part.partition('-')
        if a.isdigit() and b.isdigit():
            nums.update(range(int(a), int(b) + 1))
    return nums


def author_to_seqres(text, scheme):
    """SEQRES positions of author-numbered segments such as '0-153' or '10-60,70A-90'."""
    position = {}
    for k, (label, _) in enumerate(scheme):
        position.setdefault(label, k)
    nums = set()
    for part in filter(None, text.split(',')):
        m = RANGE.match(part)
        if m and m.group(1) in position and m.group(2) in position:
            nums.update(seq for _, seq in scheme[position[m.group(1)]:position[m.group(2)] + 1])
    return nums


def entry_domain_rows(con, pdb_id):
    """{author chain: [database row, ...]} for one entry."""
    rows = defaultdict(list)
    marks = ','.join('?' * len(DOMAIN_SOURCES))
    for row in con.execute('SELECT chain, db, accession, name, description, numbering, ranges, author_range, '
                           f'details FROM domains WHERE pdb_id = ? AND db IN ({marks})', (pdb_id, *DOMAIN_SOURCES)):
        rows[row[0]].append(row)
    return rows


def chain_records(rows, scheme):
    """Domain records of one chain with SEQRES positions; scheme = [(author number, SEQRES
    number)] in SEQRES order, used for the author-numbered CATH boundaries."""
    records = []
    for _, db, acc, name, desc, numbering, text, auth, details in rows:
        if numbering == 'author':
            seqres = author_to_seqres(text, scheme)
            if text and not seqres:
                details = f'{details}; boundaries not found in this chain'
        else:
            seqres = seqres_numbers(text)
        records.append({'database': db, 'accession': acc, 'name': name, 'description': desc,
                        'seqres': seqres, 'author_range': '' if auth == '-' else auth, 'details': details})
    return records


def domain_summary(rows):
    """'Pfam:PF00042; CATH:1.10.490.10; ...' for the summary table."""
    accs = defaultdict(list)
    for row in rows:
        if row[2] not in accs[row[1]]:
            accs[row[1]].append(row[2])
    return '; '.join(f'{db}:{",".join(a)}' for db, a in accs.items())


def read_residue_map(path):
    """(scheme [(author number, SEQRES number)], {SEQRES number: new number}) of a chain."""
    scheme, label_to_new = [], {}
    with open(path, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row['seqres_number']:
                seq = int(row['seqres_number'])
                scheme.append((row['author_number'], seq))
                if row['new_number']:
                    label_to_new[seq] = int(row['new_number'])
    return scheme, label_to_new


def map_domains(out_root, db_path):
    """Write <stem>_domains.csv for every chain folder (from its residue map) and fill the
    domains column of the summary for every protein chain."""
    summary_path = os.path.join(out_root, 'unique_chain_summary.csv')
    folders = defaultdict(list)                       # pdb id -> [(folder, author chain)]
    no_map = 0
    for name in sorted(os.listdir(out_root)):
        if name.startswith(('_', '.')) or not os.path.isdir(os.path.join(out_root, name)):
            continue
        if os.path.exists(os.path.join(out_root, name, f'{name}_residue_map.csv')):
            pdb_id, _, chain = name.partition('_')
            folders[pdb_id].append((name, chain))
        else:
            no_map += 1
    summary = []
    if os.path.exists(summary_path):
        with open(summary_path, newline='', encoding='utf-8') as fh:
            summary = list(csv.DictReader(fh))
    wanted = defaultdict(set)                         # pdb id -> protein chains in the summary
    for row in summary:
        if row.get('chain_type') == 'protein':
            wanted[row['pdb_id']].add(row['chain'])
    con = sqlite3.connect(f'file:{db_path}?mode=ro', uri=True)
    found, written = {}, 0
    try:
        for pdb_id in tqdm(sorted(set(folders) | set(wanted)), desc='domain mapping', unit=' entries'):
            rows = entry_domain_rows(con, pdb_id)
            for chain in wanted.get(pdb_id, ()):
                found[(pdb_id, chain)] = domain_summary(rows.get(chain, []))
            for stem, chain in folders.get(pdb_id, ()):
                scheme, label_to_new = read_residue_map(os.path.join(out_root, stem, f'{stem}_residue_map.csv'))
                write_domains(os.path.join(out_root, stem, f'{stem}_domains.csv'),
                              chain_records(rows.get(chain, []), scheme), label_to_new, pdb_id, chain)
                written += 1
    finally:
        con.close()
    if summary:
        for row in summary:
            if (row['pdb_id'], row['chain']) in found:
                row['domains'] = found[(row['pdb_id'], row['chain'])]
        with open(summary_path + '.part', 'w', newline='', encoding='utf-8') as fh:
            writer = csv.DictWriter(fh, SUMMARY_COLUMNS, lineterminator='\n', extrasaction='ignore')
            writer.writeheader()
            writer.writerows(summary)
        os.replace(summary_path + '.part', summary_path)
    print(f'domain files written for {written} chains; summary updated for {len(found)} protein chains'
          + (f'; {no_map} folders without a residue map (rerun the processing for them)' if no_map else ''),
          file=sys.stderr)


def ranges(nums):
    out = []
    for n in sorted(nums):
        if out and n == out[-1][1] + 1:
            out[-1][1] = n
        else:
            out.append([n, n])
    return ','.join(f'{a}-{b}' if a != b else str(a) for a, b in out)


def write_domains(path, records, label_to_new, pdb_id, chain):
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, DOMAIN_COLUMNS, lineterminator='\n')
        writer.writeheader()
        if not records:
            writer.writerow({'database': 'none', 'details': f'no Pfam/CATH/SCOP2/UniProt assignment for '
                                                            f'{pdb_id} chain {chain} in the domain database'})
            return
        merged = {}                    # one row per database/accession/detail
        for rec in records:
            key = (rec['database'], rec['accession'], rec['details'])
            if key in merged:
                merged[key]['seqres'] = merged[key]['seqres'] | rec['seqres']
                merged[key]['author_range'] = ','.join(x for x in (merged[key]['author_range'],
                                                                    rec['author_range']) if x)
            else:
                merged[key] = dict(rec)
        order = {db: k for k, db in enumerate(DOMAIN_SOURCES)}
        for rec in sorted(merged.values(),
                          key=lambda r: (order.get(r['database'], 9), min(r['seqres'], default=0), r['accession'])):
            new = [label_to_new[n] for n in rec['seqres'] if n in label_to_new]
            writer.writerow({'database': rec['database'], 'accession': rec['accession'], 'name': rec['name'],
                             'description': rec['description'], 'residues_new_numbering': ranges(new),
                             'n_residues': len(new), 'seqres_range': ranges(rec['seqres']),
                             'author_range': rec['author_range'],
                             'details': rec['details'] or ('' if new or not rec['seqres']
                                                           else 'not in the written model')})


# ---- writing a chain ------------------------------------------------------------------------

def write_residue_map(path, scheme, seqres, observed, polymer, hets, waters):
    """Old and new numbering of every residue of a written chain: all SEQRES positions
    (observed, rebuilt or missing), then the HETATM groups and waters. Call it before
    write_chain, which renumbers the residues. waters = [(residue, role), ...]."""
    new_of = {r.label_seq: (i, r.name) for i, r in enumerate(polymer, 1) if r.label_seq is not None}
    rows, seen = [], set()
    for author, seq in scheme or [(str(r.seqid), r.label_seq) for r in polymer if r.label_seq is not None]:
        if seq in seen:
            continue
        seen.add(seq)
        new, name = new_of.get(seq, ('', seqres[seq - 1] if 0 < seq <= len(seqres) else ''))
        status = 'observed' if seq in observed else 'rebuilt' if new else 'missing'
        rows.append([new, name, seq, author if status != 'rebuilt' or scheme else '', status])
    n = len(polymer)
    for k, (r, status) in enumerate([(h, 'hetatm') for h in hets] + [(w, f'water ({role})') for w, role in waters]):
        rows.append([n + k + 1, r.name, '', str(r.seqid), status])
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh, lineterminator='\n')
        writer.writerow(RESIDUE_MAP_COLUMNS)
        writer.writerows(rows)


def write_chain(folder, stem, chain_name, title, method, polymer, hets, waters, polymer_type):
    """Renumber, centre and write one chain as mmCIF and PDB. Returns notes."""
    notes = []
    out_name = chain_name if len(chain_name) == 1 else 'A'
    if out_name != chain_name:
        notes.append(f'chain written as {out_name} (PDB format allows one character)')
    st = gemmi.Structure()
    st.name = stem.upper()
    st.info['_entry.id'] = stem.upper()
    st.info['_struct.title'] = f'{title} - chain {chain_name}, renumbered, centre of mass at origin'
    if method:
        st.info['_exptl.method'] = method
    model = gemmi.Model(1)
    ch = gemmi.Chain(out_name)
    n = 0
    for r in polymer:
        n += 1
        r.seqid, r.label_seq, r.subchain, r.entity_id = gemmi.SeqId(n, ' '), n, '', ''
        ch.add_residue(r)
    for r in hets + waters:
        n += 1
        r.seqid, r.label_seq, r.subchain, r.entity_id = gemmi.SeqId(n, ' '), None, '', ''
        ch.add_residue(r)
    model.add_chain(ch)
    st.add_model(model)
    atoms = [a for r in st[0][0] for a in r]
    mass = np.array([a.element.weight for a in atoms])
    xyz = np.array([a.pos.tolist() for a in atoms])
    center = (mass[:, None] * xyz).sum(axis=0) / mass.sum()
    for a, p in zip(atoms, xyz - center):
        a.pos = gemmi.Position(*p)
    st.setup_entities()
    for ent in st.entities:
        if ent.entity_type == gemmi.EntityType.Polymer:
            ent.polymer_type = polymer_type
            ent.full_sequence = [r.name for r in polymer]
    st.make_mmcif_document().write_file(os.path.join(folder, f'{stem}.cif'))
    options = gemmi.PdbWriteOptions()
    options.cryst1_record = False
    try:
        st.write_pdb(os.path.join(folder, f'{stem}.pdb'), options)
    except (RuntimeError, ValueError) as exc:
        notes.append(f'PDB file not written: {exc}')
    return notes


# ---- one entry --------------------------------------------------------------------------

def one_letter(name):
    info = gemmi.find_tabulated_residue(name)
    code = info.one_letter_code.upper() if info is not None else 'X'
    return code if code.strip() else 'X'


def process_entry(path, do_rebuild, retry_note, out_root):
    pdb_id = entry_id(path)
    st, descriptions, scheme = read_entry(path)
    resolve_altlocs(st[0])
    chains = classify_chains(st, descriptions)
    assign_water_roles(st, chains)
    title = st.info['_struct.title'] if '_struct.title' in st.info else ''
    method = st.info['_exptl.method'] if '_exptl.method' in st.info else ''
    base = {'pdb_id': pdb_id, 'source_file': os.path.basename(path)}
    rows = {}
    groups = {}
    for info in chains:
        roles = Counter(r for r in info['roles'] if r)
        hets = Counter(r.name for r in info['hets'])
        row = dict(base, pdbid_chainid=f'{pdb_id}_{info["name"]}', chain=info['name'], chain_type=info['kind'],
                   het_count=sum(hets.values()), hets='; '.join(f'{k}:{v}' for k, v in hets.items()),
                   water_significant=sum(roles.values()), water_total=len(info['waters']),
                   sequence_length=len(info['polymer']),
                   sequence=''.join(one_letter(r.name) for r in info['polymer']), notes=retry_note or '')
        row.update({f'water_{k}': roles[k] for k in WATER_ROLES})
        if info['kind'] == 'none':
            row['status'] = 'skipped: no ATOM records'
        elif info['kind'] != 'protein':
            row['status'] = f'skipped: {info["kind"]}'
        else:
            key = (tuple(r.name for r in info['polymer']), tuple(sorted(hets.items())), tuple(sorted(roles.items())))
            groups.setdefault(key, []).append(info)
        rows[info['name']] = row
    tmp_root = os.path.join(out_root, f'.tmp_{pdb_id}_{os.getpid()}')
    os.makedirs(tmp_root, exist_ok=True)
    written = []
    try:
        for members in groups.values():
            rep = max(members, key=lambda c: (sum(len(r) for r in c['polymer']),
                                              -np.mean([a.b_iso for r in c['polymer'] for a in r] or [0]),
                                              -c['order']))
            stem = f'{pdb_id}_{rep["name"]}'
            ids = [f'{pdb_id}_{c["name"]}' for c in members]
            for c in members:
                rows[c['name']].update(status='unique' if c is rep else 'duplicate', unique_chain=stem,
                                       group_members=','.join(ids), n_members=len(ids))
            notes = []
            polymer = [r.clone() for r in rep['polymer']]
            seqres = [gemmi.Entity.first_mon(m) for m in rep['entity'].full_sequence] if rep['entity'] else []
            built = added = 0
            if do_rebuild:
                polymer, built, added = rebuild(polymer, seqres, notes)
            hets = [r.clone() for r in rep['hets']]
            kept = [(r.clone(), role) for r, role in zip(rep['waters'], rep['roles']) if role]
            folder = os.path.join(tmp_root, stem)
            os.makedirs(folder, exist_ok=True)
            write_residue_map(os.path.join(folder, f'{stem}_residue_map.csv'), scheme.get(rep['name'], []),
                              seqres, {r.label_seq for r in rep['polymer']}, polymer, hets, kept)
            notes += write_chain(folder, stem, rep['name'], title, method, polymer, hets, [w for w, _ in kept],
                                 rep['entity'].polymer_type if rep['entity'] else gemmi.PolymerType.PeptideL)
            row = rows[rep['name']]
            row.update(residues_rebuilt=built, atoms_added=added,
                       notes='; '.join(x for x in [row['notes']] + notes if x))
            written.append(stem)
        for stem in written:
            dest = os.path.join(out_root, stem)
            if os.path.exists(dest):
                shutil.rmtree(dest)
            os.replace(os.path.join(tmp_root, stem), dest)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
    return list(rows.values())


def worker(path, do_rebuild, retry_note, out_root, conn):
    signal.signal(signal.SIGINT, signal.SIG_IGN)   # the main process handles Ctrl-C
    try:
        os.nice(5)                                  # keep the desktop responsive
    except OSError:
        pass
    try:
        conn.send({'ok': True, 'rows': process_entry(path, do_rebuild, retry_note, out_root)})
    except Exception as exc:  # any failure becomes an error row, never a crash of the run
        conn.send({'ok': False, 'error': f'{type(exc).__name__}: {str(exc) or "unreadable or unexpected content"}',
                   'trace': traceback.format_exc(limit=4)})
    finally:
        conn.close()


# ---- process supervision -------------------------------------------------------------------

def mem_available_gb():
    with open('/proc/meminfo') as fh:
        for line in fh:
            if line.startswith('MemAvailable:'):
                return int(line.split()[1]) / 1024 ** 2
    return float('inf')


def process_tree(pid):
    pids, todo = [], [pid]
    while todo:
        p = todo.pop()
        pids.append(p)
        try:
            for task in os.listdir(f'/proc/{p}/task'):
                with open(f'/proc/{p}/task/{task}/children') as fh:
                    todo += [int(c) for c in fh.read().split()]
        except OSError:
            pass
    return pids


def tree_rss_gb(pid):
    total = 0
    for p in process_tree(pid):
        try:
            with open(f'/proc/{p}/status') as fh:
                total += next((int(line.split()[1]) for line in fh if line.startswith('VmRSS:')), 0)
        except OSError:
            pass
    return total / 1024 ** 2


def kill_tree(pid):
    for p in reversed(process_tree(pid)):
        try:
            os.kill(p, signal.SIGKILL)
        except OSError:
            pass


def find_files(folder):
    names = sorted(os.listdir(folder))
    return [os.path.join(folder, n) for n in names
            if n.lower().endswith(EXTENSIONS) and os.path.isfile(os.path.join(folder, n))]


def main():
    ap = argparse.ArgumentParser(description='Unique protein chains (ATOM-record sequence + listed HETATM groups + '
                                             'significant waters), rebuilt, renumbered, centred, with domains.')
    ap.add_argument('mmcif_files_location', nargs='?', default=mmcif_files_location,
                    help=f'folder with .cif.gz/.cif files (default: {mmcif_files_location})')
    ap.add_argument('-o', '--unique-chain-location', default=unique_chain_location,
                    help=f'output folder (default: {unique_chain_location})')
    ap.add_argument('--workers', type=int, default=max(1, min(4, (os.cpu_count() or 2) // 2)),
                    help='files processed at the same time (default: half the CPU threads, at most 4)')
    group = ap.add_mutually_exclusive_group()
    group.add_argument('--resume', action='store_true', help='skip files already in the summary and append')
    group.add_argument('--overwrite', action='store_true', help='start a new summary even if one exists')
    ap.add_argument('--limit', type=int, help='process only the first N files (for testing)')
    ap.add_argument('--no-rebuild', action='store_true', help='do not rebuild gaps or missing atoms')
    ap.add_argument('--domain-data', default=domain_data_location,
                    help=f'folder of the domain database (default: {domain_data_location})')
    ap.add_argument('--prepare-domain-data', action='store_true',
                    help='download the bulk domain/family files (~14 GB, resumable) and build the database, then stop')
    ap.add_argument('--without-interpro', action='store_true',
                    help='with --prepare-domain-data: do not download the 13 GB InterPro file (not used in the mapping)')
    ap.add_argument('--map-domains', action='store_true',
                    help='write the domain files of the chains already in the output folder, then stop')
    args = ap.parse_args()

    if args.prepare_domain_data:
        prepare_domain_data(os.path.abspath(args.domain_data), with_interpro=not args.without_interpro)
        return
    if args.map_domains:
        domain_db = os.path.join(os.path.abspath(args.domain_data), DOMAIN_DB)
        if not os.path.exists(domain_db):
            sys.exit(f'{domain_db} not found: run once with --prepare-domain-data')
        map_domains(os.path.abspath(args.unique_chain_location), domain_db)
        return
    out_root = os.path.abspath(args.unique_chain_location)
    os.makedirs(out_root, exist_ok=True)
    summary_path = os.path.join(out_root, 'unique_chain_summary.csv')
    log_path = os.path.join(out_root, 'py_unique_chain_selector.log')
    for name in os.listdir(out_root):                # leftovers of an interrupted run
        if name.startswith('.tmp_'):
            shutil.rmtree(os.path.join(out_root, name), ignore_errors=True)
    files = find_files(args.mmcif_files_location)
    if args.limit:
        files = files[:args.limit]
    done = set()
    if os.path.exists(summary_path):
        if args.resume:
            with open(summary_path, newline='', encoding='utf-8') as fh:
                done = {r['source_file'] for r in csv.DictReader(fh)}
        elif not args.overwrite:
            sys.exit(f'{summary_path} exists: use --resume to continue that run or --overwrite to start again')
    files = [f for f in files if os.path.basename(f) not in done]
    print(f'{len(files)} files to process ({len(done)} already done), {args.workers} workers, '
          f'output {out_root}', file=sys.stderr)

    ctx = mp.get_context('fork')
    pending = deque((f, not args.no_rebuild, False, '') for f in files)   # path, rebuild, alone, note
    running = {}
    tqdm.monitor_interval = 0                        # no monitor thread: keeps forking safe
    append = args.resume and os.path.exists(summary_path)
    with open(summary_path, 'a' if append else 'w', newline='', encoding='utf-8') as out, \
            open(log_path, 'a', encoding='utf-8') as log:
        writer = csv.DictWriter(out, SUMMARY_COLUMNS, lineterminator='\n', extrasaction='ignore')
        if not append:
            writer.writeheader()
        bar = tqdm(total=len(files), unit='file', dynamic_ncols=True,
                   bar_format='{n_fmt}/{total_fmt} files |{bar}| {elapsed}<{remaining} {desc}')

        def show(last=None):
            if last:
                show.last = last
            names = ', '.join(os.path.basename(j['path']) for j in running.values())
            bar.set_description_str(f'| done: {show.last} | running: {names or "-"}')
        show.last = '-'

        def finish(job, rows, status):
            writer.writerows(rows)
            out.flush()
            log.write(f'{time.strftime("%Y-%m-%d %H:%M:%S")}\t{os.path.basename(job["path"])}\t{status}\t'
                      f'{time.time() - job["start"]:.1f} s\n')
            log.flush()
            show(os.path.basename(job['path']))
            bar.update(1)

        def stop(conn, job):
            kill_tree(job['proc'].pid)
            job['proc'].join(5)
            conn.close()
            running.pop(conn, None)

        def retry_or_fail(job, reason):
            if job['rebuild']:
                pending.append((job['path'], False, True, f'rebuilding skipped after: {reason}'))
                log.write(f'{time.strftime("%Y-%m-%d %H:%M:%S")}\t{os.path.basename(job["path"])}\t'
                          f'retry without rebuilding: {reason}\n')
                log.flush()
            else:
                pdb_id = entry_id(job['path'])
                finish(job, [{'pdbid_chainid': pdb_id, 'pdb_id': pdb_id, 'status': 'error', 'notes': reason,
                              'source_file': os.path.basename(job['path'])}], f'error: {reason}')

        try:
            while pending or running:
                alone = any(j['alone'] for j in running.values())
                while pending and not alone and len(running) < args.workers:
                    path, do_rebuild, run_alone, note = pending[0]
                    if running and (run_alone or mem_available_gb() < MIN_FREE_MEMORY_GB):
                        break
                    pending.popleft()
                    recv, send = ctx.Pipe(duplex=False)
                    proc = ctx.Process(target=worker, args=(path, do_rebuild, note, out_root, send))
                    proc.start()
                    send.close()
                    running[recv] = {'path': path, 'rebuild': do_rebuild, 'alone': run_alone, 'proc': proc,
                                     'start': time.time()}
                    alone = run_alone
                    show()
                for conn in mp.connection.wait(list(running), timeout=1.0):
                    job = running.pop(conn)
                    try:
                        result = conn.recv()
                    except (EOFError, OSError):
                        result = None
                    conn.close()
                    job['proc'].join(30)
                    if job['proc'].is_alive():
                        kill_tree(job['proc'].pid)
                        job['proc'].join(5)
                    if result is None:
                        retry_or_fail(job, f'worker process ended unexpectedly (exit code {job["proc"].exitcode})')
                    elif result['ok']:
                        finish(job, result['rows'], 'ok')
                    else:
                        pdb_id = entry_id(job['path'])
                        finish(job, [{'pdbid_chainid': pdb_id, 'pdb_id': pdb_id, 'status': 'error',
                                      'notes': result['error'], 'source_file': os.path.basename(job['path'])}],
                               f'error: {result["error"]}')
                        log.write(''.join(f'    {line}\n' for line in result.get('trace', '').splitlines()))
                        log.flush()
                now = time.time()
                for conn, job in list(running.items()):
                    if now - job['start'] > FILE_TIMEOUT:
                        stop(conn, job)
                        retry_or_fail(job, f'time limit ({FILE_TIMEOUT} s)')
                    elif tree_rss_gb(job['proc'].pid) > WORKER_MAX_MEMORY_GB:
                        stop(conn, job)
                        retry_or_fail(job, f'memory limit ({WORKER_MAX_MEMORY_GB} GB)')
                free = mem_available_gb()
                if running and (free < LOW_MEMORY_GB and len(running) > 1 or free < LOW_MEMORY_GB / 2):
                    conn, job = max(running.items(), key=lambda kv: tree_rss_gb(kv[1]['proc'].pid))
                    stop(conn, job)
                    if len(running) or not job['alone']:
                        pending.appendleft((job['path'], job['rebuild'], True, ''))   # rerun alone
                        log.write(f'{time.strftime("%Y-%m-%d %H:%M:%S")}\t{os.path.basename(job["path"])}\t'
                                  f'stopped at {free:.1f} GB free memory; rerun alone\n')
                        log.flush()
                    else:
                        retry_or_fail(job, f'not enough free memory ({free:.1f} GB available)')
        except KeyboardInterrupt:
            for conn, job in list(running.items()):
                stop(conn, job)
            bar.close()
            print(f'\nstopped; rerun with --resume to continue ({summary_path})', file=sys.stderr)
            sys.exit(130)
        bar.close()
    print(f'done: summary in {summary_path}', file=sys.stderr)


if __name__ == '__main__':
    main()
