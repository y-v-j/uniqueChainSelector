# uniqueChainSelector

Pick the **unique protein chains** of every structure in a folder of mmCIF files, rebuild
their missing loops and side chains, renumber them from 1, move them to the origin, and
write each one as its own `.cif` / `.pdb` file with **Pfam, CATH, SCOP2 and UniProt**
annotations in the new numbering.

Two chains of the same entry count as the same chain only when all of these match:

- the **ATOM-record sequence**: residues that have coordinates, in file order, with
  residue names as written (modified residues such as `SEP` or `MSE` keep their names);
  unmodelled SEQRES residues do not count
- the **HETATM groups** the chain carries, from a fixed list of 52 biologically relevant
  codes (metal ions, heme, FAD/NAD, nucleotides, glycans, …), counted by code
- the **waters with a structural or functional role**, counted by role (see below)

The whole thing is one Python script, `py_unique_chain_selector.py`, written to run
unattended over the full PDB (~200,000 entries).

## What it does

For every `*.cif.gz` / `*.cif` file:

1. Reads the structure with [gemmi](https://gemmi.readthedocs.io), keeps the **first model**
   and drops hydrogens.
2. Resolves **alternate locations** by occupancy: at each position it keeps the residue type
   and the conformer with the highest occupancy; of two overlapping whole-residue
   alternates (for example `AGS` on altloc A and `ADP` on altloc B) it keeps the better
   occupied one.
3. Splits the entry into author chains (`auth_asym_id`). Protein chains are processed;
   chains with no ATOM records (only HETATM entries) and nucleic-acid chains are listed in
   the summary as skipped.
4. Gives each water a **role** (each water counts once, under the first role it meets):

   | Role | Rule |
   |---|---|
   | `metal` | first coordination shell of a metal of a listed group (free ion or cofactor metal such as Fe of HEM): deposited `metalc` link, or distance ≤ 2.5–3.3 Å depending on the element |
   | `ligand` | H-bond (≤ 3.5 Å) to a listed group's N/O/S **and** to the protein or a second listed group |
   | `buried` | at most 5 % of its solvent-accessible surface open (1.4 Å probe, other waters removed) and at least one H-bond to the protein or a listed group |
   | `interface` | H-bonds to the polymers of two or more chains |
   | `structural` | ≥ 3 H-bonds to protein N/O atoms of at least two residues |

   Only the asymmetric unit is used; contacts with symmetry mates are not seen.
5. Groups the protein chains by their signature (sequence + HETATM groups + water roles).
   Each group is written once; its **representative** is the member with the most observed
   heavy atoms, then the lowest mean B factor, then the first in the file.
6. For each representative:
   - **rebuilds** internal gaps of up to 30 residues (from SEQRES) and missing heavy atoms of
     standard residues with [PDBFixer](https://github.com/openmm/pdbfixer) /
     [OpenMM](https://openmm.org), one small local fragment at a time. Every rebuilt loop is
     checked (peptide C–N ≤ 1.8 Å, three tries with different seeds); a loop that does not
     close is left open and noted. Unmodelled terminal tails are not built. Added atoms get
     **occupancy 0.00** so they stay recognisable.
   - **renumbers**: protein residues 1…N, then the HETATM groups, then the waters
   - moves the **centre of mass** (all written atoms, atomic masses) to the origin
   - writes the chain with its listed HETATM groups and significant waters, plus a
     **residue map** that links old and new numbering

Domain mapping is a separate, quick step that reads the residue maps, so annotations can be
added or refreshed without processing the structures again.

## Installation

Linux, Python ≥ 3.12. The script runs on the CPU; no GPU is needed.

```bash
git clone https://github.com/y-v-j/uniqueChainSelector.git
cd uniqueChainSelector
conda env create -f environment.yml
conda activate unique-chain-selector
```

All packages come from conda-forge: gemmi, numpy, scipy, tqdm, openmm and pdbfixer.

## Usage

There are three steps. Steps 1 and 2 are independent and can run in either order.

**1. Unique chains** (the long step)

```bash
python py_unique_chain_selector.py path/to/mmcif_files -o path/to/unique_chains --workers 4
```

**2. Domain and family data** (once; downloads the bulk files and builds a local database)

```bash
python py_unique_chain_selector.py --prepare-domain-data --domain-data path/to/unique_chain_references
```

This downloads about 40 MB of mapping tables plus InterPro's `protein2ipr.dat.gz` (13 GB),
which is kept for reference but not used in the mapping. Add `--without-interpro` to skip
it. Downloads resume if interrupted.

**3. Domain mapping** (seconds to minutes)

```bash
python py_unique_chain_selector.py --map-domains -o path/to/unique_chains --domain-data path/to/unique_chain_references
```

Run step 3 again whenever the domain database has been rebuilt.

The default folders (`mmcif_files/`, `unique_chains/`, `unique_chain_references/`, relative
to the current directory) and the HETATM list are set at the top of the script.

### Options

| Option | Meaning |
|---|---|
| `mmcif_files_location` | input folder with `.cif.gz` / `.cif` files |
| `-o`, `--unique-chain-location` | output folder |
| `--workers N` | files processed at the same time (default: half the CPU threads, at most 4) |
| `--resume` | continue a stopped run: skip the files already in the summary |
| `--overwrite` | start a new summary even if one exists |
| `--limit N` | process only the first N files (for testing) |
| `--no-rebuild` | do not rebuild gaps or missing atoms (much faster) |
| `--domain-data DIR` | folder of the domain database |
| `--prepare-domain-data` | download the bulk files and build the domain database, then stop |
| `--without-interpro` | with `--prepare-domain-data`: do not download the 13 GB InterPro file |
| `--map-domains` | write the domain files for the chains in the output folder, then stop |

### Large runs

Each file runs in its own process, and every PDBFixer fragment in a short-lived child
process with a time limit, so one bad or huge file cannot stop the run:

- a file over the time limit (1 h) or the per-worker memory limit (6 GB), or whose process
  dies, is retried once without rebuilding, then recorded as an error
- new files start only while at least 2.5 GB of memory is free; below 1 GB the largest worker
  is stopped and its file rerun alone
- results are written after every file; Ctrl-C stops cleanly and `--resume` carries on

For the whole PDB, run it in the background and, if you like, under a memory cap:

```bash
nohup systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=0 \
  python py_unique_chain_selector.py path/to/mmcif_files -o path/to/unique_chains --workers 4 \
  > run.log 2>&1 &
```

As a guide, a 100-entry sample took about 6 minutes with 4 workers (median under 1 s per
file; entries with many gaps up to 3–4 minutes). Rebuilding takes most of the time;
`--no-rebuild` is much faster.

## Output

```
unique_chains/
├── unique_chain_summary.csv          one row per chain of every entry
├── py_unique_chain_selector.log      one line per file (status, time)
└── 1abc_A/
    ├── 1abc_A.cif                    the chain, renumbered and centred
    ├── 1abc_A.pdb                    same in PDB format
    ├── 1abc_A_residue_map.csv        old and new numbering of every residue
    └── 1abc_A_domains.csv            Pfam / CATH / SCOP2 / UniProt in the new numbering
```

The structure files open in PyMOL, ChimeraX, RasMol and other viewers. A chain whose ID has
more than one character is written as chain `A` (PDB format allows one character); the
summary notes it.

**`unique_chain_summary.csv`**: `pdbid_chainid`, `pdb_id`, `chain`, `chain_type`, `status`
(`unique`, `duplicate`, `skipped: …`, `error`), `unique_chain` (the representative written
for this chain), `group_members`, `n_members`, `sequence_length`, `het_count`, `hets`
(e.g. `HEM:1; PO4:1`), `water_metal`, `water_ligand`, `water_buried`, `water_interface`,
`water_structural`, `water_significant`, `water_total`, `residues_rebuilt`, `atoms_added`,
`domains` (filled by `--map-domains`), `notes` (gaps left open and why, retries),
`source_file`, `sequence` (one-letter, ATOM records).

**`<chain>_residue_map.csv`**: `new_number`, `residue`, `seqres_number`, `author_number`,
`status`. Every SEQRES position is listed as `observed`, `rebuilt` or `missing`, followed by
the HETATM groups (`hetatm`) and waters (`water (role)`) with their old and new numbers.

**`<chain>_domains.csv`**: `database`, `accession`, `name`, `description`,
`residues_new_numbering` (e.g. `1-94,120-286`), `n_residues`, `seqres_range`,
`author_range`, `details`.

## Domain and family data

`--prepare-domain-data` downloads these files once and builds one SQLite database
(`domain_annotations.sqlite`, about 650 MB); `--map-domains` only reads that database.

| Source | Files | Gives |
|---|---|---|
| UniProt | SIFTS `pdb_chain_uniprot` + one UniProt query of all entries with a 3D structure | accession, entry and protein name, ranges |
| Pfam | SIFTS `pdb_pfam_mapping` + `Pfam-A.clans` | family, clan, ranges |
| CATH | CATH-B daily boundaries + `cath-names` | superfamily code and names, domain boundaries |
| SCOP2 | SIFTS `pdb_chain_scop2b_sf_uniprot`, `pdb_chain_scop2_uniprot` + SCOP2 release files | superfamily domains with ranges, family per chain |
| InterPro | `protein2ipr.dat.gz`, `entry.list` | downloaded and kept, not used |

CATH boundaries use author numbering and are converted with each entry's own numbering
table. A few entries have no mappings in SIFTS; their domain files say so.

## Settings

Constants at the top of the script: `hetatm_list`, the water cutoffs (`HBOND_CUTOFF`,
`MIN_HBONDS`, `MAX_EXPOSED`, `METAL_CUTOFF`), rebuilding (`MAX_GAP_LENGTH`,
`MAX_PEPTIDE_BOND`, `GAP_TIMEOUT`) and the safety limits (`FILE_TIMEOUT`,
`WORKER_MAX_MEMORY_GB`, `MIN_FREE_MEMORY_GB`, `LOW_MEMORY_GB`).

## Notes and limits

- Rebuilt loops and side chains are modelled, not observed: their geometry is approximate,
  and they are marked with occupancy 0.00 and listed as `rebuilt` in the residue map.
- Chains are compared exactly. Copies of one protein that differ in their HETATM groups or
  in their significant waters (common in crystal structures) count as different chains.
- Water roles use the asymmetric unit only, so contacts made through crystal symmetry are
  not counted.

## Credits

Data from [SIFTS / PDBe](https://www.ebi.ac.uk/pdbe/docs/sifts/),
[CATH](https://www.cathdb.info), [Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/),
[SCOP2](https://www.ebi.ac.uk/pdbe/scop/), [UniProt](https://www.uniprot.org) and
[InterPro](https://www.ebi.ac.uk/interpro/); please cite them and follow their terms of
use. Structure handling by [gemmi](https://gemmi.readthedocs.io); rebuilding by
[PDBFixer](https://github.com/openmm/pdbfixer) and [OpenMM](https://openmm.org).

## License

MIT; see [LICENSE](LICENSE).
