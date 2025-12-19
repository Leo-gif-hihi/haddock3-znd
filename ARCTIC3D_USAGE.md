# Using ARCTIC3D for Restraint Generation in HADDOCK3

This document explains how to use ARCTIC3D instead of fpocket to generate high-quality interface-based restraints for protein-protein and protein-ligand docking with HADDOCK3.

## Overview

ARCTIC3D is a tool that identifies protein-protein interaction interfaces by analyzing structural databases and clustering residues that frequently appear at binding interfaces. The updated `a.sh` script now supports ARCTIC3D for generating restraints automatically.

## Key Features

1. **Automatic Uniprot ID Extraction**: Extracts uniprot IDs from PDB headers (DBREF, COMPND) or via BLAST
2. **Multi-Chain Support**: Handles complexes with multiple chains by generating restraints for each chain pair
3. **Exclusion Strategy**: Prevents finding interfaces between chains of the same partner (e.g., A+B complex docking to ligand C)
4. **Automatic Restraint Merging**: Combines restraints from multiple chain pairs into a single file for HADDOCK3

## Installation Requirements

Ensure you have the following tools installed:

```bash
# ARCTIC3D and its dependencies
pip install arctic3d

# BLAST (for sequence-based uniprot ID lookup)
sudo apt-get install ncbi-blast+  # On Ubuntu/Debian
# or
conda install -c bioconda blast

# HADDOCK3
pip install haddock3
```

## Basic Usage

### Simple Two-Partner Docking with ARCTIC3D

```bash
./script/a.sh \
  --partner A=protein_A.pdb \
  --partner B=protein_B.pdb \
  --use-arctic3d \
  --arctic3d-prob 0.4 \
  --v 3
```

### Multi-Chain Complex Docking

When docking a multi-chain complex (e.g., A+B) to a ligand (C):

```bash
./script/a.sh \
  --partner AB=complex_AB.pdb \
  --partner C=ligand.pdb \
  --use-arctic3d \
  --arctic3d-prob 0.4 \
  --v 3
```

The script will:
1. Extract uniprot IDs for chains A and B from complex_AB.pdb
2. Run ARCTIC3D for chain A (excluding interactions with chain B)
3. Run ARCTIC3D for chain B (excluding interactions with chain A)
4. Extract uniprot ID for ligand C
5. Run ARCTIC3D for ligand C
6. Generate restraints for:
   - A vs C
   - B vs C
7. Merge all restraints into single ambiguous/unambiguous files

## Command-Line Options

### ARCTIC3D-Specific Options

- `--use-arctic3d`: Enable ARCTIC3D restraint generation
- `--arctic3d-prob <float>`: Probability threshold for restraints (default: 0.4)
  - Higher values = more confident but fewer restraints
  - Lower values = more restraints but less confident
  - Recommended range: 0.3 - 0.5
- `--arctic3d-blast-db <path>`: Path to local BLAST database for faster uniprot ID lookup
  - If not provided, remote BLAST will be used (slower)
  - Example: `--arctic3d-blast-db /data/blast/swissprot`

### Other Relevant Options

- `--keep-chains-separate`: Don't force chains from same PDB to move together
- `--v {1|2|3}`: Workflow version (use 3 for best results with ARCTIC3D)
- `--pair A,B`: Restrict docking to specific partner pairs

## Workflow Details

### Step 1: Uniprot ID Extraction

The script attempts to extract uniprot IDs in this order:

1. **DBREF records** in PDB file
   ```
   DBREF  7KBJ A    1   423  UNP    Q8BHN3   Q8BHN3_MOUSE     1    423
   ```

2. **COMPND records** containing uniprot references

3. **BLAST search** (local or remote)
   - Extracts sequence from PDB
   - Searches against UniProt/SwissProt database
   - Uses top hit

### Step 2: ARCTIC3D Execution

For each chain in each partner:

```bash
arctic3d <uniprot_id> \
  --pdb_to_use=<prepared_pdb> \
  --chain=<chain_id> \
  --out_uniprot=<other_chains_to_exclude>
```

Example for chain A of a complex with chains A and B:
```bash
arctic3d Q8BHN3 \
  --pdb_to_use=partner_A_combined.pdb \
  --chain=A \
  --out_uniprot=O08795  # Exclude chain B's uniprot
```

### Step 3: Restraint Generation

For each pair of chains between two partners:

```bash
arctic3d-restraints \
  --r1 arctic3d-Q8BHN3 \
  --r2 arctic3d-P12345 \
  --prob_threshold 0.4 \
  --ch1 A \
  --ch2 C \
  --run_dir restraints_A_C
```

### Step 4: Restraint Merging

All generated restraint files are merged into two master files:

- `arctic3d_PAIR_<name>_ambig.tbl`: All ambiguous restraints
- `arctic3d_PAIR_<name>_unambig.tbl`: All unambiguous restraints

The merged files include comments indicating the source of each restraint:

```
! ARCTIC3D-generated ambiguous restraints
! Merged from 2 chain pair(s)
!
! From: AB_A_vs_C_C/restraint_1.tbl
assign (segid A and resid 25) (segid C and resid 10) 2.0 2.0 0.0
! From: AB_B_vs_C_C/restraint_1.tbl
assign (segid B and resid 42) (segid C and resid 15) 2.0 2.0 0.0
```

## Output Structure

```
result/
└── PROJECT_NAME/
    ├── prepared/
    │   ├── A/
    │   │   ├── A_combined.pdb
    │   │   ├── A_mapping.json
    │   │   └── computational_data/
    │   │       └── A/
    │   │           └── arctic3d/
    │   │               ├── chain_manifest.json
    │   │               └── chain_A/
    │   │                   └── arctic3d-Q8BHN3/
    │   └── B/
    │       └── ...
    └── PAIR_A_vs_B/
        ├── input/
        │   ├── A_0_A_combined.pdb
        │   └── B_0_B_combined.pdb
        ├── restraints/
        │   ├── arctic3d_PAIR_A_vs_B_ambig.tbl
        │   ├── arctic3d_PAIR_A_vs_B_unambig.tbl
        │   └── arctic3d_restraints_work/
        │       ├── A_A_vs_B_A/
        │       └── A_A_vs_B_B/
        └── haddock3.cfg
```

## Important Notes

### Chain Renaming and Restraints

The script handles chain renaming automatically:

1. Original PDB chains are mapped to new chain IDs
2. The mapping is stored in `*_mapping.json` files
3. Restraints use the **new chain IDs** from the mapping
4. This ensures restraints remain valid after chain renumbering

### Multi-Chain Complexes

When docking a multi-chain complex to a ligand:

- Each chain of the complex is analyzed separately by ARCTIC3D
- The `--out_uniprot` flag excludes other chains of the same complex
- This prevents finding the intra-complex interface (A-B) when you want the complex-ligand interface (AB-C)
- Restraints are generated for all chain combinations: A-C, B-C
- All restraints are merged into a single file

### Performance Considerations

1. **Local BLAST Database**: Using `--arctic3d-blast-db` significantly speeds up uniprot ID extraction
   
   ```bash
   # Download SwissProt database
   mkdir -p /data/blast
   cd /data/blast
   wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
   tar -xzf swissprot.tar.gz
   
   # Use in script
   --arctic3d-blast-db /data/blast/swissprot
   ```

2. **Parallel Execution**: The script runs ARCTIC3D sequentially for each chain, but you can run multiple pairs in parallel:

   ```bash
   find result/PROJECT_NAME -name 'haddock3.cfg' | \
     parallel -j 5 'cd {//} && haddock3 {/} > run.log 2>&1'
   ```

## Troubleshooting

### Uniprot ID Not Found

If uniprot ID extraction fails:

1. **Check PDB Headers**: Ensure DBREF or COMPND records contain uniprot information
2. **Use Remote BLAST**: Remove `--arctic3d-blast-db` to try remote BLAST
3. **Manual Specification**: Modify the script to hardcode uniprot IDs if needed

### ARCTIC3D Fails

Common issues:

1. **Uniprot ID Not in Database**: The protein may not have interaction data
2. **Network Issues**: ARCTIC3D requires internet access to query databases
3. **PDB Format Issues**: Ensure PDB is properly formatted

### No Restraints Generated

If no restraints are generated:

1. **Lower Probability Threshold**: Try `--arctic3d-prob 0.3` instead of 0.4
2. **Check ARCTIC3D Output**: Look in `arctic3d_restraints_work/*/arctic3d_*.log`
3. **Verify Interface Data**: ARCTIC3D may not have interface data for the proteins

## Advanced Example

Complete example for docking a heterodimer (A+B) to a ligand (C) with all options:

```bash
./script/a.sh \
  --partner AB=7KBJ.pdb \
  --partner C=ligand.pdb \
  --use-arctic3d \
  --arctic3d-prob 0.4 \
  --arctic3d-blast-db /data/blast/swissprot \
  --v 3 \
  --ncores 16 \
  --out result/7KBJ_docking \
  --project 7KBJ_vs_ligand \
  --reference native_complex.pdb
```

This will:
1. Prepare 7KBJ.pdb (chains A and B)
2. Prepare ligand.pdb
3. Extract uniprot IDs using local BLAST
4. Run ARCTIC3D for each chain separately
5. Generate restraints for A-C and B-C
6. Merge all restraints
7. Create HADDOCK3 config with 16 cores
8. Evaluate results against native_complex.pdb

## References

- ARCTIC3D: https://github.com/haddocking/arctic3d
- HADDOCK3: https://github.com/haddocking/haddock3
- Original fpocket-based approach: See `arctic_examples.md`

## Support

For issues related to:
- **Script functionality**: Check script comments and error logs
- **ARCTIC3D**: https://github.com/haddocking/arctic3d/issues
- **HADDOCK3**: https://github.com/haddocking/haddock3/issues
