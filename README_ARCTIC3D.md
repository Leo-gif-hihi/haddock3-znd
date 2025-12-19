# ARCTIC3D Integration - Complete Summary

## What Was Changed

The HADDOCK3 automation script (`script/a.sh`) has been updated to use **ARCTIC3D** instead of **fpocket** for generating protein-protein interface restraints. This addresses the key limitation that `arctic3d-restraints` only accepts two inputs at a time, implementing a robust "divide and conquer" strategy for multi-chain complexes.

## Key Features

✅ **Automatic Uniprot ID extraction** from PDB headers or BLAST  
✅ **Multi-chain aware** - processes each chain independently  
✅ **Exclusion strategy** - prevents unwanted intra-complex interfaces  
✅ **Automatic restraint merging** - combines all chain pairs into master files  
✅ **Fully automated** - no manual intervention required  
✅ **Backward compatible** - existing workflows still work  

## Files Modified

- `script/a.sh` - Main automation script (~500 lines of changes)

## Files Created

1. **ARCTIC3D_QUICKSTART.md** - Quick reference for common use cases
2. **ARCTIC3D_USAGE.md** - Comprehensive user documentation
3. **CHANGES_SUMMARY.md** - Technical implementation details
4. **MIGRATION_GUIDE.md** - Guide for migrating from fpocket
5. **README.md** (this file) - Complete summary

## Quick Start

### Basic Usage

```bash
./script/a.sh \
  --partner A=protein_A.pdb \
  --partner B=protein_B.pdb \
  --use-arctic3d \
  --v 3
```

### Multi-Chain Complex

```bash
./script/a.sh \
  --partner Complex=heterodimer.pdb \
  --partner Ligand=ligand.pdb \
  --use-arctic3d \
  --arctic3d-prob 0.4 \
  --v 3
```

## How It Works

### The "Divide and Conquer" Strategy

**Problem:** ARCTIC3D can only analyze two partners at a time, but protein complexes have multiple chains.

**Solution:** 

1. **Divide**: Process each chain separately
   - Extract uniprot ID for chain A
   - Extract uniprot ID for chain B
   - Run ARCTIC3D on A (excluding B)
   - Run ARCTIC3D on B (excluding A)

2. **Generate**: Create pairwise restraints
   - A vs Ligand → restraints_A_L
   - B vs Ligand → restraints_B_L

3. **Conquer**: Merge all restraints
   - Combine all .tbl files
   - Create master ambig/unambig files
   - Feed to HADDOCK3

### Example Workflow

For complex with chains A+B docking to ligand C:

```
Input: heterodimer.pdb (chains A, B)
       ligand.pdb (chain C)

Step 1: Extract Uniprot IDs
  Chain A → Q8BHN3
  Chain B → O08795
  Chain C → P12345

Step 2: Run ARCTIC3D (with exclusions)
  arctic3d Q8BHN3 --chain A --out_uniprot=O08795
  arctic3d O08795 --chain B --out_uniprot=Q8BHN3
  arctic3d P12345 --chain C

Step 3: Generate Restraints
  arctic3d-restraints Q8BHN3 P12345 → A_C_restraints.tbl
  arctic3d-restraints O08795 P12345 → B_C_restraints.tbl

Step 4: Merge
  all_restraints.tbl = A_C + B_C

Result: Single restraint file for HADDOCK3
```

## Installation

### Prerequisites

```bash
# ARCTIC3D
pip install arctic3d

# BLAST (recommended)
sudo apt-get install ncbi-blast+  # Ubuntu/Debian
# or
conda install -c bioconda blast

# HADDOCK3
pip install haddock3
```

### Optional: Local BLAST Database

```bash
mkdir -p ~/blast_db && cd ~/blast_db
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
tar -xzf swissprot.tar.gz

# Use with:
--arctic3d-blast-db ~/blast_db/swissprot
```

## Documentation Guide

| Document | Use When... |
|----------|-------------|
| **ARCTIC3D_QUICKSTART.md** | You want to start using ARCTIC3D immediately |
| **ARCTIC3D_USAGE.md** | You need detailed documentation and examples |
| **CHANGES_SUMMARY.md** | You want to understand the implementation |
| **MIGRATION_GUIDE.md** | You're migrating from fpocket to ARCTIC3D |
| **README.md** (this file) | You want an overview of everything |

## Common Commands

### Simple Docking

```bash
./script/a.sh \
  --partner A=proteinA.pdb \
  --partner B=proteinB.pdb \
  --use-arctic3d
```

### With All Options

```bash
./script/a.sh \
  --partner Receptor=receptor.pdb \
  --partner Ligand=ligand.pdb \
  --use-arctic3d \
  --arctic3d-prob 0.4 \
  --arctic3d-blast-db ~/blast_db/swissprot \
  --v 3 \
  --ncores 16 \
  --out results \
  --project my_docking
```

### Run HADDOCK3

```bash
# Single pair
cd results/my_docking/PAIR_Receptor_vs_Ligand
haddock3 haddock3.cfg

# All pairs in parallel
find results/my_docking -name 'haddock3.cfg' | \
  parallel -j 5 'cd {//} && haddock3 {/} > run.log 2>&1'
```

## Key Advantages Over fpocket

| Feature | fpocket | ARCTIC3D |
|---------|---------|----------|
| **Basis** | Geometric pockets | Known interfaces |
| **Quality** | Structure-based | Knowledge-based |
| **Multi-chain** | ❌ Single entity | ✅ Chain-aware |
| **Exclusions** | ❌ Not supported | ✅ Via --out_uniprot |
| **Specificity** | Pocket-centric | Interface-centric |
| **For PPI** | ⚠️ May miss interfaces | ✅ Designed for PPI |

## Critical Design Features

### 1. Uniprot ID Extraction

**Three-tier fallback system:**
1. DBREF records in PDB header
2. COMPND fields with uniprot references
3. BLAST search (local or remote)

### 2. Chain-Level Processing

**Each chain processed independently:**
- Extracts its own uniprot ID
- Runs its own ARCTIC3D analysis
- Excludes other chains from the same partner
- Generates separate restraints

### 3. Restraint Mapping

**Chain renaming handled automatically:**
- Original chains → mapped chains stored in JSON
- ARCTIC3D uses mapped chain IDs
- Restraints reference correct chains
- No manual intervention needed

### 4. Intelligent Merging

**Automatic classification and merging:**
- Identifies ambiguous vs unambiguous restraints
- Combines by type
- Preserves source information
- Validates format

## Output Structure

```
result/PROJECT_NAME/
├── prepared/
│   ├── PartnerA/
│   │   ├── PartnerA_combined.pdb
│   │   ├── PartnerA_mapping.json
│   │   └── computational_data/
│   │       └── PartnerA/
│   │           └── arctic3d/
│   │               ├── chain_manifest.json
│   │               └── chain_A/
│   │                   ├── arctic3d-Q8BHN3/
│   │                   │   ├── clustered_interfaces.out
│   │                   │   └── ...
│   │                   └── arctic3d_A.log
│   └── PartnerB/
│       └── ...
└── PAIR_PartnerA_vs_PartnerB/
    ├── input/
    │   ├── PartnerA_0_PartnerA_combined.pdb
    │   └── PartnerB_0_PartnerB_combined.pdb
    ├── restraints/
    │   ├── arctic3d_PAIR_PartnerA_vs_PartnerB_ambig.tbl    ← Master ambiguous
    │   ├── arctic3d_PAIR_PartnerA_vs_PartnerB_unambig.tbl  ← Master unambiguous
    │   └── arctic3d_restraints_work/
    │       ├── PartnerA_A_vs_PartnerB_A/
    │       │   └── *.tbl (individual restraints)
    │       └── PartnerA_A_vs_PartnerB_B/
    │           └── *.tbl (individual restraints)
    ├── meta/
    │   ├── PartnerA_mapping.json
    │   └── PartnerB_mapping.json
    └── haddock3.cfg                                         ← Ready to run!
```

## Troubleshooting

### "arctic3d command not found"
```bash
pip install arctic3d
```

### "Could not determine uniprot ID"
```bash
# Add to PDB or use local BLAST
--arctic3d-blast-db ~/blast_db/swissprot
```

### "No restraints generated"
```bash
# Lower threshold
--arctic3d-prob 0.3
```

### "BLAST is too slow"
```bash
# Use local database
--arctic3d-blast-db ~/blast_db/swissprot
```

## Testing

### Test 1: Simple Two-Partner

```bash
./script/a.sh \
  --partner A=examples/1ppe_E.pdb \
  --partner B=examples/1ppe_I.pdb \
  --use-arctic3d \
  --dry-run
```

Expected: Restraints generated in `result/.../PAIR_A_vs_B/restraints/`

### Test 2: Multi-Chain Complex

```bash
./script/a.sh \
  --partner AB=examples/multi_chain.pdb \
  --partner C=examples/ligand.pdb \
  --use-arctic3d \
  --dry-run
```

Expected: Chain-wise restraints merged properly

### Test 3: With Manual Restraints

```bash
./script/a.sh \
  --partner A=proteinA.pdb \
  --partner B=proteinB.pdb \
  --use-arctic3d \
  --ambig manual.tbl \
  --dry-run
```

Expected: Both ARCTIC3D and manual restraints included

## Performance Notes

- **Uniprot extraction**: 5-30 sec per chain (local BLAST faster)
- **ARCTIC3D run**: 30-120 sec per chain (depends on database)
- **Restraint generation**: 10-30 sec per chain pair
- **Total overhead**: ~2-5 min per docking pair (with local BLAST)

**Optimization tips:**
1. Use local BLAST database
2. Run multiple pairs in parallel
3. Cache ARCTIC3D results if re-running

## Known Limitations

1. **Requires Uniprot ID**: Cannot work with novel/unpublished proteins
2. **Needs Internet**: ARCTIC3D queries online databases
3. **Database Coverage**: Limited to proteins with known interfaces
4. **Slower than fpocket**: More computation required

## Future Enhancements

Potential improvements for future versions:

- [ ] Parallel ARCTIC3D execution
- [ ] Cached uniprot ID database
- [ ] Interface visualization
- [ ] Confidence-based filtering
- [ ] Integration with AlphaFold predictions
- [ ] Support for custom interface databases

## Version History

- **v1.0** (Current): Initial ARCTIC3D integration
  - Replace fpocket with ARCTIC3D
  - Multi-chain support
  - Automatic restraint merging
  - Chain-aware exclusions

## Support & Contributing

### Getting Help

- **Documentation**: See files listed above
- **ARCTIC3D issues**: https://github.com/haddocking/arctic3d/issues
- **HADDOCK3 issues**: https://github.com/haddocking/haddock3/issues
- **Script issues**: Check logs in `result/*/prepared/*/computational_data/*/arctic3d/`

### Reporting Bugs

When reporting issues, please include:
1. Full command used
2. Log files from `result/*/prepared/*/computational_data/*/arctic3d/`
3. PDB files (or sample subset)
4. ARCTIC3D version (`arctic3d --version`)
5. HADDOCK3 version (`haddock3 --version`)

### Contributing

Contributions welcome! Areas for improvement:
- Performance optimization
- Additional restraint types
- Better error handling
- More comprehensive testing

## Acknowledgments

- **ARCTIC3D**: Interface prediction tool by HADDOCK team
- **HADDOCK3**: Integrative modeling platform
- **Original fpocket integration**: Basis for this work

## License

Same as parent HADDOCK3 project

---

## Quick Reference Card

```bash
# Basic command
./script/a.sh --partner A=A.pdb --partner B=B.pdb --use-arctic3d

# With options
./script/a.sh --partner A=A.pdb --partner B=B.pdb \
  --use-arctic3d --arctic3d-prob 0.4 --arctic3d-blast-db ~/blast_db/swissprot

# Multi-chain
./script/a.sh --partner Complex=AB.pdb --partner Ligand=C.pdb --use-arctic3d

# Run HADDOCK
cd result/PROJECT/PAIR_*/
haddock3 haddock3.cfg

# Parallel run
find result -name 'haddock3.cfg' | parallel -j 5 'cd {//} && haddock3 {/}'
```

---

**For detailed information, see the respective documentation files.**
