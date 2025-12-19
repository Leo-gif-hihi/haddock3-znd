# Migration Guide: fpocket → ARCTIC3D

This guide helps you migrate existing workflows from fpocket to ARCTIC3D.

## Quick Migration

### Step 1: Update Command

Simply replace `--use-fpocket` with `--use-arctic3d`:

```bash
# Before (fpocket)
./script/a.sh \
  --partner A=proteinA.pdb \
  --partner B=proteinB.pdb \
  --use-fpocket

# After (ARCTIC3D)
./script/a.sh \
  --partner A=proteinA.pdb \
  --partner B=proteinB.pdb \
  --use-arctic3d
```

### Step 2: Install Dependencies

```bash
# Install ARCTIC3D
pip install arctic3d

# Install BLAST (optional but recommended)
sudo apt-get install ncbi-blast+  # Ubuntu/Debian
# or
conda install -c bioconda blast
```

### Step 3: (Optional) Setup Local BLAST

For faster uniprot ID extraction:

```bash
# Download SwissProt database
mkdir -p ~/blast_db
cd ~/blast_db
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
tar -xzf swissprot.tar.gz

# Use in commands
./script/a.sh ... --use-arctic3d --arctic3d-blast-db ~/blast_db/swissprot
```

## Option Mapping

| fpocket Option | ARCTIC3D Equivalent | Notes |
|----------------|---------------------|-------|
| `--use-fpocket` | `--use-arctic3d` | Enable feature |
| N/A (pocket selection was automatic) | `--arctic3d-prob 0.4` | Control restraint quality |
| N/A | `--arctic3d-blast-db <path>` | Speed up uniprot lookup |

**Removed fpocket-specific options:**
- No equivalent for fpocket top-N pockets
- No equivalent for druggability score
- No equivalent for minimum volume

These are replaced by ARCTIC3D's probability threshold which controls interface confidence.

## Behavioral Changes

### 1. Restraint Generation Logic

**fpocket:**
- Identified binding pockets geometrically
- Generated restraints for pocket residues
- Single-structure analysis

**ARCTIC3D:**
- Identifies interfaces from database of known interactions
- Generates restraints for interface residues
- Multi-structure, knowledge-based analysis

### 2. Multi-Chain Handling

**fpocket:**
- Treated all chains as single entity
- Could not distinguish chain-specific pockets

**ARCTIC3D:**
- Analyzes each chain independently
- Can exclude specific chains from analysis
- Generates chain-specific restraints

**Example:** Docking complex A+B to ligand C

```bash
# fpocket (old)
# Would find pockets in A+B combined
# Could include A-B interface pockets (unwanted)

# ARCTIC3D (new)
# Analyzes chain A (excluding B)
# Analyzes chain B (excluding A)
# Generates A-C and B-C restraints (correct!)
```

### 3. PDB Requirements

**fpocket:**
- Only needed 3D coordinates
- No sequence information required

**ARCTIC3D:**
- Needs uniprot ID extraction
- Requires PDB headers OR sequence BLAST
- May need internet access

**Best practice:** Add DBREF records to your PDBs:
```
DBREF  1ABC A    1   100  UNP    P12345   PROT_HUMAN       1    100
```

### 4. Output Structure

**fpocket:**
```
prepared/A/computational_data/A/fpocket/
├── raw/
│   └── pockets/
└── A_fpocket.json
```

**ARCTIC3D:**
```
prepared/A/computational_data/A/arctic3d/
├── chain_manifest.json
└── chain_A/
    ├── arctic3d-P12345/
    │   └── (full ARCTIC3D output)
    └── arctic3d_A.log
```

### 5. Restraint Files

**fpocket:**
- Generated passive restraints based on pocket residues
- Integrated into auto restraints pipeline

**ARCTIC3D:**
- Generates explicit ambiguous/unambiguous restraint files
- Separate from auto restraints
- Named `arctic3d_PAIR_<name>_ambig.tbl`

## Common Migration Scenarios

### Scenario 1: Simple Protein-Protein Docking

**Before:**
```bash
./script/a.sh \
  --partner A=receptor.pdb \
  --partner B=ligand.pdb \
  --use-fpocket \
  --v 3
```

**After:**
```bash
./script/a.sh \
  --partner A=receptor.pdb \
  --partner B=ligand.pdb \
  --use-arctic3d \
  --arctic3d-prob 0.4 \
  --v 3
```

**Changes:**
- Replace `--use-fpocket` with `--use-arctic3d`
- Add `--arctic3d-prob` to control quality
- Everything else stays the same

### Scenario 2: Batch Processing

**Before:**
```bash
for receptor in receptors/*.pdb; do
  for ligand in ligands/*.pdb; do
    ./script/a.sh \
      --partner A=$receptor \
      --partner B=$ligand \
      --use-fpocket \
      --out results/$(basename $receptor .pdb)_vs_$(basename $ligand .pdb)
  done
done
```

**After:**
```bash
# Add one-time BLAST DB setup
BLAST_DB=~/blast_db/swissprot

for receptor in receptors/*.pdb; do
  for ligand in ligands/*.pdb; do
    ./script/a.sh \
      --partner A=$receptor \
      --partner B=$ligand \
      --use-arctic3d \
      --arctic3d-blast-db $BLAST_DB \
      --out results/$(basename $receptor .pdb)_vs_$(basename $ligand .pdb)
  done
done
```

**Changes:**
- Setup local BLAST DB once
- Use `--arctic3d-blast-db` for speed
- Otherwise identical

### Scenario 3: Complex with Manual Restraints

**Before:**
```bash
./script/a.sh \
  --partner A=complex.pdb \
  --partner B=ligand.pdb \
  --use-fpocket \
  --ambig manual_restraints.tbl \
  --v 3
```

**After:**
```bash
./script/a.sh \
  --partner A=complex.pdb \
  --partner B=ligand.pdb \
  --use-arctic3d \
  --ambig manual_restraints.tbl \
  --v 3
```

**Changes:**
- Only the flag changes
- Manual restraints still work
- Will be combined with ARCTIC3D restraints

## Troubleshooting Migration Issues

### Issue 1: "arctic3d command not found"

**Solution:**
```bash
pip install arctic3d
# or
conda install -c conda-forge arctic3d
```

### Issue 2: "Could not determine uniprot ID"

**Possible causes:**
1. PDB lacks DBREF records
2. BLAST database not accessible
3. Network issues (remote BLAST)

**Solutions:**

**Add DBREF to PDB:**
```python
# Python script to add DBREF
import sys
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

pdb_file = sys.argv[1]
# ... extract sequence, run BLAST, add DBREF ...
```

**Use local BLAST:**
```bash
--arctic3d-blast-db ~/blast_db/swissprot
```

**Manual override** (modify script):
```bash
# In extract_uniprot_id function, add hardcoded mapping:
case "$label" in
  A) echo "P12345" ;;
  B) echo "Q98765" ;;
esac
```

### Issue 3: "ARCTIC3D takes too long"

**Cause:** Remote database queries

**Solutions:**
1. Use local BLAST DB (faster uniprot extraction)
2. Pre-compute ARCTIC3D runs separately
3. Run pairs in parallel

### Issue 4: "No restraints generated"

**Possible causes:**
1. Probability threshold too high
2. No interface data in ARCTIC3D database
3. ARCTIC3D analysis failed

**Solutions:**

**Lower threshold:**
```bash
--arctic3d-prob 0.3  # or even 0.2
```

**Check ARCTIC3D output:**
```bash
ls result/PROJECT/prepared/*/computational_data/*/arctic3d/chain_*/arctic3d-*/
```

**Check logs:**
```bash
cat result/PROJECT/prepared/A/computational_data/A/arctic3d/chain_A/arctic3d_A.log
```

**Fallback to manual restraints:**
```bash
--ambig my_manual_restraints.tbl
```

### Issue 5: "Different results than fpocket"

**Expected!** ARCTIC3D and fpocket use different methods:

- **fpocket**: Geometric pocket detection → may find buried pockets
- **ARCTIC3D**: Interface database → only known interaction sites

**Recommendations:**
1. ARCTIC3D generally more specific for protein-protein interfaces
2. For small molecule binding, consider both methods
3. Can combine both approaches:
   ```bash
   # Run both separately
   ./script/a.sh ... --use-arctic3d --out results_arctic3d
   
   # (After installing old version)
   # ./old_script/a.sh ... --use-fpocket --out results_fpocket
   
   # Compare and combine restraints manually
   ```

## Performance Comparison

| Metric | fpocket | ARCTIC3D |
|--------|---------|----------|
| **Speed** | ~5-10 sec/structure | ~30-60 sec/chain |
| **Memory** | Low (~100 MB) | Moderate (~500 MB) |
| **Network** | None | Required (can cache) |
| **Accuracy** | Geometry-based | Knowledge-based |

**Recommendations:**
- **Speed critical?** Setup local BLAST DB
- **High-throughput?** Pre-compute ARCTIC3D runs
- **Offline work?** Not possible with ARCTIC3D (needs DB access)

## Best Practices After Migration

1. **Always use local BLAST:**
   ```bash
   --arctic3d-blast-db ~/blast_db/swissprot
   ```

2. **Start with default threshold (0.4):**
   ```bash
   --arctic3d-prob 0.4
   ```

3. **Check logs on first run:**
   ```bash
   tail -f result/PROJECT/prepared/*/computational_data/*/arctic3d/chain_*/arctic3d_*.log
   ```

4. **Verify restraints generated:**
   ```bash
   ls -lh result/PROJECT/PAIR_*/restraints/arctic3d_*.tbl
   ```

5. **Compare with manual restraints:**
   ```bash
   # Use both if available
   --use-arctic3d --ambig experimental_data.tbl
   ```

## Rollback Plan

If you need to revert to fpocket:

1. **Restore old script:**
   ```bash
   git checkout <old_commit> -- script/a.sh
   ```

2. **Or keep both versions:**
   ```bash
   cp script/a.sh script/a_arctic3d.sh
   git checkout <old_commit> -- script/a.sh
   mv script/a.sh script/a_fpocket.sh
   ```

3. **Use appropriate version:**
   ```bash
   ./script/a_arctic3d.sh ... --use-arctic3d
   ./script/a_fpocket.sh ... --use-fpocket
   ```

## Getting Help

- **Migration issues**: Check logs in `result/PROJECT/prepared/*/computational_data/*/arctic3d/`
- **ARCTIC3D problems**: https://github.com/haddocking/arctic3d/issues
- **Script bugs**: Check `CHANGES_SUMMARY.md` for known issues
- **HADDOCK3 issues**: https://github.com/haddocking/haddock3/issues

## Summary Checklist

- [ ] Install ARCTIC3D (`pip install arctic3d`)
- [ ] Install BLAST (`apt-get install ncbi-blast+`)
- [ ] Download BLAST database (optional but recommended)
- [ ] Replace `--use-fpocket` with `--use-arctic3d`
- [ ] Add `--arctic3d-blast-db` for speed
- [ ] Test on small dataset first
- [ ] Check logs for errors
- [ ] Verify restraints generated
- [ ] Compare results with fpocket (optional)
- [ ] Update scripts/pipelines
- [ ] Document changes in your workflow

**Migration should take: ~15 minutes (+ BLAST DB download time)**
