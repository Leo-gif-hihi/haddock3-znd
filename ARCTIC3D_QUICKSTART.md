# ARCTIC3D Quick Start Guide

## TL;DR

Replace fpocket with ARCTIC3D for better restraint generation:

```bash
# Old way (fpocket)
./script/a.sh --partner A=proteinA.pdb --partner B=proteinB.pdb --use-fpocket

# New way (ARCTIC3D)
./script/a.sh --partner A=proteinA.pdb --partner B=proteinB.pdb --use-arctic3d
```

## Common Use Cases

### 1. Simple Protein-Protein Docking

```bash
./script/a.sh \
  --partner A=receptor.pdb \
  --partner B=ligand.pdb \
  --use-arctic3d \
  --v 3
```

### 2. Multi-Chain Complex + Ligand

Docking a heterodimer (chains A+B) to a ligand (C):

```bash
./script/a.sh \
  --partner AB=complex.pdb \
  --partner C=ligand.pdb \
  --use-arctic3d \
  --v 3
```

**What happens:**
- ARCTIC3D runs on chain A (excluding B)
- ARCTIC3D runs on chain B (excluding A)
- ARCTIC3D runs on chain C
- Generates restraints for A-C and B-C
- Merges all restraints automatically

### 3. With Local BLAST (Faster)

```bash
# First, download BLAST database
mkdir -p /data/blast && cd /data/blast
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
tar -xzf swissprot.tar.gz

# Then use it
./script/a.sh \
  --partner A=proteinA.pdb \
  --partner B=proteinB.pdb \
  --use-arctic3d \
  --arctic3d-blast-db /data/blast/swissprot
```

### 4. Adjust Restraint Stringency

```bash
# More restraints, less confident (default: 0.4)
./script/a.sh --partner A=A.pdb --partner B=B.pdb --use-arctic3d --arctic3d-prob 0.3

# Fewer restraints, more confident
./script/a.sh --partner A=A.pdb --partner B=B.pdb --use-arctic3d --arctic3d-prob 0.5
```

### 5. Complete Example with All Options

```bash
./script/a.sh \
  --partner Receptor=7KBJ.pdb \
  --partner Ligand=compound.pdb \
  --use-arctic3d \
  --arctic3d-prob 0.4 \
  --arctic3d-blast-db /data/blast/swissprot \
  --v 3 \
  --ncores 16 \
  --out results/my_docking \
  --project 7KBJ_docking \
  --reference native.pdb
```

## Key Options

| Option | Description | Default |
|--------|-------------|---------|
| `--use-arctic3d` | Enable ARCTIC3D restraints | disabled |
| `--arctic3d-prob <float>` | Confidence threshold | 0.4 |
| `--arctic3d-blast-db <path>` | Local BLAST database | remote BLAST |
| `--v {1\|2\|3}` | Workflow version | 3 |
| `--ncores <int>` | CPU cores for HADDOCK | 10 |

## Output Files

```
result/PROJECT_NAME/
└── PAIR_A_vs_B/
    ├── restraints/
    │   ├── arctic3d_PAIR_A_vs_B_ambig.tbl     ← Use this in HADDOCK
    │   └── arctic3d_PAIR_A_vs_B_unambig.tbl   ← Use this in HADDOCK
    └── haddock3.cfg                             ← Already configured!
```

## Running HADDOCK3

The script auto-generates `haddock3.cfg` with restraints included.

### Run Single Pair

```bash
cd result/PROJECT_NAME/PAIR_A_vs_B
haddock3 haddock3.cfg
```

### Run All Pairs Sequentially

```bash
for dir in result/PROJECT_NAME/PAIR_*; do
    (cd $dir && haddock3 haddock3.cfg)
done
```

### Run All Pairs in Parallel (Recommended)

```bash
# 5 jobs at a time
find result/PROJECT_NAME -name 'haddock3.cfg' | \
  parallel -j 5 'cd {//} && haddock3 {/} > run.log 2>&1'
```

## Troubleshooting

### "arctic3d command not found"

```bash
pip install arctic3d
```

### "Could not determine uniprot ID"

**Option 1:** Add DBREF records to your PDB
```
DBREF  1ABC A    1   100  UNP    P12345   PROT_HUMAN       1    100
```

**Option 2:** Use local BLAST
```bash
--arctic3d-blast-db /path/to/blast/db
```

### "No restraints generated"

**Try lower threshold:**
```bash
--arctic3d-prob 0.3  # instead of default 0.4
```

**Check logs:**
```bash
# Look for ARCTIC3D errors
grep -r "arctic3d" result/PROJECT_NAME/PAIR_*/restraints/arctic3d_restraints_work/*/arctic3d_*.log
```

### "BLAST takes too long"

Use local BLAST database (see "With Local BLAST" example above)

## Comparison: fpocket vs ARCTIC3D

| Feature | fpocket | ARCTIC3D |
|---------|---------|----------|
| **Basis** | Geometric pockets | Interface databases |
| **Output** | Binding pockets | Protein-protein interfaces |
| **Multi-chain** | ❌ Treats as single unit | ✅ Chain-aware |
| **Exclusions** | ❌ Not supported | ✅ `--out_uniprot` |
| **Quality** | Pocket-based | Interface-based |
| **Speed** | Fast | Moderate (requires DB queries) |

## Need More Help?

- **Full Documentation**: See `ARCTIC3D_USAGE.md`
- **Technical Details**: See `CHANGES_SUMMARY.md`
- **ARCTIC3D Issues**: https://github.com/haddocking/arctic3d/issues
- **HADDOCK3 Issues**: https://github.com/haddocking/haddock3/issues
