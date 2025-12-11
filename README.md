# HADDOCK3 automation

This repository contains a Bash orchestration script (`script/a.sh`) that prepares and runs protein–protein docking workflows using HADDOCK3. It supports both multichain and split-chain inputs, integrates optional fpocket pocket detection, and auto-generates restraints when possible.

Key features:

- Clean PDB partners: removes alternative conformations (B, C, etc.), keeps only first conformation.
- **Optional HETATM removal**: Use `--remove-hetatm` to exclude all HETATM records (default: keep all valid molecules).
- **Multi-chain file handling**: PDB files are kept with multiple chains in a single file, with proper TER statements between chains and END statement at the end.
- Auto body-restraints by default (via `haddock3-restraints restrain_bodies`) to keep chains from the same PDB together; disable with `--keep-chains-separate`.
- Optional fpocket integration (`--use-fpocket`) producing pocket-based annotations stored under `computational_data/<label>/`.
- Flexible restraint ingestion: manual tables, fpocket, user-provided computational/experimental metadata.
- Blind docking safeguards (`cmrest` + automatic `ranair` whenever no restraints survive).
- Config writer covering `[topoaa]`, `[rigidbody]`, `[seletop]`, `[flexref]`, `[mdref]`, `[emref]`; analysis executed by HADDOCK3 after completion.
- Multi-pair execution preparation (generates configs for parallel execution).

Refer to **us manual.md** for detailed usage, CLI options, workflow examples, and fpocket notes.

## Quick start

1. **Generate configurations:**
   ```bash
   conda activate haddock3  # ensure haddock3 & fpocket on PATH
   ./script/a.sh --auto-partners protein_test
   ```

2. **Run jobs in parallel (e.g., 5 at a time):**
   ```bash
   ls result/<project_name>/PAIR_*/haddock3.cfg | parallel -j 5 "cd \$(dirname {}) && haddock3 \$(basename {})"
   ```

Results appear under `result/<project>/PAIR_*`, with HADDOCK outputs in `run_*`, analysis reports, and a `SUMMARY.txt` in the project root.

See `us manual.md` for more information.
