# script/a.sh — HADDOCK3 automation + fpocket

This repository contains a Bash orchestration script (`script/a.sh`) that prepares and runs protein–protein docking workflows using HADDOCK3. It supports both multichain and split-chain inputs, integrates optional fpocket pocket detection, and auto-generates restraints when possible.

Key features:

- Clean & split PDB partners with automated chain renumbering.
- Auto body-restraints for split chains (via `haddock3-restraints restrain_bodies`).
- Optional fpocket integration (`--use-fpocket`) producing pocket-based annotations stored under `computational_data/<label>/`.
- Flexible restraint ingestion: manual tables, fpocket, user-provided computational/experimental metadata.
- Blind docking safeguards (`cmrest` + automatic `ranair` whenever no restraints survive).
- Config writer covering `[topoaa]`, `[rigidbody]`, `[seletop]`, `[flexref]`, `[mdref]`, `[emref]`; analysis executed by HADDOCK3 after completion.
- Multi-pair execution, dry-run support, and summary logging.

Refer to **us manual.md** for detailed usage, CLI options, workflow examples, and fpocket notes.

## Quick start

```bash
conda activate haddock3  # ensure haddock3 & fpocket on PATH
./script/a.sh --auto-partners protein_test
./script/a.sh --use-fpocket --auto-partners protein_test
./script/a.sh --input-mode split --group-bodies P1CLL=ref.pdb --group-bodies P1ZG4=ref2.pdb --use-fpocket
```

Results appear under `result/<project>/PAIR_*`, with HADDOCK outputs in `run_*`, analysis reports, and a `SUMMARY.txt` in the project root.

See `us manual.md` for more information.
