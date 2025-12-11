# HADDOCK3 Automation Script — User Manual

This manual describes how to drive `script/a.sh` to prepare and execute HADDOCK3 runs, including the optional fpocket integration. It covers CLI options, pipeline stages, and common workflows.

---

## 1. Prerequisites
- **HADDOCK3** installed and available in the current shell (`haddock3` & `haddock3-restraints`).
- **Python 3**, bash, and common UNIX utilities.
- **fpocket** (optional, only if you use `--use-fpocket`). Verify via `fpocket -v`.
- Ensure input PDB files are accessible. Split and multi-chain modes expect different formats.

---

## 2. Command line interface (summary)

```
./script/a.sh [options]
```

Key options:

| Category | Option | Description |
|----------|--------|-------------|
| **Partner selection** | `--partner A=file.pdb` | Define a partner label and its PDB/manifest. Repeat per partner. |
| | `--auto-partners dir` | Auto-load every `*.pdb` in `dir` (label derived from filename). |
| | `--input-mode {split,multichain}` | Interpret partner specs as per-chain files or multi-chain PDB. |
| | `--group-bodies label=path` | Provide reference PDB/manifest for body restraints (split mode). |
| **Workflow tiers** | `--v {1|2|3}` | Version: V1 blind, V2 computational, V3 experimental priority (default 3). |
| | `--abinitio` | Shortcut: set version 1, high sampling, `ranair=true`. |
| | `--ranair` | Force random AIR generation during rigidbody stage. |
| **Restraints** | `--ambig files` | Manual ambiguous tables (comma-separated). |
| | `--unambig files` | Manual unambiguous tables. |
| | `--computational-dir dir` | Additional annotation directory tree (still supported). |
| | `--experimental-dir dir` | Experimental annotation directory tree (V3). |
| | `--use-fpocket` | Run fpocket on prepared structures to seed computational annotations. |
| **Scheduling** | `--pair L1,L2` | Restrict docking to specific pairs (repeatable). |
| | `--pairs-file file` | Text/CSV file listing pairs. |
| **Execution** | `--out dir`, `--project name` | Control result location. |
| | `--ncores`, `--sampling`, `--reference pdb` | Standard HADDOCK settings. |
| | `--dry-run` | Prepare inputs/config only (no `haddock3`). |
| | `--no-run` | Alias for `--dry-run`. |

Readable help: `./script/a.sh --help`.

---

## 3. Pipeline at a glance

1. **Argument parsing** → collects partners, flags, directories.
2. **Preparation phase** → per partner:
   - Clean PDB (`rename` mode even for multichain to avoid chain collisions).
   - Split by chain (if requested) and generate `*_mapping.json`.
   - Run `haddock3-restraints restrain_bodies` for split partners where applicable.
   - With `--group-bodies`, copy reference PDBs and create extra body locks.
   - If `--use-fpocket`, run fpocket on the cleaned combined PDB, parse top N pockets and export `computational_data/<label>/fpocket/<label>_fpocket.json`.

3. **Pair scheduling** → either all-vs-all or filtered by `--pair`/`--pairs-file`.
4. **Per pair**:
   - Copy sanitized molecules into `PAIR_*/input/` (with unique names).
   - Copy mapping metadata and aggregated restraints.
   - Collect manual, fpocket, and user-supplied computational data → generate ambiguous/unambiguous restraints automatically (versions ≥ 2).
   - **Blind docking safety**: if no restraints remain, `ranair = true` ensures random AIRs (applies for any version once the run is restraint-free).
   - Compose `haddock3.cfg` with `[topoaa]`, `[rigidbody]` (always `cmrest = true`, `ranair` determined at runtime), `[seletop]`, `[flexref]`, `[mdref]`, `[emref]`. CAPRI analysis runs post-hoc.
   - **Note**: The script generates the configuration and directory structure but **does not execute** HADDOCK3 immediately. This allows for efficient parallel execution of multiple pairs.

5. **Outputs**: per pair, you’ll find configs and input files under `PAIR_*`.
6. Final summary printed to stdout and saved as `result/<project>/SUMMARY.txt`, including the command to run the jobs.

---

## 4. Running the jobs (Parallel Execution)

Since the script only prepares the run directories, you need to execute the docking jobs. The script outputs a command at the end, typically using GNU Parallel:

```bash
ls result/<project_name>/PAIR_*/haddock3.cfg | parallel -j 5 "cd \$(dirname {}) && haddock3 \$(basename {})"
```

- Adjust `-j 5` to the number of concurrent jobs you want to run.
- This command finds all generated config files, changes into their directory, and runs `haddock3`.

---

## 5. fpocket integration details

- Enabled with `--use-fpocket`.
- Script checks for `fpocket` in `PATH`. If missing, prints a warning and continues without pockets.
- Runs on each partner’s cleaned combined PDB (`prepared/<label>/<label>_combined.pdb`).
- Keeps top 3 pockets with druggability ≥ 0.5 and volume ≥ 50 Å³ (adjusted in code). The residues inside those pockets are written to JSON:

```json
{
  "label": "P1CLL",
  "pockets": [...],
  "residue_annotations": [
    {"chain": "A", "residue": "123", "confidence": 0.72, "source": "fpocket:pocket1", "score": 0.82, "volume": 76.5},
    ...
  ]
}
```

- Auto-restraint generator consumes these annotations: residues from pocket files behave like computational “active” hints (confidence is derived from fpocket `druggability` or general score).
- Data is stored under `PAIR_*/computational_data/<label>/fpocket/` for transparency.
- The script still accepts user-provided annotations (`--computational-dir`). All sources are merged.

---

## 6. Modes & performance knobs

- `--input-mode split` → supply each chain separately. The script auto-creates body restraints to keep chain groups rigid.
- `--input-mode multichain` → each partner is a multichain PDB.
- `--ncores` influences parallelism in HADDOCK (default 10).
- `--sampling` adjusts rigid-body sampling (versions: V1 default 10000, others 4000 unless overridden). Lowering it speeds up first stage.
- `--skip-flexref` (if you add it) or editing `[flexref]/[mdref]/[emref]` blocks can shorten runs, but by default all refinement stages run sequentially.
- Even when higher versions fall back to blind docking (no restraints available), `ranair` switches on automatically.

---

## 7. Typical workflows

### 7.1 All-vs-all, blind docking
1. Generate configs:
   ```bash
   ./script/a.sh --auto-partners protein_test --v 1 --use-fpocket
   ```
2. Run:
   ```bash
   ls result/<project>/PAIR_*/haddock3.cfg | parallel -j 5 "cd \$(dirname {}) && haddock3 \$(basename {})"
   ```

### 7.2 Guided docking with external annotations
```bash
./script/a.sh --auto-partners protein_test --computational-dir annotations/ --experimental-dir exp_data/
# Then run with parallel command
```

### 7.3 Split chains with body locks
```bash
./script/a.sh \
  --input-mode split \
  --partner P1CLL=split/P1CLL_A.pdb,split/P1CLL_B.pdb \
  --partner P1ZG4=split/P1ZG4_A.pdb,split/P1ZG4_B.pdb \
  --group-bodies P1CLL=reference/P1CLL_all.pdb \
  --group-bodies P1ZG4=reference/P1ZG4_all.pdb \
  --use-fpocket
# Then run with parallel command
```

### 7.4 Dry-run (inspect configs)
```bash
./script/a.sh --auto-partners protein_test --use-fpocket --dry-run
```
Check `result/<project>/PAIR_*/haddock3.cfg`, computational data, etc. Since execution is decoupled, `--dry-run` effectively just skips the final "ready to run" log message or specific setup checks if added back.

---

## 8. Notes & best practices

- Always validate that `haddock3` and `fpocket` are available in the activated shell.
- If you want the traditional CAPRI steps between every stage, you can extend `create_config` to insert `[caprieval]` etc.; the current setup relies on `haddock3-analyse` after completion.
- For comparisons or debug, consider running baseline (multichain PDB) + split without body locks + split with locks to check the effect on docking quality.
- Summary stats (success/failure counts, guided vs blind counts, number of auto restraints) are printed at the end.
- `computational_data/<label>/fpocket/raw/` retains original fpocket output for reference.
- You can still supply manual restraints (they take priority when present).

---

## 8. Troubleshooting

| Symptom | Possible cause & fix |
|---------|----------------------|
| `Missing required command 'haddock3'` | Activate the HADDOCK3 environment (`conda activate haddock3`). |
| `fpocket output missing summary` warning | Old runs; rerun after the parser update, or check fpocket executable. |
| `haddock3 --setup` fails | Inspect `PAIR_*/setup.log` for parameter issues (e.g. bad `cmrest` values). Note: Setup is now run during the parallel execution phase. |
| Docking falls back to blind | Logs show “No restraints”; script will set `ranair=true` automatically. Consider providing more data. |
| Long runtimes | Reduce `--sampling`, tweak `select` in `[seletop]`, disable `mdref`/`emref` (requires editing config writer). |

---

## 9. References

---

## 9. References
- [HADDOCK3 user manual](https://www.bonvinlab.org/haddock3-user-manual/)
- [fpocket GitHub](https://github.com/Discngine/fpocket)
- Additional script logic and pipeline notes live in `README.md`.

---

Happy docking!
