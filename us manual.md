# HADDOCK3 Automation Script — User Manual

## 1. Purpose & Scope
The Bash script at `script/a.sh` automates multi-stage HADDOCK3 protein–protein docking. It builds clean PDB inputs, organises restraints, writes run-ready configuration files, and optionally executes `haddock3`. It supports the three workflow “versions” described by Bonvin Lab (blind, computationally guided, experimentally driven) and covers both split-chain and multi-chain input conventions. This manual explains every option, the processing pipeline, and provides examples for both single-pair and batch experiments.

---

## 2. Quick Start Checklist
1. **Install dependencies**: Python 3, HADDOCK3 (`haddock3`, `haddock3-restraints`).
2. **Prepare input structures**: place PDB files (split or multi-chain) in a working directory. Ensure the script can access relevant annotation folders (`ambig.tbl`, `unambig.tbl`, CSV/JSON predictions, etc.).
3. **Run a dry rehearsal**:
   ```bash
   ./script/a.sh --auto-partners protein_test --v 3 \
                 --computational-dir data/predict --experimental-dir data/exp \
                 --project demo_v3 --dry-run
   ```
   Inspect the generated run folders under `result/demo_v3/`.
4. **Launch production run**: re-run without `--dry-run`.

---

## 3. Input Concepts
### 3.1 Partner definitions
- Each *partner* corresponds to one biological entity (protein, complex) participating in docking.
- A partner can be supplied as:
  - A single multi-chain PDB file (`proteinA_multichain.pdb`), or
  - A set of per-chain PDBs (split mode), or
  - A directory / CSV list enumerating PDBs.
- Labels are arbitrary but must be unique (e.g. `A`, `B`, `Kinase`, `Ligand`).

### 3.2 Split vs multichain modes
- `--input-mode multichain` (default): treat each provided file as a multi-chain molecule. The script preserves original chain IDs unless collisions force renumbering.
- `--input-mode split`: treat each file as an individual chain. The script renumbers and assigns unique chain IDs per partner. Split mode **requires** body-restraint manifests (see §5.3) to keep chains grouped.

### 3.3 Automatic partner discovery
`--auto-partners <dir>` scans `<dir>` for `*.pdb` and registers each file as a partner (label = uppercase basename). You can still mix auto-discovered partners with explicit `--partner` declarations.

---

## 4. Command-line Parameters
The script accepts repeatable `--partner` arguments and optional modifiers. You can also use `./script/a.sh --help` to print this summary.

| Category | Flag | Description |
| --- | --- | --- |
| **Partner setup** | `--partner <label>=<spec>` | Register a partner. `<spec>` can be a file, directory, comma-separated file list, or `.lst` file with one path per line. |
| | `--auto-partners <dir>` | Auto-discover partners by scanning `<dir>` for `*.pdb`. Labels derive from filenames (uppercased, non-alphanumerics → `_`). |
| | `--input-mode {split,multichain}` | Choose how partner files are interpreted (see §3). Default: `multichain`. |
| | `--group-bodies <label>=<path>` | Provide a reference PDB or JSON mapping for partner `<label>` to generate body restraints in split mode. Repeat per partner. |
| | `--group-bodies <json>` | If `<json>` omits `label=`, it is parsed as a manifest mapping labels to source files and optional chain lists. |
| **Workflow versioning** | `--v`, `--version {1|2|3}` | Select pipeline: V1 blind, V2 computational, V3 experimental priority. Default 3. |
| | `--abinitio` | Convenience flag for version 1 with high sampling and random AIR seeding (`ranair=true`). |
| | `--ranair` | Force `[rigidbody] ranair = true` regardless of version. Useful for sparse/blind docking. |
| **Restraints** | `--ambig <tbl[,tbl…]>` | Manual ambiguous restraints copied into each run directory. |
| | `--unambig <tbl[,tbl…]>` | Manual unambiguous restraints (e.g. experimental, body locks). |
| | `--computational-dir <dir>` | Root directory containing computational predictions (CSV/TSV/JSON). Used in V2/V3. |
| | `--experimental-dir <dir>` | Root directory with experimental annotations (JSON/CSV). Prioritised in V3. |
| | `--confidence-threshold <float>` | Minimum average confidence (0–1) to accept computational predictions. Default 0.6. |
| | `--max-active <int>` | Cap active residues per partner when building ambiguous restraints (default 40). |
| | `--max-pairs <int>` | Maximum ambiguous restraints generated per pair (default 4000). |
| **Pair selection** | `--pair <label1,label2>` | Restrict runs to specific partner pairs. Repeatable. |
| | `--pairs-file <file>` | Provide text/CSV file (one pair per line) selecting pairings. Supports `label1,label2`, `label1 label2`, or `label1|label2`. |
| **Execution** | `--out <dir>` | Root output directory (`$PWD/result` by default). |
| | `--project <name>` | Folder created inside `--out`. Defaults to `dock_<timestamp>`. |
| | `--ncores <int>` | Number of cores to pass to HADDOCK3. Default 4 (set to host capacity in script header). |
| | `--sampling <int>` | Fixed rigid-body sampling for all runs. If omitted, V1 defaults to 10000, V2/V3 to 4000. |
| | `--reference <pdb>` | Native complex for CAPRI evaluation (`[caprieval]` stage). Copied into each pair directory. |
| | `--dry-run` / `--no-run` | Prepare inputs and configs but skip `haddock3 --setup` and execution. |
| | `--help` | Print usage text.

You may combine multiple flags; any unspecified option falls back to the values defined at the top of `script/a.sh`.

---

## 5. Processing Pipeline
### 5.1 Pre-processing
For each partner:
1. All input PDB lines are normalised to 80-column format; non-standard residues are skipped.
2. Chain IDs are reassigned (avoid collisions across partners). Residues are renumbered sequentially if `--input-mode split` or when the original numbering is inconsistent.
3. TER/END records are enforced, SEGIDs filled, atom serials re-indexed.
4. Chain-specific PDBs are emitted as `prepared/<label>/<label>_chain_<ID>.pdb` alongside a combined file and JSON mapping (`<label>_mapping.json`).
5. If split mode yields multiple chains, `haddock3-restraints restrain_bodies` autogenerates “body lock” restraints to hold chains together.
6. Optional `--group-bodies` references are copied to `prepared/<label>/group_source/` and converted into additional body restraints.

### 5.2 Pair scheduling
- If `--pair/--pairs-file` are absent, the script builds an all-vs-all list of unique partner pairs (combinations). With five partners you obtain ten separate runs.
- Pair selection honours order in `--pair` arguments; duplicates are ignored.
- For each pair, a dedicated run directory is created:
  ```
  <out>/<project>/PAIR_<Label1>_vs_<Label2>/
    ├── input/                # Partner-specific copies of cleaned PDBs
    ├── restraints/           # Manual, body, and auto-generated restraints
    ├── meta/                 # Residue mappings for both partners
    ├── haddock3.cfg          # Auto-written configuration
    ├── setup.log, haddock3.log (if executed)
    └── reference.pdb (optional)
  ```

### 5.3 Restraint sources and precedence
1. **Manual**: `--ambig`/`--unambig` files are copied verbatim into each pair directory.
2. **Body locks**: Automatically generated from split inputs and optional group manifests. These accumulate in the unambiguous list.
3. **Automatic predictions** (V2/V3):
   - Computational annotations become ambiguous restraints if above the confidence cutoff.
   - Experimental annotations (V3) populate unambiguous restraints; computational residues supplement ambiguous lists.
   - Limits obey `--max-active` (per partner) and `--max-pairs` (total restraints).
4. The final `[flexref]/[mdref]/[emref]` sections receive comma-separated lists of ambiguous/unambiguous tables.

### 5.4 Version-specific behaviour
- **Version 1 (blind)**: emphasises sampling (`sampling ≥ 10000`) and optionally random AIRs via `ranair=true`. No mandatory restraints.
- **Version 2 (computational)**: auto-generates ambiguous restraints from predictions. `randremoval=true` is set when loading ambiguous tables to avoid overfitting.
- **Version 3 (hybrid)**: prioritises experimental inputs as unambiguous restraints, while computational data remains ambiguous.

---

## 6. Execution Modes
### 6.1 Dry-run vs full execution
- **Dry-run (`--dry-run`)**: Generates cleaned inputs, restraints, and configs, then stops. Summary still reports planned pairs (counted as “skipped”).
- **Full run**: After config creation, the script runs `haddock3 --setup` followed by `haddock3`. Failures are logged per pair without halting subsequent pairs.

### 6.2 Parallelism
Runs currently execute sequentially to simplify logging. If you need parallel execution, launch multiple instances manually or extend the script with job control wrappers (e.g. GNU parallel).

---

## 7. Practical Scenarios
### 7.1 All-vs-all docking from a directory
```
./script/a.sh --auto-partners protein_test --v 1 --abinitio \
              --input-mode multichain --out result --project screen_v1
```
- Discovers every `*.pdb` under `protein_test`.
- Schedules all unique pairs using blind presets (10k sampling, random AIRs).
- Results appear under `result/screen_v1/PAIR_*`.

### 7.2 Split-chain workflow with explicit body groups
```
./script/a.sh \
  --input-mode split \
  --partner A=chains/proteinA --partner B=chains/proteinB \
  --group-bodies A=manifests/proteinA.json \
  --group-bodies B=manifests/proteinB.json \
  --v 3 \
  --computational-dir annotations/predict \
  --experimental-dir annotations/exp \
  --ambig restraints/manual_calc.tbl \
  --project split_v3
```
- JSON manifests list chain files belonging to each protein and (optionally) desired ordering.
- Body restraints are generated both from the split inputs and the manifest reference.
- Version 3 loads experimental tables as unambiguous and supplements with computational data.

### 7.3 Curated pair list
Create `pairs.csv`:
```
A,B
A,C
C,D
```
Then run:
```
./script/a.sh --auto-partners protein_test --pairs-file pairs.csv --v 2 \
              --computational-dir predictions --project guided_v2
```
Only specified combinations are processed.

### 7.4 Blind docking fallback when predictions fail
If automatic restraint generation finds no trustworthy residues (e.g. low confidence), the script logs a warning and falls back to blind sampling for that pair. Increase `--confidence-threshold`, supply manual restraints, or provide better annotations.

---

## 8. Output & Logs
- `prepared/` holds partner-specific cleaned PDBs and mapping metadata.
- Each pair directory contains:
  - `input/*.pdb`: ready-to-use molecules referenced by `haddock3.cfg`.
  - `restraints/*.tbl`: merged manual, body, and auto-generated tables.
  - `meta/<label>_mapping.json`: residue mapping dictionary used during restraint projection.
  - `haddock3.cfg`: configuration referencing relative paths within the pair directory.
  - `setup.log`, `haddock3.log`: output from HADDOCK3 setup and run phases.
- A final summary (printed to stdout) reports success/failure counts, guided vs blind runs, and the number of automatically generated restraint tables.

---

## 9. Troubleshooting
| Symptom | Likely Cause | Remedy |
| --- | --- | --- |
| `Duplicate partner label` error | Two inputs resolved to same label | Use explicit `--partner Label=...` or rename files. |
| `split mode requires --group-bodies` | Body grouping not supplied | Provide `--group-bodies` per partner or a JSON manifest. |
| `Automatic restraint generation failed` warnings | Missing/invalid annotations | Check paths, reduce confidence cutoff, or switch to version 1. |
| `haddock3` command not found | HADDOCK3 not on PATH | Activate HADDOCK3 environment or install per Bonvin Lab instructions. |
| Generated config references wrong paths | Running script from unexpected directory | Use absolute paths (`--out`, `--project`, `--partner` specs) consistently. |

---

## 10. Best Practices & Notes
- Keep PDBs as clean as possible; although the script sanitises columns, severe formatting issues may still propagate.
- For multi-chain partners with natural breaks, consider providing a `--group-bodies` manifest even in multichain mode to enforce realistic rigidity during refinement.
- When running large batches, prefer `--dry-run` first to verify scheduling and restraint assembly.
- Store computational/experimental annotations in structured directories (one JSON/CSV per tool) to simplify re-use.
- Remove large `restraints` or `result` folders periodically to reclaim disk space.

---

## 11. Change Log Highlights (relative to previous script)
- Supports **multiple partners** with automatic all-vs-all pairing or custom manifests.
- Adds `--auto-partners`, `--pair`, `--pairs-file` for batch scheduling.
- Generates body restraints both from split PDBs and optional external references.
- Writes summary statistics even during dry-runs.

---

For questions or suggestions, update the repository README or reach out through project communication channels.
