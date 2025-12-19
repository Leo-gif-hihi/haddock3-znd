#!/usr/bin/env bash
# ===========================================================================
# HADDOCK3 Automated Docking Workflow Script
# ===========================================================================
# This script automates the setup and execution of HADDOCK3 docking runs.
# It handles:
#   1. Input preparation (PDB cleaning, chain renaming, merging).
#   2. Partner discovery and pairing (manual or automatic).
#   3. Restraint generation (from experimental or computational data).
#   4. Configuration generation (haddock3.cfg).
#   5. Execution of HADDOCK3 (or preparation for parallel execution).
#
# Usage:
#   ./a.sh --partner A=proteinA.pdb --partner B=proteinB.pdb [options]
#
# See 'usage()' function or run with --help for more details.
# ===========================================================================

set -euo pipefail
IFS=$'\n\t'

SCRIPT_NAME=$(basename "$0")
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# ---------------------------------------------------------------------------
# Global Configuration & Defaults
# ---------------------------------------------------------------------------
FORCE_CHAINS_TOGETHER=true

# Partner & Grouping Data Structures
declare -a PARTNER_ORDER=()          # Ordered list of partner labels
declare -A PARTNER_SPECS=()          # Input specifications (file paths/patterns)
declare -A PARTNER_FILES=()          # Resolved absolute file paths
declare -A PARTNER_COMBINED=()       # Path to combined/cleaned PDB
declare -A PARTNER_MAPPING=()        # Path to JSON mapping (original -> new chains)
declare -A PARTNER_BODY_TBL=()       # Path to rigid body restraint table
declare -A PARTNER_SPLITFILES=()     # List of split chain files (if applicable)
declare -A GROUP_BODY_SOURCE=()      # Source PDB for body grouping (split mode)
declare -A GROUP_BODY_CHAINLIST=()   # Chain order for grouping
declare -a GROUP_BODY_JSON_FILES=()  # JSON manifests for grouping
declare -A PARTNER_GROUP_BODY=()     # Generated body restraints from group source

# Restraints & Pairs
declare -a AMBIG_MANUAL=()           # Manual ambiguous restraint files
declare -a UNAMBIG_MANUAL=()         # Manual unambiguous restraint files
declare -a PAIR_OVERRIDES=()         # Specific pairs requested via CLI
declare -A PAIR_SET=()               # Set of unique pairs to avoid duplicates
declare -a PAIR_LIST=()              # Final list of pairs to dock
PAIR_FILE=""                         # File containing list of pairs

# Directories & Paths
AUTO_PARTNER_DIR=""
COMPUTATIONAL_DIR=""
EXPERIMENTAL_DIR=""
OUTPUT_ROOT="$PWD/result"
PROJECT_NAME=""

# Run Parameters
NCORES="auto"                   # "auto" = detect system cores - 2, or user-specified integer
SAMPLING_OVERRIDE=""
CONF_THRESHOLD=0.6
MAX_ACTIVE_RESIDUES=40
MAX_RESTRAINT_PAIRS=4000
ABINITIO_PRESET=false
FORCE_RANAIR=false
DRY_RUN=false
SKIP_RUN=false
REFERENCE_PDB=""
VERSION=3  # Automatic mode

# Run Mode Settings
RUN_MODE="sequential"           # "sequential" (default, for local/weak computers) or "parallel" (for powerful/cloud)
PARALLEL_JOBS=5                  # Number of parallel jobs when using parallel mode
EXECUTE_JOBS=true                # If true, automatically execute HADDOCK3 jobs after generation (default: true)

# ARCTIC3D Settings
USE_ARCTIC3D=false
ARCTIC3D_PROB_THRESHOLD=0.4
ARCTIC3D_AVAILABLE="unknown"
ARCTIC3D_BLAST_DB=""
ARCTIC3D_POSITION_TOLERANCE=1.0  # Angstrom tolerance for CA position matching

# ARCTIC3D Data Structures (populated before cleaning)
declare -A PARTNER_ORIGINAL_FILES=()      # Original PDB files (before cleaning)
declare -A PARTNER_CHAIN_UNIPROTS=()      # Chain -> Uniprot ID mapping (JSON string)
declare -A PARTNER_PDB_IDS=()             # PDB ID extracted from original files
declare -A PARTNER_ARCTIC3D_MANIFEST=()   # Path to ARCTIC3D manifest for each partner

# PDB Cleaning Options
REMOVE_HETATM=false

# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------

# Print usage information
usage() {
    cat <<'USAGE'
Usage: a.sh --partner A=<path> --partner B=<path> [options]

Required:
  --partner <label>=<spec>    Define docking partner input. Provide at least two partners.
                             <spec> can be a PDB file, a directory with *.pdb, a comma-separated
                             list of PDB files, or a .lst file with one path per line.

Modes & presets:
  --keep-chains-separate          Don't force chains from the same PDB together during docking.
  --group-bodies <label>=<pdb|json>  Mapping for body restraints. Repeat per partner or provide a JSON file.
  --auto-partners <dir>       Discover partners automatically from PDBs inside <dir>.
  --abinitio                   Blind docking mode: high sampling, ranair, no restraint processing.
  --ranair                     Force [rigidbody] ranair = true.

Pair selection (optional):
  --pair <label1,label2>       Restrict docking to the specified labelled pair. Repeatable.
  --pairs-file <file>          Text/CSV file with one pair per line (label1,label2).

Restraints & data sources:
  --ambig <tbl[,tbl..]>        Manual ambiguous restraint files (copied into the run directory).
  --unambig <tbl[,tbl..]>      Manual unambiguous restraint files.
  --computational-dir <dir>    Directory tree with computational annotations (for V2/V3).
  --experimental-dir <dir>     Directory tree with experimental annotations (for V3 priority).
  --use-arctic3d               Run ARCTIC3D on prepared structures to generate interface restraints.
  --arctic3d-prob <float>      Probability threshold for ARCTIC3D restraints (default: 0.4).
  --arctic3d-blast-db <path>   Path to local BLAST database for ARCTIC3D (optional).
  --confidence-threshold <f>   Confidence cutoff for computational data (default 0.6).
  --max-active <int>           Max active residues per partner when building restraints (default 40).
  --max-pairs <int>            Cap number of ambiguous restraints (default 4000).

Execution & I/O:
  --out <dir>                  Root directory for generated runs (default: $PWD/result).
  --project <name>             Name of the run folder inside --out (default: auto timestamp).
  --ncores <int|auto>          Number of cores passed to HADDOCK3 (default: auto).
                               'auto' = system cores - 2 (minimum 1). Specify integer to override.
  --sampling <int>             Override rigid-body sampling value.
  --reference <pdb>            Native complex for CAPRI evaluation.
  --dry-run                    Prepare config/artifacts but skip haddock3 execution.
  --remove-hetatm              Remove all HETATM records (default: keep all valid molecules).
  --run-mode <mode>            Execution mode: 'sequential' (default, for local/weak computers)
                               or 'parallel' (for powerful computers/cloud). 
  --parallel-jobs <int>        Number of parallel jobs when using --run-mode parallel (default: 5).
  --help                       Print this message.

Examples:
  # Multi-chain docking using experimental + computational restraints (Version 3)
  ./a.sh --partner A=proteinA.pdb --partner B=proteinB.pdb --v 3 \
        --computational-dir data/predict --experimental-dir data/exp --project demo_v3

  # Docking with chains moving independently (no body restraints)
  ./a.sh --partner A=proteinA.pdb --partner B=proteinB.pdb \
        --keep-chains-separate --ambig manual_air.tbl

See README for details on manifest formats accepted by --group-bodies.
USAGE
}

# Logging helpers
log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
warn() { log "WARN: $*"; }
error() { log "ERROR: $*"; }
die() { error "$*"; exit 1; }

# Detect the number of CPU cores on the system
get_system_cores() {
    local cores=1
    if [[ -f /proc/cpuinfo ]]; then
        # Linux: count processor entries
        cores=$(grep -c '^processor' /proc/cpuinfo 2>/dev/null || echo 1)
    elif command -v nproc >/dev/null 2>&1; then
        # Linux fallback: use nproc
        cores=$(nproc 2>/dev/null || echo 1)
    elif command -v sysctl >/dev/null 2>&1; then
        # macOS: use sysctl
        cores=$(sysctl -n hw.ncpu 2>/dev/null || echo 1)
    elif [[ -n "${NUMBER_OF_PROCESSORS:-}" ]]; then
        # Windows (Git Bash, WSL passthrough)
        cores="$NUMBER_OF_PROCESSORS"
    fi
    # Ensure we return at least 1
    [[ "$cores" -gt 0 ]] 2>/dev/null || cores=1
    echo "$cores"
}

# Resolve NCORES: if "auto", set to system cores - 2 (minimum 1)
resolve_ncores() {
    local ncores_setting="$1"
    if [[ "$ncores_setting" == "auto" ]]; then
        local system_cores
        system_cores=$(get_system_cores)
        local resolved=$(( system_cores - 2 ))
        # Ensure minimum of 1 core
        [[ "$resolved" -lt 1 ]] && resolved=1
        log "Auto-detected $system_cores CPU cores, using $resolved cores for HADDOCK3"
        echo "$resolved"
    else
        # User specified explicit value
        echo "$ncores_setting"
    fi
}

# Check if a command exists in PATH
require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Missing required command '$1'"
}

# Resolve absolute path using Python
abs_path() {
    python3 - "$1" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
}

# Split a comma-separated string into lines
split_csv() {
    local list="$1"
    local IFS=','
    read -r -a _tmp <<< "$list"
    printf '%s\n' "${_tmp[@]}"
}

# Select a pool of chain IDs based on index to avoid collisions
select_chain_pool() {
    local index="$1"
    local base="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    local length=${#base}
    local offset=$(( index % length ))
    printf '%s%s\n' "${base:offset}" "${base:0:offset}"
}

# Append a value to an array only if it's not already present
append_unique() {
    local -n _arr=$1
    local value="$2"
    local existing
    for existing in "${_arr[@]}"; do
        if [[ "$existing" == "$value" ]]; then
            return 0
        fi
    done
    _arr+=("$value")
}

# Generate a canonical key for a pair of partners (sorted)
canonical_pair_key() {
    local a="$1"
    local b="$2"
    if [[ "$a" < "$b" ]]; then
        printf '%s|%s\n' "$a" "$b"
    else
        printf '%s|%s\n' "$b" "$a"
    fi
}

# ---------------------------------------------------------------------------
# Chain Count Detection for Pair Validation
# ---------------------------------------------------------------------------
MAX_CHAIN_LIMIT=36
declare -a SKIPPED_PAIRS_INFO=()  # Stores "lhs|rhs|total_chains|reason" for skipped pairs

# Count chains from a partner's mapping JSON file
count_partner_chains() {
    local mapping_json="$1"
    if [[ ! -f "$mapping_json" ]]; then
        echo "0"
        return
    fi
    python3 -c "import json; print(len(json.load(open('$mapping_json')).get('chain_mapping', {})))"
}

# Check if a pair exceeds the maximum chain limit
# Returns 0 (true) if pair should be SKIPPED, 1 (false) if OK to process
check_pair_chain_limit() {
    local lhs="$1"
    local rhs="$2"
    local lhs_mapping="${PARTNER_MAPPING[$lhs]:-}"
    local rhs_mapping="${PARTNER_MAPPING[$rhs]:-}"
    
    local lhs_chains=0
    local rhs_chains=0
    
    if [[ -n "$lhs_mapping" && -f "$lhs_mapping" ]]; then
        lhs_chains=$(count_partner_chains "$lhs_mapping")
    fi
    
    if [[ -n "$rhs_mapping" && -f "$rhs_mapping" ]]; then
        rhs_chains=$(count_partner_chains "$rhs_mapping")
    fi
    
    local total_chains=$((lhs_chains + rhs_chains))
    
    if (( total_chains >= MAX_CHAIN_LIMIT )); then
        # Return the total chain count for logging
        echo "$total_chains"
        return 0  # Should skip
    fi
    
    echo "$total_chains"
    return 1  # OK to process
}

# Add a pair to the list if it hasn't been added yet
add_pair_unique() {
    local lhs="$1"
    local rhs="$2"

    [[ -z "$lhs" || -z "$rhs" ]] && return 1
    [[ "$lhs" == "$rhs" ]] && die "Pair definition requires two distinct partners (got '$lhs' twice)."

    local canonical
    canonical=$(canonical_pair_key "$lhs" "$rhs")
    if [[ -z "${PAIR_SET[$canonical]:-}" ]]; then
        PAIR_SET["$canonical"]=1
        PAIR_LIST+=("$lhs|$rhs")
    fi
    return 0
}

# Parse a JSON file defining body groups for split mode
parse_group_json() {
    local json_path="$1"
    [[ -f "$json_path" ]] || die "Group manifest not found: $json_path"
    mapfile -t _entries < <(python3 - "$json_path" <<'PY'
import json
import os
import sys
path = sys.argv[1]
with open(path, 'r', encoding='utf-8') as handle:
    data = json.load(handle)
if isinstance(data, dict) and 'groups' in data and isinstance(data['groups'], dict):
    data = data['groups']
if not isinstance(data, dict):
    raise SystemExit('Expected mapping at top level of %s' % path)
for key, value in data.items():
    if isinstance(value, str):
        source = value
        chains = []
    elif isinstance(value, dict):
        source = value.get('source') or value.get('pdb') or value.get('file') or ''
        chains = value.get('chains') or value.get('files') or []
    else:
        source = ''
        chains = []
    if isinstance(chains, str):
        chains = [chain.strip() for chain in chains.split(',') if chain.strip()]
    elif isinstance(chains, list):
        chains = [str(chain).strip() for chain in chains if str(chain).strip()]
    else:
        chains = []
    print(f"{key}|{source}|{','.join(chains) if chains else '-'}")
PY
    )
    local line label source chains
    for line in "${_entries[@]}"; do
        [[ -z "$line" ]] && continue
        label=${line%%|*}
        source=${line#*|}
        source=${source%%|*}
        chains=${line##*|}
        if [[ -n "$source" ]]; then
            GROUP_BODY_SOURCE["$label"]="$source"
        fi
        if [[ -n "$chains" && "$chains" != '-' ]]; then
            GROUP_BODY_CHAINLIST["$label"]="$chains"
        fi
    done
}

# Resolve partner specifications (files, directories, lists) to absolute paths
resolve_partner_spec() {
    local label="$1"
    local spec="$2"
    local -a files=()

    if [[ "$spec" == *,* ]]; then
        mapfile -t files < <(split_csv "$spec")
    elif [[ -d "$spec" ]]; then
        mapfile -t files < <(find "$spec" -maxdepth 1 -type f -name '*.pdb' | sort)
    elif [[ -f "$spec" && "$spec" == *.lst ]]; then
        mapfile -t files < <(grep -v '^#' "$spec" | sed '/^\s*$/d')
    elif [[ -f "$spec" ]]; then
        files=("$spec")
    else
        die "Partner '$label' spec does not resolve to files: $spec"
    fi

    (( ${#files[@]} == 0 )) && die "Partner '$label' has no PDB files"

    local idx
    for idx in "${!files[@]}"; do
        files[$idx]=$(abs_path "${files[$idx]}")
        [[ -f "${files[$idx]}" ]] || die "File not found for partner '$label': ${files[$idx]}"
    done

    PARTNER_FILES["$label"]=$(printf '%s\n' "${files[@]}")
}

# Reorder partner files based on group manifest if present
resolve_group_chain_order() {
    local label="$1"
    local -a files=()
    mapfile -t files < <(printf '%s\n' "${PARTNER_FILES[$label]}")
    local manifest="${GROUP_BODY_CHAINLIST[$label]:-}"
    [[ -z "$manifest" ]] && return 0
    local IFS=','
    read -r -a desired <<< "$manifest"
    (( ${#desired[@]} == 0 )) && return 0

    local -a reordered=()
    local name file found
    for name in "${desired[@]}"; do
        found=""
        for file in "${files[@]}"; do
            if [[ "$(basename "$file")" == "$name" ]]; then
                found="$file"
                break
            fi
        done
        [[ -z "$found" ]] && die "Group manifest lists '$name' for '$label' but file not provided"
        reordered+=("$found")
    done
    if (( ${#reordered[@]} != ${#files[@]} )); then
        warn "Group manifest for '$label' does not cover all files; keeping unmatched afterwards."
        local extras=()
        for file in "${files[@]}"; do
            local base="$(basename "$file")"
            local present=false
            for name in "${desired[@]}"; do
                [[ "$base" == "$name" ]] && { present=true; break; }
            done
            $present || extras+=("$file")
        done
        reordered+=("${extras[@]}")
    fi
    PARTNER_FILES["$label"]=$(printf '%s\n' "${reordered[@]}")
}

# ---------------------------------------------------------------------------
# Dynamic Rechaining Helpers
# ---------------------------------------------------------------------------

rechain_partner() {
    local json_file="$1"
    local start_char="$2"
    local arctic_dir="$3"
    shift 3
    local pdb_files=("$@")
    
    python3 - "$json_file" "$start_char" "$arctic_dir" "${pdb_files[@]}" <<'PY'
import sys
import json
import os
import shutil

json_path = sys.argv[1]
start_char = sys.argv[2]
arctic_dir = sys.argv[3]
pdb_files = sys.argv[4:]

chain_pool = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
try:
    start_idx = chain_pool.index(start_char)
except ValueError:
    start_idx = 0

with open(json_path, 'r') as f:
    mapping = json.load(f)

chain_map = mapping.get("chain_mapping", {})
# Sort by current new_chain to ensure stability
sorted_chains = sorted(chain_map.items(), key=lambda x: x[1]["new_chain"])

remap_dict = {}
current_idx = start_idx

for original_chain, info in sorted_chains:
    old_chain = info["new_chain"]
    if current_idx >= len(chain_pool):
        raise ValueError("Ran out of chain identifiers during rechaining")
    new_chain = chain_pool[current_idx]
    info["new_chain"] = new_chain
    info["segid"] = new_chain
    remap_dict[old_chain] = new_chain
    current_idx += 1

# Update JSON
with open(json_path, 'w') as f:
    json.dump(mapping, f, indent=2)

# Update PDBs
for pdb_path in pdb_files:
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            chain = line[21]
            if chain in remap_dict:
                new_chain = remap_dict[chain]
                # Update chain ID (col 22) and segid (col 73-76)
                line = line[:21] + new_chain + line[22:72] + new_chain.rjust(4) + line[76:]
        new_lines.append(line)
    with open(pdb_path, 'w') as f:
        f.writelines(new_lines)

# Rename Arctic Folders
if arctic_dir and os.path.isdir(arctic_dir):
    moves = []
    for old_c, new_c in remap_dict.items():
        old_p = os.path.join(arctic_dir, f"chain_{old_c}")
        new_p = os.path.join(arctic_dir, f"chain_{new_c}")
        if os.path.exists(old_p):
            moves.append((old_p, new_p))
            
    # Use temp moves to avoid collisions (e.g. A->B when B exists)
    temp_moves = []
    final_moves = []
    for i, (old_p, new_p) in enumerate(moves):
        temp_p = f"{old_p}_tmp_{i}"
        temp_moves.append((old_p, temp_p))
        final_moves.append((temp_p, new_p))
        
    for src, dst in temp_moves:
        os.rename(src, dst)
    for src, dst in final_moves:
        if os.path.exists(dst):
            shutil.rmtree(dst)
        os.rename(src, dst)

    # Update manifest
    manifest_path = os.path.join(arctic_dir, "chain_manifest.json")
    if os.path.exists(manifest_path):
        with open(manifest_path, 'r') as f:
            manifest = json.load(f)
        new_chains = {}
        for old_c, new_c in remap_dict.items():
            if old_c in manifest.get("chains", {}):
                new_chains[new_c] = manifest["chains"][old_c]
        manifest["chains"] = new_chains
        with open(manifest_path, 'w') as f:
            json.dump(manifest, f, indent=2)

# Print mapping for shell usage: A:C,B:D
pairs = [f"{k}:{v}" for k, v in remap_dict.items()]
print(",".join(pairs))
PY
}

remap_tbl_chain() {
    local tbl_file="$1"
    local chain_map="$2" # e.g. "A:C,B:D"
    
    python3 - "$tbl_file" "$chain_map" <<'PY'
import sys
import re

tbl_path = sys.argv[1]
chain_map_str = sys.argv[2]

if not chain_map_str:
    sys.exit(0)

mapping = {}
for pair in chain_map_str.split(','):
    if ':' in pair:
        k, v = pair.split(':')
        mapping[k] = v

with open(tbl_path, 'r') as f:
    content = f.read()

# Regex to find segid "X" or chain "X" in TBL
# HADDOCK TBL format usually: assign (segid A and resid 10) ...
# We need to be careful not to replace other things.
# Usually segid is used.

def replace_chain(match):
    chain = match.group(1)
    return f'segid {mapping.get(chain, chain)}'

# Replace 'segid X'
pattern = re.compile(r'segid\s+([A-Za-z0-9])')
new_content = pattern.sub(replace_chain, content)

with open(tbl_path, 'w') as f:
    f.write(new_content)
PY
}

# ---------------------------------------------------------------------------
# Source ARCTIC3D Functions
# ---------------------------------------------------------------------------
# Load ARCTIC3D-related functions from separate file for better code organization
ARCTIC3D_FUNCTIONS_FILE="$SCRIPT_DIR/arctic3d_functions.sh"
if [[ -f "$ARCTIC3D_FUNCTIONS_FILE" ]]; then
    source "$ARCTIC3D_FUNCTIONS_FILE"
elif [[ "$USE_ARCTIC3D" == true ]]; then
    warn "ARCTIC3D functions file not found: $ARCTIC3D_FUNCTIONS_FILE"
    warn "ARCTIC3D features will not be available."
    USE_ARCTIC3D=false
fi

# ---------------------------------------------------------------------------
# Core Logic Functions
# ---------------------------------------------------------------------------

# Prepare a partner structure: clean, rename chains, merge if needed
run_prepare_partner() {
    local label="$1"
    local index="$2"
    local -a files=()
    mapfile -t files < <(printf '%s\n' "${PARTNER_FILES[$label]}")
    local prep_dir="$3"
    mkdir -p "$prep_dir"

    local rename_mode="rename"

    local chain_pool="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    # local pool_offset=$(( index * 8 ))
    # local pool="${chain_pool:pool_offset}${chain_pool:0:pool_offset}"
    local pool="$chain_pool"

    mapfile -t PREP_INFO < <(
        python3 - "$rename_mode" "$prep_dir" "$label" "$pool" "$FORCE_CHAINS_TOGETHER" "$REMOVE_HETATM" "${files[@]}" <<'PY'
import os
import sys
import json
import shutil
import subprocess
from collections import OrderedDict, defaultdict

rename_mode = sys.argv[1]
output_dir = sys.argv[2]
label = sys.argv[3]
chain_pool = sys.argv[4]
force_together = sys.argv[5].lower() == 'true'
remove_hetatm = sys.argv[6].lower() == 'true'
input_files = sys.argv[7:]

if rename_mode not in {"rename", "preserve"}:
    raise SystemExit(f"Unsupported rename mode: {rename_mode}")

if not input_files:
    raise SystemExit("No input files provided")

STANDARD_RESIDUES = {
    # Standard amino acids
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "MSE",
    # Protonation states and variants (already includes ASH, GLH from modified amino acids)
    "HID", "HIE", "HIP", "ASH", "GLH", "LYN",
    # Carbohydrates
    "A2G", "ABE", "BDP", "BGC", "BMA", "FCA", "FCB", "FUC",
    "FUL", "GAL", "GLB", "GLA", "GLC", "GXL", "MAG", "MAN",
    "MMA", "NAG", "NDG", "NGA", "RAM", "SIA", "SIB", "XYP", "XYS",
    # Ions (Single)
    "AG", "AL", "AU", "BR", "CA", "CD", "CL", "CO", "CR", "CS",
    "CU", "F", "FE", "HG", "HO", "I", "IR", "K", "KR", "LI",
    "MG", "MN", "MO", "NA", "NI", "OS", "PB", "PT", "SR", "U",
    "V", "YB", "ZN",
    # Ions (Multi-atom)
    "PO4", "SO4", "WO4",
    # Water
    "TIP", "WAT",
    # Co-factors
    "HEB", "HEC",
    # Nucleic Acids (DNA)
    "DA", "DC", "DG", "DT",
    # Nucleic Acids (RNA)
    "A", "C", "G", "U",
    # Modified Amino Acids (ASH, GLH, MSE already listed above)
    "ACE", "ALY", "CFE", "CSP", "CTN", "CYC", "CYF", "CYM",
    "DDZ", "HY3", "HYP", "M3L", "MLY", "MLZ", "NEP", "NME",
    "PCA", "PNS", "PTR", "QSR", "SEC", "SEP", "TOP", "TYP", "TYS"
}

pool_chars = [c for c in chain_pool if not c.isspace()]

raw_lines = []
for path in input_files:
    if not os.path.isfile(path):
        raise SystemExit(f"File not found: {path}")
    with open(path, 'r', encoding='utf-8', errors='ignore') as handle:
        lines = handle.readlines()
    for line in lines:
        if line.startswith('END'):  # avoid duplicate END
            continue
        raw_lines.append(line.rstrip('\n'))
    if raw_lines and not raw_lines[-1].startswith('TER'):
        raw_lines.append('TER')

if not raw_lines:
    raise SystemExit("Input files produced no content")

chain_order = []
for line in raw_lines:
    if not line.startswith(('ATOM', 'HETATM')):
        continue
    chain_id = line[21].strip() or '_'
    if chain_id not in chain_order:
        chain_order.append(chain_id)

if not chain_order:
    chain_order.append('A')

if len(pool_chars) < len(chain_order):
    raise SystemExit("Not enough unique chain identifiers provided")

chain_mapping = OrderedDict()
for idx, original in enumerate(chain_order):
    if rename_mode == "rename":
        target = pool_chars[idx]
    else:
        target = original if original != '_' else pool_chars[idx]
    chain_mapping[original] = {"new_chain": target, "segid": target}

residue_map = {}
residue_counter = defaultdict(int)
combined_lines = []
atom_serial = 0
current_chain = None

for line in raw_lines:
    record = line[:6].strip()
    if record == 'TER':
        continue  # Will add TER properly between chains
    if record not in {'ATOM', 'HETATM'}:
        continue
    # Skip HETATM records if --remove-hetatm is enabled
    if remove_hetatm and record == 'HETATM':
        continue
    # Handle alternative conformations - only keep first conformation (A or blank)
    altloc = line[16].strip()
    if altloc and altloc != 'A':
        continue  # Skip alternative conformations B, C, etc.
    resn = line[17:20].strip()
    if resn not in STANDARD_RESIDUES:
        continue
    atom_serial += 1
    original_chain = line[21].strip() or '_'
    mapping = chain_mapping.get(original_chain)
    if mapping is None:
        mapping = {"new_chain": pool_chars[len(chain_mapping)], "segid": pool_chars[len(chain_mapping)]}
        chain_mapping[original_chain] = mapping
    
    # Add TER statement when chain changes
    new_chain = mapping["new_chain"]
    if current_chain is not None and current_chain != new_chain:
        combined_lines.append('TER')
    current_chain = new_chain
    
    resseq = line[22:26].strip() or '0'
    inscode = line[26].strip()
    key = (original_chain, resseq, inscode or '_')
    if key not in residue_map:
        if rename_mode == "rename":
            residue_counter[original_chain] += 1
            new_res = residue_counter[original_chain]
        else:
            try:
                new_res = int(resseq)
            except Exception:
                residue_counter[original_chain] += 1
                new_res = residue_counter[original_chain]
            residue_counter[original_chain] = max(residue_counter[original_chain], new_res)
        residue_map[key] = {"chain": mapping["new_chain"], "residue": new_res, "icode": inscode}
    mapped = residue_map[key]

    chars = list(f"{line:<80}")
    chars[6:11] = list(f"{atom_serial:5d}")
    chars[16] = ' '  # Clear alternative location indicator
    chars[21] = mapped['chain'][0]
    chars[22:26] = list(f"{mapped['residue']:4d}")
    chars[26] = mapped['icode'][:1] if mapped['icode'] else ' '
    segid = mapping['segid']
    segid_fmt = segid[:4].rjust(4)
    chars[72:76] = list(segid_fmt)
    new_line = ''.join(chars[:80])
    combined_lines.append(new_line)

# Add final TER before END
if combined_lines and current_chain is not None:
    combined_lines.append('TER')

def ensure_end(lines):
    if not lines:
        return lines
    # Ensure END statement at the end
    if lines[-1] != 'END':
        lines.append('END')
    return lines

combined_path = os.path.join(output_dir, f"{label}_combined.pdb")
with open(combined_path, 'w', encoding='utf-8') as handle:
    ensure_end(combined_lines)
    handle.write('\n'.join(combined_lines))
    handle.write('\n')

mapping_path = os.path.join(output_dir, f"{label}_mapping.json")
residue_dump = {}
for (chain, res, icode), info in residue_map.items():
    residue_dump[f"{chain}:{res}:{icode}"] = info

payload = {
    "label": label,
    "rename_mode": rename_mode,
    "force_chains_together": force_together,
    "inputs": input_files,
    "output_path": os.path.abspath(combined_path),
    "chain_mapping": chain_mapping,
    "residue_map": residue_dump,
}

with open(mapping_path, 'w', encoding='utf-8') as handle:
    json.dump(payload, handle, indent=2)

# No longer creating split chain files - keep multi-chain file
# Body restraints still generated if needed
body_tbl = ''
if force_together and len(chain_mapping) > 1:
    restrain_bin = shutil.which('haddock3-restraints')
    if restrain_bin:
        try:
            target_tbl = os.path.join(output_dir, f"{label}_restrain_bodies.tbl")
            result = subprocess.run(
                [restrain_bin, 'restrain_bodies', combined_path],
                cwd=output_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            candidate = os.path.join(output_dir, 'restrain_bodies.tbl')
            if result.returncode == 0 and result.stdout.strip():
                with open(target_tbl, 'w', encoding='utf-8') as handle:
                    handle.write(result.stdout)
                    if not result.stdout.endswith('\n'):
                        handle.write('\n')
                body_tbl = target_tbl
            elif os.path.exists(candidate):
                os.replace(candidate, target_tbl)
                body_tbl = target_tbl
        except Exception as exc:  # pragma: no cover
            print(f"[WARN] restrain_bodies failed for {label}: {exc}", file=sys.stderr)
    else:
        print("[WARN] haddock3-restraints not available; skipping restrain_bodies", file=sys.stderr)

# Output: combined_path, mapping_path, body_tbl, and always 1 output file (the combined file)
print(os.path.abspath(combined_path))
print(os.path.abspath(mapping_path))
print(body_tbl)
print("1")  # Always output 1 file (combined)
print(os.path.abspath(combined_path))
PY
    )

    (( ${#PREP_INFO[@]} >= 4 )) || die "prepare_structure failed for partner '$label'"

    local combined="${PREP_INFO[0]}"
    local mapping="${PREP_INFO[1]}"
    local body_tbl="${PREP_INFO[2]}"
    local split_count="${PREP_INFO[3]}"

    PARTNER_COMBINED["$label"]="$combined"
    PARTNER_MAPPING["$label"]="$mapping"
    PARTNER_BODY_TBL["$label"]="$body_tbl"

    local -a split_files=()
    local idx
    for (( idx=0; idx<split_count; idx++ )); do
        split_files+=("${PREP_INFO[4+idx]}")
    done
    if (( ${#split_files[@]} == 0 )); then
        split_files=("$combined")
    fi
    PARTNER_SPLITFILES["$label"]=$(printf '%s\n' "${split_files[@]}")
    
    # Run ARCTIC3D analysis if enabled
    if [[ "$USE_ARCTIC3D" == true ]]; then
        run_arctic3d_analysis "$label" "$prep_dir" "$mapping"
    fi
}

# Generate ambiguous and unambiguous restraints from computational/experimental data
generate_auto_restraints() {
    local out_dir="$1"
    local pair_label="$2"
    local mapping_a="$3"
    local mapping_b="$4"
    local combined_a="$5"
    local combined_b="$6"
    local compute_paths="${7:-$COMPUTATIONAL_DIR}"
    local experimental_paths="${8:-$EXPERIMENTAL_DIR}"

    mkdir -p "$out_dir"
    local ambig_path="$out_dir/auto_${pair_label}_ambig.tbl"
    local unambig_path="$out_dir/auto_${pair_label}_unambig.tbl"

    local result
    if ! result=$(
        AUT_MAX_ACTIVE="$MAX_ACTIVE_RESIDUES" \
        AUT_MAX_PAIRS="$MAX_RESTRAINT_PAIRS" \
    python3 - "$compute_paths" "$experimental_paths" "$CONF_THRESHOLD" \
        "$mapping_a" "$mapping_b" "$ambig_path" "$unambig_path" "$combined_a" "$combined_b" <<'PY'
import json
import math
import os
import sys
from collections import defaultdict

comput_dir = sys.argv[1].strip()
exp_dir = sys.argv[2].strip()
conf_threshold = float(sys.argv[3])
map_a_path = sys.argv[4]
map_b_path = sys.argv[5]
ambig_out = sys.argv[6]
unambig_out = sys.argv[7]
combined_a = sys.argv[8]
combined_b = sys.argv[9]

max_active = max(1, int(os.environ.get('AUT_MAX_ACTIVE', '40')))
max_pairs = max(1, int(os.environ.get('AUT_MAX_PAIRS', '4000')))

log_lines = []

def log(msg):
    log_lines.append(msg)

path_sep = os.pathsep
comput_paths = [p for p in comput_dir.split(path_sep) if p.strip()]
exp_paths = [p for p in exp_dir.split(path_sep) if p.strip()]

# Auto-detect: has_experimental determines if we treat exp data as unambiguous
#              has_computational determines if we use randremoval
has_computational = any(os.path.isdir(p) for p in comput_paths)
has_experimental = any(os.path.isdir(p) for p in exp_paths)

# If no data directories at all, skip restraint processing
if not has_computational and not has_experimental:
    print(json.dumps({"ambig": "", "unambig": "", "status": "skip"}))
    sys.exit(0)

with open(map_a_path, 'r', encoding='utf-8') as handle:
    map_a = json.load(handle)
with open(map_b_path, 'r', encoding='utf-8') as handle:
    map_b = json.load(handle)

ident_a = set(map_a.get('inputs', []))
ident_b = set(map_b.get('inputs', []))

# Utility helpers ---------------------------------------------------------

def normalize_score(value):
    if value is None:
        return None
    try:
        score = float(value)
    except Exception:
        return None
    if math.isnan(score) or math.isinf(score):
        return None
    if score > 1.0 and score <= 100.0:
        score /= 100.0
    return max(0.0, min(1.0, score))

PAIR_KEYS = [
    ('chain_a', 'residue_a', 'chain_b', 'residue_b'),
    ('chainA', 'residueA', 'chainB', 'residueB'),
    ('chain1', 'residue1', 'chain2', 'residue2'),
]

SINGLE_KEYS = ['chain', 'chain_id', 'chainId']
RES_KEYS = ['residue', 'resseq', 'resSeq', 'resid']
SCORE_KEYS = ['confidence', 'score', 'probability', 'weight']
DATASET_KEYS = ['dataset', 'source', 'tool']


def parse_record(record):
    """Return tuple (pair, single) where pair is ((chain_a,res_a),(chain_b,res_b))."""
    if not isinstance(record, dict):
        return None, None, None
    pair = None
    for keys in PAIR_KEYS:
        if all(k in record for k in keys):
            pair = (
                (str(record[keys[0]]).strip(), str(record[keys[1]]).strip()),
                (str(record[keys[2]]).strip(), str(record[keys[3]]).strip()),
            )
            break
    chain = None
    residue = None
    for key in SINGLE_KEYS:
        if key in record:
            chain = str(record[key]).strip()
            break
    for key in RES_KEYS:
        if key in record:
            residue = str(record[key]).strip()
            break
    score = None
    for key in SCORE_KEYS:
        if key in record:
            score = normalize_score(record[key])
            break
    dataset = None
    for key in DATASET_KEYS:
        if key in record:
            dataset = str(record[key]).strip()
            break
    return pair, (chain, residue, score, dataset), record


def load_annotations(roots):
    entries = []
    for root in roots:
        if not root or not os.path.isdir(root):
            continue
        for dirpath, _, filenames in os.walk(root):
            for filename in filenames:
                path = os.path.join(dirpath, filename)
                try:
                    if filename.lower().endswith('.json'):
                        with open(path, 'r', encoding='utf-8', errors='ignore') as handle:
                            data = json.load(handle)
                        if isinstance(data, dict):
                            if 'results' in data and isinstance(data['results'], list):
                                data = data['results']
                            elif 'data' in data and isinstance(data['data'], list):
                                data = data['data']
                            elif 'residue_annotations' in data and isinstance(data['residue_annotations'], list):
                                data = data['residue_annotations']
                            else:
                                data = [data]
                    else:
                        with open(path, 'r', encoding='utf-8', errors='ignore') as handle:
                            lines = [line.strip() for line in handle if line.strip() and not line.startswith('#')]
                        sep = ',' if any(',' in line for line in lines) else None
                        if sep:
                            parts = [line.split(sep) for line in lines]
                            data = [dict(enumerate(chunk)) for chunk in parts]
                        else:
                            data = [line.split() for line in lines]
                            data = [dict(enumerate(chunk)) for chunk in data]
                    for entry in data:
                        pair, single, raw = parse_record(entry)
                        if pair is None and (single is None or not single[0] or not single[1]):
                            continue
                        entries.append((path, pair, single, raw))
                except Exception as exc:  # pragma: no cover
                    log(f"Failed to parse {filename}: {exc}")
    return entries


def build_map(meta, chain, residue):
    residue_map = meta.get('residue_map', {}) or {}
    for key, mapped in residue_map.items():
        parts = key.split(':')
        if len(parts) < 2:
            continue
        orig_chain, orig_res = parts[0], parts[1]
        if chain and orig_chain.strip() != chain.strip():
            continue
        if residue and orig_res.strip() != residue.strip():
            continue
        return mapped
    return None

comput_entries = load_annotations(comput_paths)
exp_entries = load_annotations(exp_paths)

ambig_pairs = []
unambig_pairs = []
residues_a = set()
residues_b = set()

with open(map_a_path, 'r', encoding='utf-8') as handle:
    meta_a = json.load(handle)
with open(map_b_path, 'r', encoding='utf-8') as handle:
    meta_b = json.load(handle)

# Process experimental entries as unambiguous restraints (if present)
if has_experimental and exp_entries:
    for path, pair, single, raw in exp_entries:
        if pair:
            mapped_a = build_map(meta_a, pair[0][0], pair[0][1])
            mapped_b = build_map(meta_b, pair[1][0], pair[1][1])
            if mapped_a and mapped_b:
                unambig_pairs.append((mapped_a, mapped_b))
        elif single:
            mapped = build_map(meta_a, single[0], single[1])
            if mapped:
                residues_a.add((mapped['chain'], mapped['residue']))
            mapped = build_map(meta_b, single[0], single[1])
            if mapped:
                residues_b.add((mapped['chain'], mapped['residue']))

comput_residues = defaultdict(set)
for path, pair, single, raw in comput_entries:
    if single:
        chain, residue, score, dataset = single
        if score is not None and score < conf_threshold:
            continue
        mapped = build_map(meta_a, chain, residue)
        if mapped:
            comput_residues['A'].add((mapped['chain'], mapped['residue']))
        mapped = build_map(meta_b, chain, residue)
        if mapped:
            comput_residues['B'].add((mapped['chain'], mapped['residue']))
    if pair:
        mapped_a = build_map(meta_a, pair[0][0], pair[0][1])
        mapped_b = build_map(meta_b, pair[1][0], pair[1][1])
        if mapped_a and mapped_b:
            ambig_pairs.append((mapped_a, mapped_b))

if not ambig_pairs:
    if comput_residues['A'] and comput_residues['B']:
        for a in list(comput_residues['A'])[:max_active]:
            for b in list(comput_residues['B'])[:max_active]:
                ambig_pairs.append(({'chain': a[0], 'residue': a[1]}, {'chain': b[0], 'residue': b[1]}))

ambig_pairs = ambig_pairs[:max_pairs]

if residues_a and residues_b and not unambig_pairs:
    for a in list(residues_a)[:max_active]:
        for b in list(residues_b)[:max_active]:
            unambig_pairs.append(({'chain': a[0], 'residue': a[1]}, {'chain': b[0], 'residue': b[1]}))

ambig_written = ''
if ambig_pairs:
    with open(ambig_out, 'w', encoding='utf-8') as handle:
        handle.write('! Auto-generated ambiguous restraints\n')
        for lhs, rhs in ambig_pairs:
            handle.write(
                f"assign (segid {lhs['chain']} and resid {lhs['residue']}) "
                f"(segid {rhs['chain']} and resid {rhs['residue']}) 2.0 2.0 0.0\n"
            )
    ambig_written = os.path.abspath(ambig_out)
else:
    if os.path.exists(ambig_out):
        os.remove(ambig_out)

unambig_written = ''
if unambig_pairs:
    with open(unambig_out, 'w', encoding='utf-8') as handle:
        handle.write('! Auto-generated unambiguous restraints\n')
        for lhs, rhs in unambig_pairs:
            handle.write(
                f"assign (segid {lhs['chain']} and resid {lhs['residue']}) "
                f"(segid {rhs['chain']} and resid {rhs['residue']}) 0.0 0.0 0.0\n"
            )
    unambig_written = os.path.abspath(unambig_out)
else:
    if os.path.exists(unambig_out):
        os.remove(unambig_out)

print(json.dumps({
    "ambig": ambig_written,
    "unambig": unambig_written,
    "status": "ok"
}))
PY
    ); then
        warn "Automatic restraint generation failed for $pair_label"
        echo ""
        return 1
    fi
    if [[ -z "$result" ]]; then
        echo ""
        return 0
    fi
    readarray -t result_lines <<<"$result"
    if (( ${#result_lines[@]} > 1 )); then
        for (( idx=0; idx<${#result_lines[@]}-1; idx++ )); do
            line="${result_lines[$idx]}"
            [[ -n "$line" ]] && warn "$line"
        done
    fi
    last_index=$(( ${#result_lines[@]} - 1 ))
    printf '%s\n' "${result_lines[$last_index]}"
}

# Create the haddock3.cfg configuration file
create_config() {
    local config_path="$1"
    local run_dir_name="$2"
    local sampling="$3"
    local ncores="$4"
    local ranair="$5"
    local version="$6"
    local ambig_value="$7"
    local unambig_value="$8"
    local reference="$9"
    shift 9
    local -a molecules=("$@")

    python3 - "$config_path" "$run_dir_name" "$sampling" "$ncores" "$ranair" "$version" "$ambig_value" "$unambig_value" "$reference" "${molecules[@]}" <<'PY'
import os
import sys

config_path = sys.argv[1]
run_dir = sys.argv[2]
sampling = int(sys.argv[3])
ncores = int(sys.argv[4])
ranair = sys.argv[5] == 'true'
version = int(sys.argv[6])
ambig_value = sys.argv[7]
unambig_value = sys.argv[8]
reference = sys.argv[9]
molecules = sys.argv[10:]

def parse_files(value):
    return [item for item in value.split(",") if item]

ambig_files = parse_files(ambig_value)
unambig_files = parse_files(unambig_value)

with open(config_path, 'w', encoding='utf-8') as fh:
    fh.write(f"run_dir = \"{run_dir}\"\n")
    fh.write("mode = \"local\"\n")
    fh.write(f"ncores = {ncores}\n\n")

    fh.write("molecules = [\n")
    for mol in molecules:
        fh.write(f"  \"{mol}\",\n")
    fh.write("]\n\n")

    fh.write("[topoaa]\n\n")

    fh.write("[rigidbody]\n")
    fh.write("tolerance = 5\n")
    fh.write(f"sampling = {sampling}\n")
    # Only enable cmrest for ab initio docking (no ambig files)
    if not ambig_files:
        fh.write("cmrest = true\n")
    if ranair:
        fh.write("ranair = true\n")
    if ambig_files:
        if len(ambig_files) == 1:
            fh.write(f"ambig_fname = \"{ambig_files[0]}\"\n")
        else:
            fh.write("ambig_fname = [\n")
            for path in ambig_files:
                fh.write(f"  \"{path}\",\n")
            fh.write("]\n")
    fh.write("\n")

    fh.write("[seletop]\nselect = 50\n\n")  # Changed from 200

    def write_restraints(section):
        fh.write(f"[{section}]\n")
        fh.write("tolerance = 5\n")
        if ambig_files:
            if len(ambig_files) == 1:
                fh.write(f"ambig_fname = \"{ambig_files[0]}\"\n")
            else:
                fh.write("ambig_fname = [\n")
                for path in ambig_files:
                    fh.write(f"  \"{path}\",\n")
                fh.write("]\n")
            # Enable randremoval when using computational predictions
            if has_computational:
                fh.write("randremoval = true\n")
        if unambig_files:
            if len(unambig_files) == 1:
                fh.write(f"unambig_fname = \"{unambig_files[0]}\"\n")
            else:
                fh.write("unambig_fname = [\n")
                for path in unambig_files:
                    fh.write(f"  \"{path}\",\n")
                fh.write("]\n")
        fh.write("\n")

    write_restraints('flexref')
    # write_restraints('mdref')  # Commented out for speedup
    write_restraints('emref')

    if reference:
        fh.write("[caprieval]\n")
        fh.write(f"reference_fname = \"{reference}\"\n\n")

    # Select top 1 conformation and calculate binding affinity
    fh.write("[seletop]\n")
    fh.write("select = 1\n\n")

    fh.write("[prodigyprotein]\n")
PY
}

# ---------------------------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------------------------

require_command python3

ARGS=("$@")
while [[ $# -gt 0 ]]; do
    case "$1" in
        --keep-chains-separate)
            FORCE_CHAINS_TOGETHER=false
            shift
            ;;
        --auto-partners)
            [[ $# -lt 2 ]] && die "--auto-partners expects a directory"
            AUTO_PARTNER_DIR="$2"
            shift 2
            ;;
        --partner)
            [[ $# -lt 2 ]] && die "--partner expects label=spec"
            local_entry="$2"
            [[ "$local_entry" == *=* ]] || die "--partner value must be label=spec"
            local_label="${local_entry%%=*}"
            local_spec="${local_entry#*=}"
            [[ -z "$local_label" ]] && die "Partner label cannot be empty"
            PARTNER_SPECS["$local_label"]="$local_spec"
            PARTNER_ORDER+=("$local_label")
            shift 2
            ;;
        --group-bodies)
            [[ $# -lt 2 ]] && die "--group-bodies expects label=path or json path"
            gb_value="$2"
            if [[ "$gb_value" == *=* ]]; then
                gb_label="${gb_value%%=*}"
                gb_path="${gb_value#*=}"
                GROUP_BODY_SOURCE["$gb_label"]="$gb_path"
            else
                GROUP_BODY_JSON_FILES+=("$gb_value")
            fi
            shift 2
            ;;
        --pair)
            [[ $# -lt 2 ]] && die "--pair expects label1,label2"
            PAIR_OVERRIDES+=("$2")
            shift 2
            ;;
        --pairs-file)
            [[ $# -lt 2 ]] && die "--pairs-file expects path"
            PAIR_FILE="$2"
            shift 2
            ;;
        --ambig)
            [[ $# -lt 2 ]] && die "--ambig expects file list"
            mapfile -t tmp < <(split_csv "$2")
            AMBIG_MANUAL+=("${tmp[@]}")
            shift 2
            ;;
        --unambig)
            [[ $# -lt 2 ]] && die "--unambig expects file list"
            mapfile -t tmp < <(split_csv "$2")
            UNAMBIG_MANUAL+=("${tmp[@]}")
            shift 2
            ;;
        --computational-dir)
            [[ $# -lt 2 ]] && die "--computational-dir expects path"
            COMPUTATIONAL_DIR=$(abs_path "$2")
            shift 2
            ;;
        --experimental-dir)
            [[ $# -lt 2 ]] && die "--experimental-dir expects path"
            EXPERIMENTAL_DIR=$(abs_path "$2")
            shift 2
            ;;
        --out)
            [[ $# -lt 2 ]] && die "--out expects directory"
            OUTPUT_ROOT=$(abs_path "$2")
            shift 2
            ;;
        --project)
            [[ $# -lt 2 ]] && die "--project expects name"
            PROJECT_NAME="$2"
            shift 2
            ;;
        --ncores)
            [[ $# -lt 2 ]] && die "--ncores expects integer"
            NCORES="$2"
            shift 2
            ;;
        --sampling)
            [[ $# -lt 2 ]] && die "--sampling expects integer"
            SAMPLING_OVERRIDE="$2"
            shift 2
            ;;
        --confidence-threshold)
            [[ $# -lt 2 ]] && die "--confidence-threshold expects float"
            CONF_THRESHOLD="$2"
            shift 2
            ;;
        --max-active)
            [[ $# -lt 2 ]] && die "--max-active expects integer"
            MAX_ACTIVE_RESIDUES="$2"
            shift 2
            ;;
        --max-pairs)
            [[ $# -lt 2 ]] && die "--max-pairs expects integer"
            MAX_RESTRAINT_PAIRS="$2"
            shift 2
            ;;
        --reference)
            [[ $# -lt 2 ]] && die "--reference expects path"
            REFERENCE_PDB=$(abs_path "$2")
            shift 2
            ;;
        --abinitio)
            ABINITIO_PRESET=true
            FORCE_RANAIR=true
            shift
            ;;
        --ranair)
            FORCE_RANAIR=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            SKIP_RUN=true
            shift
            ;;
        --use-arctic3d)
            USE_ARCTIC3D=true
            shift
            ;;
        --arctic3d-prob)
            ARCTIC3D_PROB_THRESHOLD="$2"
            shift 2
            ;;
        --arctic3d-blast-db)
            ARCTIC3D_BLAST_DB="$2"
            shift 2
            ;;
        --arctic3d-position-tolerance)
            ARCTIC3D_POSITION_TOLERANCE="$2"
            shift 2
            ;;
        --remove-hetatm)
            REMOVE_HETATM=true
            shift
            ;;
        --run-mode)
            [[ $# -lt 2 ]] && die "--run-mode expects 'sequential' or 'parallel'"
            case "$2" in
                sequential|parallel)
                    RUN_MODE="$2"
                    ;;
                *)
                    die "--run-mode must be 'sequential' or 'parallel' (got '$2')"
                    ;;
            esac
            shift 2
            ;;
        --parallel-jobs)
            [[ $# -lt 2 ]] && die "--parallel-jobs expects an integer"
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Partner Discovery & Validation
# ---------------------------------------------------------------------------

if [[ -n "$AUTO_PARTNER_DIR" ]]; then
    AUTO_PARTNER_DIR=$(abs_path "$AUTO_PARTNER_DIR")
    [[ -d "$AUTO_PARTNER_DIR" ]] || die "Auto partner directory not found: $AUTO_PARTNER_DIR"
    mapfile -t auto_candidates < <(find "$AUTO_PARTNER_DIR" -type f -name "*.pdb" | sort)
    (( ${#auto_candidates[@]} > 0 )) || die "No PDB files found in $AUTO_PARTNER_DIR"
    for file in "${auto_candidates[@]}"; do
        label="$(basename "${file%.*}")"
        label="${label//[^A-Za-z0-9]/_}"
        label="${label^^}"
        [[ "$label" =~ ^[0-9] ]] && label="P${label}"
        suffix=0
        base_label="$label"
        while [[ -n "${PARTNER_SPECS[$label]:-}" ]]; do
            suffix=$((suffix + 1))
            label="${base_label}_${suffix}"
        done
        PARTNER_SPECS["$label"]="$file"
        PARTNER_ORDER+=("$label")
    done
fi

(( ${#PARTNER_ORDER[@]} >= 2 )) || die "Provide at least two partners (via --partner or --auto-partners)."

declare -A _label_seen=()
for label in "${PARTNER_ORDER[@]}"; do
    if [[ -n "${_label_seen[$label]:-}" ]]; then
        die "Duplicate partner label detected: $label"
    fi
    _label_seen["$label"]=1
done

case "$VERSION" in
    1|2|3) ;;
    *) die "Unsupported version '$VERSION'";
        ;;
esac

if (( ${#GROUP_BODY_JSON_FILES[@]} )); then
    for json_file in "${GROUP_BODY_JSON_FILES[@]}"; do
        parse_group_json "$json_file"
    done
fi

for label in "${!GROUP_BODY_SOURCE[@]}"; do
    GROUP_BODY_SOURCE["$label"]=$(abs_path "${GROUP_BODY_SOURCE[$label]}")
done

# Resolve partner specs into absolute file lists
for label in "${PARTNER_ORDER[@]}"; do
    resolve_partner_spec "$label" "${PARTNER_SPECS[$label]}"
    resolve_group_chain_order "$label"
    mapfile -t _files < <(printf '%s\n' "${PARTNER_FILES[$label]}")
    (( ${#_files[@]} >= 1 )) || die "At least one PDB file required for partner '$label'"
    if [[ -n "${GROUP_BODY_SOURCE[$label]:-}" && ! -f "${GROUP_BODY_SOURCE[$label]}" ]]; then
        die "Group source for '$label' not found: ${GROUP_BODY_SOURCE[$label]}"
    fi
done

# Resolve manual restraints to absolute paths
for idx in "${!AMBIG_MANUAL[@]}"; do
    path="${AMBIG_MANUAL[$idx]}"
    [[ -f "$path" ]] || die "Ambiguous restraint file not found: $path"
    AMBIG_MANUAL[$idx]=$(abs_path "$path")
done

for idx in "${!UNAMBIG_MANUAL[@]}"; do
    path="${UNAMBIG_MANUAL[$idx]}"
    [[ -f "$path" ]] || die "Unambiguous restraint file not found: $path"
    UNAMBIG_MANUAL[$idx]=$(abs_path "$path")
done

if [[ -n "$REFERENCE_PDB" ]]; then
    [[ -f "$REFERENCE_PDB" ]] || die "Reference PDB not found: $REFERENCE_PDB"
fi

if [[ -n "$COMPUTATIONAL_DIR" && ! -d "$COMPUTATIONAL_DIR" ]]; then
    die "Computational annotation directory not found: $COMPUTATIONAL_DIR"
fi

if [[ -n "$EXPERIMENTAL_DIR" && ! -d "$EXPERIMENTAL_DIR" ]]; then
    die "Experimental annotation directory not found: $EXPERIMENTAL_DIR"
fi

# ---------------------------------------------------------------------------
# Main Execution: Preparation & Configuration Generation
# ---------------------------------------------------------------------------

# Resolve NCORES (auto-detect or use user-specified value)
NCORES=$(resolve_ncores "$NCORES")

OUTPUT_ROOT=$(abs_path "$OUTPUT_ROOT")
mkdir -p "$OUTPUT_ROOT"

if [[ -z "$PROJECT_NAME" ]]; then
    timestamp=$(date +%Y%m%d_%H%M%S)
    base_project="dock_${timestamp}"
    PROJECT_NAME="$base_project"
    counter=1
    while [[ -d "$OUTPUT_ROOT/$PROJECT_NAME" ]]; do
        PROJECT_NAME="${base_project}_${counter}"
        ((counter+=1))
    done
fi

RUN_DIR="$OUTPUT_ROOT/$PROJECT_NAME"
mkdir -p "$RUN_DIR"

PREP_DIR="$RUN_DIR/prepared"
mkdir -p "$PREP_DIR"

log "Preparing partners in $PREP_DIR"
log "  Partners: ${PARTNER_ORDER[*]}"
log "  Mode: multi-chain files | Version: $VERSION | Force chains together: $FORCE_CHAINS_TOGETHER"
if [[ "$USE_ARCTIC3D" == true ]]; then
    log "  ARCTIC3D: enabled (probability threshold: $ARCTIC3D_PROB_THRESHOLD)"
fi

    idx=0
for label in "${PARTNER_ORDER[@]}"; do
    run_prepare_partner "$label" "$idx" "$PREP_DIR/$label"
    (( idx += 1 ))
done

for label in "${PARTNER_ORDER[@]}"; do
    source_path="${GROUP_BODY_SOURCE[$label]:-}"
    [[ -n "$source_path" ]] || continue
    if [[ ! -f "$source_path" ]]; then
        warn "Group source for '$label' not found: $source_path"
        continue
    fi
    group_dir="$PREP_DIR/$label/group_source"
    mkdir -p "$group_dir"
    group_copy="$group_dir/$(basename "$source_path")"
    cp "$source_path" "$group_copy"
    if command -v haddock3-restraints >/dev/null 2>&1; then
        group_tbl="$group_dir/${label}_restrain_bodies.tbl"
        if haddock3-restraints restrain_bodies "$group_copy" > "$group_tbl" 2> "$group_dir/${label}_restrain_bodies.log"; then
            PARTNER_GROUP_BODY["$label"]="$group_tbl"
        else
            warn "haddock3-restraints failed for '$label' group source ($source_path); see $group_dir/${label}_restrain_bodies.log"
        fi
    else
        warn "haddock3-restraints not available; cannot process --group-bodies source for '$label'"
    fi
done

if [[ -n "$PAIR_FILE" ]]; then
    PAIR_FILE=$(abs_path "$PAIR_FILE")
    [[ -f "$PAIR_FILE" ]] || die "Pairs file not found: $PAIR_FILE"
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line%%#*}"
        line="${line//,/ }"
        line="${line//|/ }"
        line="${line//:/ }"
        read -r lhs rhs _ <<< "$line"
        [[ -z "$lhs" || -z "$rhs" ]] && continue
        [[ -n "${PARTNER_SPECS[$lhs]:-}" ]] || die "Pairs file references unknown partner '$lhs'"
        [[ -n "${PARTNER_SPECS[$rhs]:-}" ]] || die "Pairs file references unknown partner '$rhs'"
        add_pair_unique "$lhs" "$rhs"
    done < "$PAIR_FILE"
fi

for entry in "${PAIR_OVERRIDES[@]}"; do
    cleaned="${entry//,/ }"
    cleaned="${cleaned//|/ }"
    cleaned="${cleaned//:/ }"
    read -r lhs rhs _ <<< "$cleaned"
    [[ -n "$lhs" && -n "$rhs" ]] || die "--pair expects label1,label2 (got '$entry')"
    [[ -n "${PARTNER_SPECS[$lhs]:-}" ]] || die "--pair references unknown partner '$lhs'"
    [[ -n "${PARTNER_SPECS[$rhs]:-}" ]] || die "--pair references unknown partner '$rhs'"
    add_pair_unique "$lhs" "$rhs"
done

if (( ${#PAIR_LIST[@]} == 0 )); then
    for ((i=0; i<${#PARTNER_ORDER[@]}-1; i++)); do
        for ((j=i+1; j<${#PARTNER_ORDER[@]}; j++)); do
            add_pair_unique "${PARTNER_ORDER[i]}" "${PARTNER_ORDER[j]}"
        done
    done
fi

TOTAL_PAIRS=${#PAIR_LIST[@]}
(( TOTAL_PAIRS > 0 )) || die "No partner pairs scheduled for docking."

log "Planned docking pairs:"
for pair in "${PAIR_LIST[@]}"; do
    IFS='|' read -r lhs rhs <<< "$pair"
    log "  - $lhs vs $rhs"
done

if [[ "$SKIP_RUN" != true ]]; then
    require_command haddock3
fi

relative_path() {
    local path="$1"
    echo "${path#$PAIR_DIR/}"
}

join_for_config() {
    local -n arr=$1
    local -a rels=()
    local item
    for item in "${arr[@]}"; do
        rels+=("$(relative_path "$item")")
    done
    local IFS=','
    echo "${rels[*]}"
}

combine_restraint_files() {
    local -n arr=$1
    local dest="$2"
    local header="${3:-}"
    if (( ${#arr[@]} == 0 )); then
        return 0
    fi
    if (( ${#arr[@]} == 1 )); then
        return 0
    fi
    mkdir -p "$(dirname "$dest")"
    python3 - "$dest" "$header" "${arr[@]}" <<'PY'
import sys
from pathlib import Path

dest = Path(sys.argv[1])
header = sys.argv[2]
sources = [Path(arg) for arg in sys.argv[3:]]

with dest.open('w', encoding='utf-8') as handle:
    if header:
        handle.write(f'! {header}\n')
    for src in sources:
        if not src.exists():
            continue
        text = src.read_text(encoding='utf-8')
        if not text:
            continue
        handle.write(text)
        if not text.endswith('\n'):
            handle.write('\n')
if dest.stat().st_size == 0 and header:
    dest.write_text(f'! {header}\n', encoding='utf-8')
PY
    arr=("$dest")
}

# ---------------------------------------------------------------------------
# Pairwise Docking Setup Loop
# ---------------------------------------------------------------------------

success_count=0
failure_count=0
skipped_count=0
guided_count=0
blind_count=0
auto_restraints_count=0
arctic3d_restraints_count=0

pair_idx=0
start_global=$(date +%s)

for pair in "${PAIR_LIST[@]}"; do
    pair_idx=$((pair_idx + 1))
    IFS='|' read -r lhs rhs <<< "$pair"
    pair_label="PAIR_${lhs}_vs_${rhs}"
    
    # ---------------------------------------------------------------------------
    # Check chain count limit BEFORE creating directories
    # ---------------------------------------------------------------------------
    total_chains=$(check_pair_chain_limit "$lhs" "$rhs") || true
    if (( total_chains >= MAX_CHAIN_LIMIT )); then
        skip_reason="Combined chain count ($total_chains) exceeds maximum allowed ($MAX_CHAIN_LIMIT). Chain ID pool (A-Z, 0-9) would be exhausted."
        warn "[$pair_idx/$TOTAL_PAIRS] SKIPPING $lhs vs $rhs: $skip_reason"
        SKIPPED_PAIRS_INFO+=("${lhs}|${rhs}|${total_chains}|${skip_reason}")
        skipped_count=$((skipped_count + 1))
        continue
    fi
    
    PAIR_DIR="$RUN_DIR/$pair_label"
    pair_input_dir="$PAIR_DIR/input"
    pair_rest_dir="$PAIR_DIR/restraints"
    pair_meta_dir="$PAIR_DIR/meta"
    mkdir -p "$PAIR_DIR" "$pair_input_dir" "$pair_rest_dir" "$pair_meta_dir"

    log ""
    log "[$pair_idx/$TOTAL_PAIRS] Preparing $lhs vs $rhs (chains: $total_chains)"

    declare -a molecules=()
    
    # --- Process LHS ---
    partner_label="$lhs"
    mapfile -t split_files < <(printf '%s\n' "${PARTNER_SPLITFILES[$partner_label]}")
    seq_idx=0
    for file in "${split_files[@]}"; do
        base="$(basename "$file")"
        dest_rel="input/${partner_label}_${seq_idx}_${base}"
        dest_path="$PAIR_DIR/$dest_rel"
        mkdir -p "$(dirname "$dest_path")"
        cp "$file" "$dest_path"
        molecules+=("$dest_rel")
        seq_idx=$((seq_idx + 1))
    done
    
    mapping_lhs_src="${PARTNER_MAPPING[$lhs]}"
    mapping_lhs_dest="$pair_meta_dir/${lhs}_mapping.json"
    cp "$mapping_lhs_src" "$mapping_lhs_dest"

    # Calculate LHS chain count to determine RHS start
    lhs_chain_count=$(python3 -c "import json; print(len(json.load(open('$mapping_lhs_dest'))['chain_mapping']))")
    chain_pool="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    rhs_start_char="${chain_pool:$lhs_chain_count:1}"

    # --- Process RHS ---
    partner_label="$rhs"
    mapfile -t split_files < <(printf '%s\n' "${PARTNER_SPLITFILES[$partner_label]}")
    seq_idx=0
    rhs_pdb_files=()
    for file in "${split_files[@]}"; do
        base="$(basename "$file")"
        dest_rel="input/${partner_label}_${seq_idx}_${base}"
        dest_path="$PAIR_DIR/$dest_rel"
        mkdir -p "$(dirname "$dest_path")"
        cp "$file" "$dest_path"
        molecules+=("$dest_rel")
        rhs_pdb_files+=("$dest_path")
        seq_idx=$((seq_idx + 1))
    done

    mapping_rhs_src="${PARTNER_MAPPING[$rhs]}"
    mapping_rhs_dest="$pair_meta_dir/${rhs}_mapping.json"
    cp "$mapping_rhs_src" "$mapping_rhs_dest"

    pair_comput_dir="$PAIR_DIR/computational_data"
    mkdir -p "$pair_comput_dir"
    
    if [[ "$USE_ARCTIC3D" == true ]]; then
        # Copy LHS ARCTIC3D data (chains are unchanged)
        if [[ -d "$PREP_DIR/$lhs/computational_data/$lhs/arctic3d" ]]; then
             mkdir -p "$pair_comput_dir/$lhs"
             cp -r "$PREP_DIR/$lhs/computational_data/$lhs/arctic3d" "$pair_comput_dir/$lhs/"
        fi
        
        # Copy RHS ARCTIC3D data (will be renamed by rechain_partner)
        if [[ -d "$PREP_DIR/$rhs/computational_data/$rhs/arctic3d" ]]; then
             mkdir -p "$pair_comput_dir/$rhs"
             cp -r "$PREP_DIR/$rhs/computational_data/$rhs/arctic3d" "$pair_comput_dir/$rhs/"
        fi
    fi

    declare -a pair_ambig_files=()
    declare -a pair_unambig_files=()

    for path in "${AMBIG_MANUAL[@]}"; do
        base="$(basename "$path")"
        dest="$pair_rest_dir/$base"
        cp "$path" "$dest"
        pair_ambig_files+=("$dest")
    done

    for path in "${UNAMBIG_MANUAL[@]}"; do
        base="$(basename "$path")"
        dest="$pair_rest_dir/$base"
        cp "$path" "$dest"
        pair_unambig_files+=("$dest")
    done

    if [[ "$FORCE_CHAINS_TOGETHER" == true ]]; then
        for partner_label in "$lhs" "$rhs"; do
            body_src="${PARTNER_BODY_TBL[$partner_label]}"
            if [[ -n "$body_src" && -f "$body_src" ]]; then
                body_dest="$pair_rest_dir/${partner_label}_restrain_bodies.tbl"
                cp "$body_src" "$body_dest"
                append_unique pair_unambig_files "$body_dest"
            fi
            group_tbl="${PARTNER_GROUP_BODY[$partner_label]:-}"
            if [[ -n "$group_tbl" && -f "$group_tbl" ]]; then
                group_dest="$pair_rest_dir/${partner_label}_restrain_bodies_source.tbl"
                cp "$group_tbl" "$group_dest"
                append_unique pair_unambig_files "$group_dest"
            fi
        done
    fi

    # Generate ARCTIC3D restraints if enabled
    arctic3d_json=""
    arctic3d_ambig=""
    arctic3d_unambig=""
    if [[ "$USE_ARCTIC3D" == true ]]; then
        log "   Generating ARCTIC3D restraints for $lhs vs $rhs..."
        arctic3d_json=$(generate_arctic3d_restraints_with_footprint "$pair_rest_dir" "$pair_label" \
            "$lhs" "$rhs" "$PREP_DIR" \
            "${PARTNER_COMBINED[$lhs]}" "${PARTNER_COMBINED[$rhs]}")
        
        if [[ -n "$arctic3d_json" ]]; then
            arctic3d_ambig=$(python3 - "$arctic3d_json" <<'PY'
import json, sys
print(json.loads(sys.argv[1]).get("ambig", ""))
PY
)
            arctic3d_unambig=$(python3 - "$arctic3d_json" <<'PY'
import json, sys
print(json.loads(sys.argv[1]).get("unambig", ""))
PY
)
            if [[ -n "$arctic3d_ambig" && -f "$arctic3d_ambig" ]]; then
                append_unique pair_ambig_files "$arctic3d_ambig"
                log "   → ARCTIC3D ambiguous restraints added"
                arctic3d_restraints_count=$((arctic3d_restraints_count + 1))
            fi
            if [[ -n "$arctic3d_unambig" && -f "$arctic3d_unambig" ]]; then
                append_unique pair_unambig_files "$arctic3d_unambig"
                log "   → ARCTIC3D unambiguous restraints added"
                arctic3d_restraints_count=$((arctic3d_restraints_count + 1))
            fi
        fi
    fi

    # --- Rechain RHS (PDBs, JSON, and Arctic Folders) ---
    rhs_arctic_dir=""
    if [[ -d "$pair_comput_dir/$rhs/arctic3d" ]]; then
        rhs_arctic_dir="$pair_comput_dir/$rhs/arctic3d"
    fi
    
    rhs_chain_map=$(rechain_partner "$mapping_rhs_dest" "$rhs_start_char" "$rhs_arctic_dir" "${rhs_pdb_files[@]}")
    
    # --- Remap ARCTIC3D Restraints ---
    if [[ -n "$rhs_chain_map" ]]; then
        if [[ -n "${arctic3d_ambig:-}" && -f "$arctic3d_ambig" ]]; then
            remap_tbl_chain "$arctic3d_ambig" "$rhs_chain_map"
        fi
        if [[ -n "${arctic3d_unambig:-}" && -f "$arctic3d_unambig" ]]; then
            remap_tbl_chain "$arctic3d_unambig" "$rhs_chain_map"
        fi
        
        # Remap Body Restraints for RHS
        if [[ "$FORCE_CHAINS_TOGETHER" == true ]]; then
             body_dest="$pair_rest_dir/${rhs}_restrain_bodies.tbl"
             if [[ -f "$body_dest" ]]; then
                 remap_tbl_chain "$body_dest" "$rhs_chain_map"
             fi
             group_dest="$pair_rest_dir/${rhs}_restrain_bodies_source.tbl"
             if [[ -f "$group_dest" ]]; then
                 remap_tbl_chain "$group_dest" "$rhs_chain_map"
             fi
        fi
    fi

    auto_json=""
    if [[ "$VERSION" -ge 2 ]]; then
        declare -a compute_sources=()
        if [[ -d "$pair_comput_dir" ]]; then
            compute_sources+=("$pair_comput_dir")
        fi
        if [[ -n "$COMPUTATIONAL_DIR" ]]; then
            compute_sources+=("$COMPUTATIONAL_DIR")
        fi
        compute_arg=""
        if (( ${#compute_sources[@]} )); then
            compute_arg=$(IFS=$':'; printf '%s' "${compute_sources[*]}")
        fi
        # Prepare PDBs for auto restraints
        lhs_pdb_auto="${PARTNER_COMBINED[$lhs]}"
        rhs_pdb_auto="$PAIR_DIR/${rhs}_rechained_combined.pdb"
        
        # Combine RHS files (removing END lines except for the last one)
        : > "$rhs_pdb_auto"
        for ((i=0; i<${#rhs_pdb_files[@]}; i++)); do
            if (( i < ${#rhs_pdb_files[@]} - 1 )); then
                grep -v "^END" "${rhs_pdb_files[i]}" >> "$rhs_pdb_auto"
            else
                cat "${rhs_pdb_files[i]}" >> "$rhs_pdb_auto"
            fi
        done

        auto_json=$(generate_auto_restraints "$pair_rest_dir" "$pair_label" \
            "$mapping_lhs_dest" "$mapping_rhs_dest" \
            "$lhs_pdb_auto" "$rhs_pdb_auto" \
            "$compute_arg" "$EXPERIMENTAL_DIR")
        if [[ -n "$auto_json" ]]; then
            auto_ambig=$(python3 - "$auto_json" <<'PY'
import json, sys
print(json.loads(sys.argv[1]).get("ambig", ""))
PY
)
            auto_unambig=$(python3 - "$auto_json" <<'PY'
import json, sys
print(json.loads(sys.argv[1]).get("unambig", ""))
PY
)
            if [[ -n "$auto_ambig" && -f "$auto_ambig" ]]; then
                append_unique pair_ambig_files "$auto_ambig"
                auto_restraints_count=$((auto_restraints_count + 1))
            fi
            if [[ -n "$auto_unambig" && -f "$auto_unambig" ]]; then
                append_unique pair_unambig_files "$auto_unambig"
                auto_restraints_count=$((auto_restraints_count + 1))
            fi
        fi
    fi

    if (( ${#pair_ambig_files[@]} > 1 )); then
        combine_restraint_files pair_ambig_files "$pair_rest_dir/${pair_label}_combined_ambig.tbl" \
            "Combined ambiguous restraints for $pair_label"
    fi
    if (( ${#pair_unambig_files[@]} > 1 )); then
        combine_restraint_files pair_unambig_files "$pair_rest_dir/${pair_label}_combined_unambig.tbl" \
            "Combined unambiguous restraints for $pair_label"
    fi

    if (( ${#pair_ambig_files[@]} > 0 )); then
        guided_count=$((guided_count + 1))
    elif (( ${#pair_unambig_files[@]} > 0 )); then
        guided_count=$((guided_count + 1))
    else
        blind_count=$((blind_count + 1))
    fi

    if [[ -n "$SAMPLING_OVERRIDE" ]]; then
        sampling_value="$SAMPLING_OVERRIDE"
    else
        # Set sampling based on whether ambig files are present
        if [[ ${#pair_ambig_files[@]} -eq 0 ]]; then
            sampling_value=2000  # Ab initio docking (no ambig files)
        else
            sampling_value=1000  # Restraint-guided docking (with ambig files)
        fi
    fi

    ranair_value=false
    if [[ "$FORCE_RANAIR" == true ]]; then
        ranair_value=true
    elif [[ ${#pair_ambig_files[@]} -eq 0 && ${#pair_unambig_files[@]} -eq 0 ]]; then
        ranair_value=true
    fi

    if (( ${#pair_ambig_files[@]} )); then
        AMBIG_VALUE=$(join_for_config pair_ambig_files)
    else
        AMBIG_VALUE=""
    fi
    if (( ${#pair_unambig_files[@]} )); then
        UNAMBIG_VALUE=$(join_for_config pair_unambig_files)
    else
        UNAMBIG_VALUE=""
    fi

    if [[ -n "$REFERENCE_PDB" ]]; then
        ref_dest="$PAIR_DIR/reference.pdb"
        cp "$REFERENCE_PDB" "$ref_dest"
        REFERENCE_VALUE="reference.pdb"
    else
        REFERENCE_VALUE=""
    fi

    CONFIG_PATH="$PAIR_DIR/haddock3.cfg"
    timestamp=$(date +%Y%m%d_%H%M%S)
    base_run="run_${lhs}_vs_${rhs}_${timestamp}"
    RUN_SUBDIR="$base_run"
    counter=1
    while [[ -d "$PAIR_DIR/$RUN_SUBDIR" ]]; do
        RUN_SUBDIR="${base_run}_${counter}"
        ((counter+=1))
    done

    create_config "$CONFIG_PATH" "$RUN_SUBDIR" "$sampling_value" "$NCORES" "$ranair_value" "$VERSION" "$AMBIG_VALUE" "$UNAMBIG_VALUE" "$REFERENCE_VALUE" "${molecules[@]}"

    log "  Config generated at $CONFIG_PATH"
    success_count=$((success_count + 1))
done

# ---------------------------------------------------------------------------
# Execute HADDOCK3 Jobs (if --execute is set)
# ---------------------------------------------------------------------------

if [[ "$EXECUTE_JOBS" == true && "$SKIP_RUN" != true ]]; then
    log ""
    log "======================================"
    log "Executing HADDOCK3 jobs ($RUN_MODE mode)"
    log "======================================"
    
    # Collect all config files
    mapfile -t CONFIG_FILES < <(find "$RUN_DIR" -name 'haddock3.cfg' -type f | sort)
    TOTAL_CONFIGS=${#CONFIG_FILES[@]}
    
    if (( TOTAL_CONFIGS == 0 )); then
        warn "No config files found to execute."
    else
        log "Found $TOTAL_CONFIGS job(s) to execute."
        log ""
        
        executed_count=0
        exec_failure_count=0
        
        if [[ "$RUN_MODE" == "parallel" ]]; then
            # Check if GNU Parallel is available
            if ! command -v parallel >/dev/null 2>&1; then
                warn "GNU Parallel not found. Falling back to sequential mode."
                warn "Install with: apt install parallel (Linux) or brew install parallel (macOS)"
                RUN_MODE="sequential"
            fi
        fi
        
        if [[ "$RUN_MODE" == "sequential" ]]; then
            log "Running jobs sequentially (one at a time)..."
            log ""
            
            job_num=0
            for config in "${CONFIG_FILES[@]}"; do
                job_num=$((job_num + 1))
                pair_dir=$(dirname "$config")
                pair_name=$(basename "$pair_dir")
                
                log "[$job_num/$TOTAL_CONFIGS] Running $pair_name..."
                
                pushd "$pair_dir" >/dev/null
                if haddock3 "$(basename "$config")" > haddock3_run.log 2>&1; then
                    log "  ✓ Completed successfully"
                    executed_count=$((executed_count + 1))
                else
                    error "  ✗ Failed. See $pair_dir/haddock3_run.log"
                    exec_failure_count=$((exec_failure_count + 1))
                fi
                popd >/dev/null
            done
        else
            log "Running jobs in parallel ($PARALLEL_JOBS jobs at a time)..."
            log "Progress will be logged to individual run.log files in each pair directory."
            log ""
            
            # Create a temporary script for parallel execution
            PARALLEL_SCRIPT="$RUN_DIR/.run_haddock3.sh"
            cat > "$PARALLEL_SCRIPT" <<'PSCRIPT'
#!/usr/bin/env bash
config="$1"
pair_dir=$(dirname "$config")
cd "$pair_dir" || exit 1
haddock3 "$(basename "$config")" > haddock3_run.log 2>&1
exit_code=$?
if [[ $exit_code -eq 0 ]]; then
    echo "SUCCESS: $(basename "$pair_dir")"
else
    echo "FAILED: $(basename "$pair_dir") (see $pair_dir/haddock3_run.log)"
fi
exit $exit_code
PSCRIPT
            chmod +x "$PARALLEL_SCRIPT"
            
            # Run with GNU Parallel and capture results
            parallel_output=$(printf '%s\n' "${CONFIG_FILES[@]}" | parallel -j "$PARALLEL_JOBS" --halt never "$PARALLEL_SCRIPT" {} 2>&1) || true
            
            # Parse results
            while IFS= read -r line; do
                if [[ "$line" == SUCCESS:* ]]; then
                    log "  ✓ ${line#SUCCESS: }"
                    executed_count=$((executed_count + 1))
                elif [[ "$line" == FAILED:* ]]; then
                    error "  ✗ ${line#FAILED: }"
                    exec_failure_count=$((exec_failure_count + 1))
                fi
            done <<< "$parallel_output"
            
            rm -f "$PARALLEL_SCRIPT"
        fi
        
        log ""
        log "Execution complete: $executed_count succeeded, $exec_failure_count failed"
    fi
fi

end_global=$(date +%s)
duration=$((end_global - start_global))
hours=$((duration / 3600))
mins=$(((duration % 3600) / 60))

# ---------------------------------------------------------------------------
# Generate Skipped Pairs Summary File
# ---------------------------------------------------------------------------
SKIPPED_SUMMARY_FILE="$RUN_DIR/skipped_pairs_summary.txt"
if (( ${#SKIPPED_PAIRS_INFO[@]} > 0 )); then
    {
        echo "=============================================================================="
        echo "SKIPPED PAIRS SUMMARY"
        echo "Generated: $(date)"
        echo "=============================================================================="
        echo ""
        echo "Total skipped pairs: ${#SKIPPED_PAIRS_INFO[@]}"
        echo ""
        echo "------------------------------------------------------------------------------"
        printf "%-20s %-20s %-10s %s\n" "Partner A" "Partner B" "Chains" "Reason"
        echo "------------------------------------------------------------------------------"
        for entry in "${SKIPPED_PAIRS_INFO[@]}"; do
            IFS='|' read -r skip_lhs skip_rhs skip_chains skip_reason <<< "$entry"
            printf "%-20s %-20s %-10s %s\n" "$skip_lhs" "$skip_rhs" "$skip_chains" "$skip_reason"
        done
        echo "------------------------------------------------------------------------------"
        echo ""
        echo "NOTE: These pairs were skipped because the combined number of chains"
        echo "      exceeds the maximum limit of $MAX_CHAIN_LIMIT chains (A-Z, 0-9)."
        echo "      HADDOCK3 uses single-character chain identifiers, so pairs with"
        echo "      more chains cannot be processed without chain ID collisions."
        echo ""
        echo "POSSIBLE SOLUTIONS:"
        echo "  1. Split large partners into smaller subunits"
        echo "  2. Process chains separately and combine results"
        echo "  3. Use a different docking approach for very large complexes"
    } > "$SKIPPED_SUMMARY_FILE"
    log ""
    log "WARNING: ${#SKIPPED_PAIRS_INFO[@]} pair(s) were skipped due to chain limit."
    log "         See: $SKIPPED_SUMMARY_FILE"
else
    # Create empty summary file to indicate no pairs were skipped
    echo "No pairs were skipped due to chain count limits." > "$SKIPPED_SUMMARY_FILE"
fi

log ""
log "======================================"
log "Workflow summary"
log "======================================"
log "Total runtime: ${hours}h${mins}m"
log "Pairs scheduled: $TOTAL_PAIRS"
log "Successful runs: $success_count"
log "Failed runs:     $failure_count"
log "Skipped (chain limit): $skipped_count"
if [[ "$SKIP_RUN" == true ]]; then
    log "Dry-run only; configurations generated (not executed)."
fi
log "Guided docking:  $guided_count"
log "Blind docking:   $blind_count"
if [[ "$VERSION" -ge 2 ]]; then
    log "Auto restraints: $auto_restraints_count"
fi
if [[ "$USE_ARCTIC3D" == true ]]; then
    log "ARCTIC3D restraints: $arctic3d_restraints_count"
fi
log "Results root:    $RUN_DIR"

if [[ "$EXECUTE_JOBS" == true && "$SKIP_RUN" != true ]]; then
    log ""
    log "Jobs were executed automatically. Check individual pair directories for results."
    log "To re-run a failed job:"
    log "  cd $RUN_DIR/PAIR_<name>"
    log "  haddock3 haddock3.cfg"
else
    log ""
    log "====================================="
    log "How to run the generated jobs:"
    log "====================================="
    log ""
    log "Selected run mode: $RUN_MODE"
    log ""

    if [[ "$RUN_MODE" == "sequential" ]]; then
        log "SEQUENTIAL MODE (recommended for local/weak computers):"
        log "  This mode runs one job at a time, using fewer resources."
        log ""
        log "  Command to run all jobs:"
        log "  for dir in $RUN_DIR/PAIR_*; do"
        log "    echo \"Processing \$dir...\""
        log "    (cd \$dir && haddock3 haddock3.cfg)"
        log "  done"
        log ""
        log "  Or run a single pair:"
        log "  cd $RUN_DIR/PAIR_<name>"
        log "  haddock3 haddock3.cfg"
    else
        log "PARALLEL MODE (recommended for powerful computers/cloud):"
        log "  This mode runs $PARALLEL_JOBS jobs simultaneously using GNU Parallel."
        log "  Requires: GNU Parallel (install via 'apt install parallel' or 'brew install parallel')"
        log ""
        log "  Command to run all jobs in parallel:"
        log "  find $RUN_DIR -name 'haddock3.cfg' | parallel -j $PARALLEL_JOBS 'cd {//} && haddock3 {/} > run.log 2>&1'"
        log ""
        log "  Monitor progress:"
        log "  watch -n 30 'find $RUN_DIR -name \"run.log\" -exec tail -n 1 {} \\;'"
        log ""
        log "  Or run a single pair:"
        log "  cd $RUN_DIR/PAIR_<name>"
        log "  haddock3 haddock3.cfg"
    fi

fi

log ""
log "IMPORTANT: Always cd into the pair directory before running haddock3!"
log ""

exit 0
