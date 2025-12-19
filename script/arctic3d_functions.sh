#!/usr/bin/env bash
# ===========================================================================
# ARCTIC3D Integration Functions for HADDOCK3
# ===========================================================================
# This file contains all ARCTIC3D-related functions used by the main
# HADDOCK3 workflow script. These functions handle:
#   - Checking ARCTIC3D availability (including BioPython dependency)
#   - Extracting UniProt IDs and PDB metadata
#   - Running ARCTIC3D analysis on protein structures
#   - Extracting sequences with residue numbering from PDB files
#   - Performing global sequence alignment for residue mapping
#   - Remapping ARCTIC3D binding site residues to local PDB numbering
#   - Generating interface restraints from ARCTIC3D predictions
#
# Key approach: Uses sequence-based global alignment (BioPython) instead of
# coordinate-based matching to map ARCTIC3D residue numbers to local cleaned
# PDB residue numbers, ensuring accurate binding site restraints for HADDOCK3.
#
# This file should be sourced by the main workflow script.
# ===========================================================================

# ---------------------------------------------------------------------------
# Fallback logging helpers
# ---------------------------------------------------------------------------
# When this file is sourced by HADDOCK3, the parent script typically provides
# `log` and `warn`. For standalone testing, define minimal versions if missing.
# Test command lines:
# script/arctic3d_functions.sh uniprot --pdb protein_test/7KBJ.pdb --chain G
# script/arctic3d_functions.sh extract-seq --pdb protein_test/7KBJ.pdb --chain G
if ! declare -F log >/dev/null 2>&1; then
    log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }
fi
if ! declare -F warn >/dev/null 2>&1; then
    warn() { log "WARN: $*"; }
fi

# Check if ARCTIC3D is installed and available
ensure_arctic3d_available() {
    if [[ "$ARCTIC3D_AVAILABLE" != "unknown" ]]; then
        return
    fi
    if command -v arctic3d >/dev/null 2>&1 && command -v arctic3d-restraints >/dev/null 2>&1; then
        # Check BioPython availability
        if python3 -c "from Bio.Align import PairwiseAligner" 2>/dev/null || python3 -c "from Bio import pairwise2" 2>/dev/null; then
            ARCTIC3D_AVAILABLE="yes"
        else
            warn "BioPython not found; required for ARCTIC3D residue mapping. Install with: pip install biopython"
            ARCTIC3D_AVAILABLE="no"
        fi
    else
        warn "arctic3d or arctic3d-restraints not found in PATH; skipping ARCTIC3D-based restraints."
        ARCTIC3D_AVAILABLE="no"
    fi
}

# Extract uniprot ID from PDB header or sequence
extract_uniprot_id() {
    local pdb_file="$1"
    local chain="$2"
    local label="$3"
    
    [[ -f "$pdb_file" ]] || return 1
    
    # Try to extract from DBREF records first, filtering by chain
    # DBREF format: DBREF  <pdb_id> <chain> <start> <end> UNP <uniprot_id> ...
    # Chain is at column 13 (field position varies due to spacing)
    local uniprot_id
    # First try chain-specific DBREF
    uniprot_id=$(grep "^DBREF" "$pdb_file" | grep " $chain " | grep "UNP" | head -n1 | awk '{for(i=1;i<=NF;i++) if($i=="UNP") print $(i+1)}')
    
    # Debug output to stderr
    echo "DEBUG extract_uniprot_id: file=$pdb_file chain=$chain" >&2
    echo "DEBUG extract_uniprot_id: DBREF lines with chain $chain:" >&2
    grep "^DBREF" "$pdb_file" | grep " $chain " | head -n3 >&2 || echo "  (none found)" >&2
    echo "DEBUG extract_uniprot_id: chain-specific result='$uniprot_id'" >&2
    
    # If no chain-specific match, try any DBREF with UNP (fallback for single-chain files)
    if [[ -z "$uniprot_id" ]]; then
        uniprot_id=$(grep "^DBREF" "$pdb_file" | grep "UNP" | head -n1 | awk '{for(i=1;i<=NF;i++) if($i=="UNP") print $(i+1)}')
        echo "DEBUG extract_uniprot_id: fallback result='$uniprot_id'" >&2
    fi
    
    if [[ -z "$uniprot_id" ]]; then
        # Try COMPND field
        uniprot_id=$(grep "^COMPND" "$pdb_file" | grep -i "uniprot" | head -n1 | sed -n 's/.*[Uu][Nn][Ii][Pp][Rr][Oo][Tt]:\s*\([A-Z0-9]\+\).*/\1/p')
        echo "DEBUG extract_uniprot_id: COMPND result='$uniprot_id'" >&2
    fi
    
    if [[ -z "$uniprot_id" ]]; then
        # Try extracting sequence and running BLAST
        local seq_file="$(dirname "$pdb_file")/${label}_${chain}.fasta"
        python3 - "$pdb_file" "$chain" "$seq_file" "$label" <<'PY'
import sys
from collections import OrderedDict

pdb_file = sys.argv[1]
chain = sys.argv[2]
seq_file = sys.argv[3]
label = sys.argv[4]

sequences = OrderedDict()
with open(pdb_file, 'r', encoding='utf-8', errors='ignore') as handle:
    for line in handle:
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        chain_id = line[21].strip() or '_'
        if chain != 'ALL' and chain_id != chain:
            continue
        resname = line[17:20].strip()
        resseq = line[22:26].strip()
        key = (chain_id, resseq)
        if key not in sequences:
            sequences[key] = resname

# Convert 3-letter to 1-letter
aa_map = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M'
}

seq = ''.join([aa_map.get(resname, 'X') for (_, _), resname in sequences.items()])

if seq:
    with open(seq_file, 'w') as out:
        out.write(f'>{label}_{chain}\n')
        out.write(seq + '\n')
    print(seq_file)
else:
    sys.exit(1)
PY
        
        if [[ -f "$seq_file" ]]; then
            # Run BLAST to find uniprot ID
            if [[ -n "$ARCTIC3D_BLAST_DB" && -f "$ARCTIC3D_BLAST_DB" ]]; then
                uniprot_id=$(blastp -query "$seq_file" -db "$ARCTIC3D_BLAST_DB" -outfmt "6 sacc" -max_target_seqs 1 2>/dev/null | head -n1 | cut -d'|' -f2)
            else
                # Try remote BLAST (slower)
                warn "No local BLAST database specified; attempting remote BLAST for $label chain $chain (this may be slow)"
                uniprot_id=$(blastp -query "$seq_file" -db swissprot -remote -outfmt "6 sacc" -max_target_seqs 1 2>/dev/null | head -n1 | cut -d'|' -f2)
            fi
        fi
    fi
    
    echo "$uniprot_id"
}

# Extract PDB ID and chain info from the original PDB file
extract_pdb_metadata() {
    local pdb_file="$1"
    local chain="$2"
    local label="$3"
    
    [[ -f "$pdb_file" ]] || return 1
    
    # Try to extract PDB ID from header
    local pdb_id
    pdb_id=$(grep "^HEADER" "$pdb_file" | head -n1 | awk '{print $NF}')
    
    if [[ -z "$pdb_id" || ${#pdb_id} != 4 ]]; then
        # Try filename if it looks like a PDB ID
        local filename
        filename=$(basename "$pdb_file" .pdb)
        if [[ "$filename" =~ ^[0-9][a-zA-Z0-9]{3}$ ]]; then
            pdb_id="$filename"
        fi
    fi
    
    # If still no PDB ID, we can't use remote ARCTIC3D easily
    # But we might be able to use local file mode if supported
    echo "$pdb_id"
}

# Run ARCTIC3D on a partner structure using UniProt ID and Original Chain
run_arctic3d_analysis() {
    local label="$1"
    local target_root="$2"
    local mapping_json="$3"

    [[ "$USE_ARCTIC3D" == true ]] || return 0

    ensure_arctic3d_available
    [[ "$ARCTIC3D_AVAILABLE" == "yes" ]] || return 0

    [[ -f "$mapping_json" ]] || return 0
    
    # Extract input files from mapping JSON
    local input_files_str
    input_files_str=$(python3 -c "import sys, json; print('|'.join(json.load(open('$mapping_json'))['inputs']))")
    local -a input_files
    IFS='|' read -r -a input_files <<< "$input_files_str"
    
    if (( ${#input_files[@]} == 0 )); then
        warn "No input files found in mapping for $label"
        return 0
    fi

    local target_dir="$target_root/computational_data/$label/arctic3d"
    rm -rf "$target_dir"
    mkdir -p "$target_dir"
    
    log "   ARCTIC3D: analyzing partner $label..."
    
    # Get chain mapping (Original -> New) from mapping JSON
    local chain_pairs
    chain_pairs=$(python3 - "$mapping_json" <<'PY'
import json, sys
with open(sys.argv[1], 'r') as f:
    data = json.load(f)
mapping = data.get('chain_mapping', {})
for orig, info in mapping.items():
    new_chain = info.get('new_chain', '')
    if new_chain:
        print(f"{orig}|{new_chain}")
PY
    )
    
    if [[ -z "$chain_pairs" ]]; then
        warn "No chain mapping found for $label; skipping ARCTIC3D analysis"
        return 0
    fi
    
    # Store chain-specific outputs (Keyed by NEW chain for downstream compatibility)
    declare -A CHAIN_UNIPROT_IDS
    declare -A CHAIN_ARCTIC3D_DIRS
    declare -A CHAIN_STATUS
    
    # Process each chain pair
    while IFS='|' read -r orig_chain new_chain; do
        log "      Chain $new_chain (Original: $orig_chain): extracting UniProt ID..."
        local uniprot_id=""
        local input_pdb=""
        
        # Try to find UniProt ID from any of the input files
        # Also identify which input file contains this chain
        for infile in "${input_files[@]}"; do
            # Check if chain exists in this file (PDB column 22 is chain ID)
            # Pattern: ATOM (4 chars) + 17 chars = position 21, then chain at column 22
            if grep -q "^ATOM.\{17\}$orig_chain" "$infile"; then
                input_pdb="$infile"
                log "         DEBUG: Found chain $orig_chain in $infile"
                uniprot_id=$(extract_uniprot_id "$infile" "$orig_chain" "$label")
                log "         DEBUG: extract_uniprot_id returned: '$uniprot_id'"
                [[ -n "$uniprot_id" ]] && break
            else
                log "         DEBUG: Chain $orig_chain not found in $infile"
            fi
        done
        
        if [[ -z "$uniprot_id" ]]; then
            warn "Could not determine UniProt ID for $label chain $new_chain (orig: $orig_chain); skipping"
            CHAIN_STATUS["$new_chain"]="missing_uniprot_id"
            continue
        fi
        
        log "      Chain $new_chain: UniProt ID = $uniprot_id (from $input_pdb)"
        CHAIN_UNIPROT_IDS["$new_chain"]="$uniprot_id"
        
        # Extract PDB ID from the input file (ARCTIC3D expects PDB ID, not file path)
        local pdb_id
        pdb_id=$(extract_pdb_metadata "$input_pdb" "$orig_chain" "$label")
        
        if [[ -z "$pdb_id" ]]; then
            warn "Could not determine PDB ID for $label chain $new_chain; skipping ARCTIC3D"
            CHAIN_STATUS["$new_chain"]="missing_pdb_id"
            continue
        fi
        
        log "      Chain $new_chain: PDB ID = $pdb_id"
        
        # Run ARCTIC3D for this chain
        local chain_work_dir="$target_dir/chain_${new_chain}"
        mkdir -p "$chain_work_dir"
        cd "$chain_work_dir"
        
        # Use UniProt ID as main argument
        # Use PDB ID (not path) for --pdb_to_use - ARCTIC3D will fetch it
        # Use Original Chain ID as --chain
        local -a arctic3d_cmd=(arctic3d "$uniprot_id" "--pdb_to_use=$pdb_id" "--chain=$orig_chain")
        
        if [[ -n "$ARCTIC3D_BLAST_DB" ]]; then
            arctic3d_cmd+=("--db" "$ARCTIC3D_BLAST_DB")
        fi
        
        log "      Running: ${arctic3d_cmd[*]}"
        
        if "${arctic3d_cmd[@]}" > "arctic3d_${new_chain}.log" 2>&1; then
            # Find the output directory
            local arctic3d_rundir
            arctic3d_rundir=$(find "$chain_work_dir" -maxdepth 1 -type d -name "arctic3d-*" | head -n1)
            
            if [[ -n "$arctic3d_rundir" && -d "$arctic3d_rundir" ]]; then
                CHAIN_ARCTIC3D_DIRS["$new_chain"]="$arctic3d_rundir"
                CHAIN_STATUS["$new_chain"]="success"
                log "      Chain $new_chain: ARCTIC3D completed → $arctic3d_rundir"
            else
                warn "ARCTIC3D completed but output directory not found for $label chain $new_chain"
                CHAIN_STATUS["$new_chain"]="output_missing"
            fi
        else
            warn "ARCTIC3D failed for $label chain $new_chain; see $chain_work_dir/arctic3d_${new_chain}.log"
            CHAIN_STATUS["$new_chain"]="failed"
        fi
        
        cd - >/dev/null
    done <<< "$chain_pairs"
    
    # Store the chain information for later use in restraint generation
    local manifest_file="$target_dir/chain_manifest.json"
    
    # Prepare arrays as space-separated strings to ensure fixed number of arguments
    local p_keys="${!CHAIN_UNIPROT_IDS[*]}"
    local p_vals="${CHAIN_UNIPROT_IDS[*]}"
    local d_keys="${!CHAIN_ARCTIC3D_DIRS[*]}"
    local d_vals="${CHAIN_ARCTIC3D_DIRS[*]}"
    local s_keys="${!CHAIN_STATUS[*]}"
    local s_vals="${CHAIN_STATUS[*]}"
    
    python3 - "$manifest_file" "$label" "$p_keys" "$p_vals" "$d_keys" "$d_vals" "$s_keys" "$s_vals" <<'PY'
import json, sys

manifest_file = sys.argv[1]
label = sys.argv[2]
uniprot_keys = sys.argv[3].split()
uniprot_values = sys.argv[4].split()
dir_keys = sys.argv[5].split()
dir_values = sys.argv[6].split()
status_keys = sys.argv[7].split()
status_values = sys.argv[8].split()

data = {
    'label': label,
    'chains': {}
}

all_chains = set(uniprot_keys + dir_keys + status_keys)

for chain in all_chains:
    data['chains'][chain] = {
        'uniprot_id': uniprot_values[uniprot_keys.index(chain)] if chain in uniprot_keys else None,
        'arctic3d_dir': dir_values[dir_keys.index(chain)] if chain in dir_keys else None,
        'status': status_values[status_keys.index(chain)] if chain in status_keys else 'unknown'
    }

with open(manifest_file, 'w') as f:
    json.dump(data, f, indent=2)

print(manifest_file)
PY
    
    if [[ -f "$manifest_file" ]]; then
        log "   ARCTIC3D: stored manifest for $label → $manifest_file"
    fi
}

# Extract sequence with residue numbering from PDB file
extract_sequence_with_numbering() {
    local pdb_file="$1"
    local chain="$2"
    local output_json="$3"
    
    python3 - "$pdb_file" "$chain" "$output_json" <<'PY'
import sys, json
from collections import OrderedDict

pdb_file = sys.argv[1]
target_chain = sys.argv[2]
output_file = sys.argv[3]

# 3-letter to 1-letter amino acid mapping
aa_map = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M'
}

residues = OrderedDict()

with open(pdb_file, 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        if not line.startswith('ATOM'):
            continue
        
        chain_id = line[21]
        if chain_id != target_chain:
            continue
        
        res_name = line[17:20].strip()
        res_seq = line[22:26].strip()
        
        try:
            res_num = int(res_seq)
        except ValueError:
            continue
        
        if res_num not in residues:
            residues[res_num] = res_name

# Build sequence and residue list
sequence = ''.join([aa_map.get(res_name, 'X') for res_name in residues.values()])
res_numbers = list(residues.keys())
res_names = list(residues.values())

data = {
    'sequence': sequence,
    'residue_numbers': res_numbers,
    'residue_names': res_names
}

with open(output_file, 'w') as f:
    json.dump(data, f)

print(output_file)
PY
}

# Align sequences and build residue mapping using BioPython
align_and_map_residues() {
    local arctic3d_seq_json="$1"
    local local_seq_json="$2"
    local output_map_json="$3"
    local min_identity="${4:-0.8}"  # Default 80% identity threshold
    
    python3 - "$arctic3d_seq_json" "$local_seq_json" "$output_map_json" "$min_identity" <<'PY'
import sys, json
try:
    from Bio.Align import PairwiseAligner
except Exception:
    PairwiseAligner = None

if PairwiseAligner is None:
    # Fallback for older Biopython versions.
    try:
        from Bio import pairwise2
    except Exception:
        print("ERROR: BioPython not available", file=sys.stderr)
        sys.exit(1)

arctic3d_file = sys.argv[1]
local_file = sys.argv[2]
output_file = sys.argv[3]
min_identity = float(sys.argv[4])

# Load sequence data
with open(arctic3d_file, 'r') as f:
    arctic3d_data = json.load(f)
with open(local_file, 'r') as f:
    local_data = json.load(f)

arctic3d_seq = arctic3d_data['sequence']
local_seq = local_data['sequence']
arctic3d_resnums = arctic3d_data['residue_numbers']
local_resnums = local_data['residue_numbers']

mapping = {}  # arctic3d_resnum -> local_resnum

if PairwiseAligner is not None:
    # Use the non-deprecated API.
    aligner = PairwiseAligner()
    aligner.mode = "global"
    # Match-only scoring (similar to pairwise2.align.globalxx):
    # - match = 1
    # - mismatch = 0
    # - gaps have 0 penalty
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = 0.0
    aligner.extend_gap_score = 0.0

    alignments = aligner.align(arctic3d_seq, local_seq)
    if len(alignments) == 0:
        print("ERROR: No alignment found", file=sys.stderr)
        sys.exit(1)

    aln = alignments[0]
    coords = aln.coordinates

    # coords is a 2xN array of path coordinates.
    a_coords = list(map(int, coords[0]))
    b_coords = list(map(int, coords[1]))

    alignment_length = len(arctic3d_seq)
    matches = 0

    for k in range(len(a_coords) - 1):
        a0, a1 = a_coords[k], a_coords[k + 1]
        b0, b1 = b_coords[k], b_coords[k + 1]
        da = a1 - a0
        db = b1 - b0

        if da > 0 and db > 0:
            if da != db:
                print("ERROR: Unexpected alignment path (diagonal segment length mismatch)", file=sys.stderr)
                sys.exit(1)
            for t in range(da):
                if arctic3d_seq[a0 + t] == local_seq[b0 + t]:
                    matches += 1
                mapping[str(arctic3d_resnums[a0 + t])] = local_resnums[b0 + t]
        # da>0,db==0 => gap in local; da==0,db>0 => gap in arctic; nothing to map

    identity = matches / alignment_length if alignment_length > 0 else 0.0

else:
    # Legacy fallback using pairwise2.
    alignments = pairwise2.align.globalxx(arctic3d_seq, local_seq)

    if not alignments:
        print("ERROR: No alignment found", file=sys.stderr)
        sys.exit(1)

    best_alignment = alignments[0]
    aligned_arctic3d = best_alignment.seqA
    aligned_local = best_alignment.seqB

    matches = sum(1 for a, b in zip(aligned_arctic3d, aligned_local) if a == b and a != '-')
    alignment_length = len([a for a in aligned_arctic3d if a != '-'])
    identity = matches / alignment_length if alignment_length > 0 else 0

    arctic3d_idx = 0
    local_idx = 0
    for a_res, l_res in zip(aligned_arctic3d, aligned_local):
        if a_res != '-' and l_res != '-':
            arctic3d_resnum = arctic3d_resnums[arctic3d_idx]
            local_resnum = local_resnums[local_idx]
            mapping[str(arctic3d_resnum)] = local_resnum
            arctic3d_idx += 1
            local_idx += 1
        elif a_res != '-':
            arctic3d_idx += 1
        elif l_res != '-':
            local_idx += 1

if identity < min_identity:
    print(f"ERROR: Sequence identity {identity:.2%} below threshold {min_identity:.2%}", file=sys.stderr)
    sys.exit(1)

result = {
    'mapping': mapping,
    'identity': identity,
    'alignment_length': alignment_length,
    'matches': matches
}

with open(output_file, 'w') as f:
    json.dump(result, f, indent=2)

print(output_file)
PY
}

# Remap clustered_residues_probs.out file using residue mapping
remap_clustered_residues_probs() {
    local arctic3d_dir="$1"
    local mapping_json="$2"
    local output_dir="$3"
    
    python3 - "$arctic3d_dir" "$mapping_json" "$output_dir" <<'PY'
import sys, json, os, shutil

arctic3d_dir = sys.argv[1]
mapping_file = sys.argv[2]
output_dir = sys.argv[3]

# Load mapping
with open(mapping_file, 'r') as f:
    mapping_data = json.load(f)
mapping = mapping_data['mapping']

# Read clustered_residues_probs.out
input_file = os.path.join(arctic3d_dir, 'clustered_residues_probs.out')
if not os.path.exists(input_file):
    print(f"ERROR: {input_file} not found", file=sys.stderr)
    sys.exit(1)

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Copy mapping file for documentation
shutil.copy(mapping_file, os.path.join(output_dir, 'residue_mapping.json'))

# Remap residue numbers
output_file = os.path.join(output_dir, 'clustered_residues_probs.out')
total_residues = 0
mapped_residues = 0
skipped_residues = []

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
        # Preserve headers and cluster lines
        if line.startswith('Cluster') or line.startswith('rank'):
            f_out.write(line)
            continue
        
        # Skip empty lines
        if not line.strip():
            f_out.write(line)
            continue
        
        # Parse data line: rank\tresid\tresname\tprobability
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            total_residues += 1
            rank = parts[0]
            resid = parts[1]
            resname = parts[2]
            probability = parts[3]
            
            # Remap residue number
            if resid in mapping:
                new_resid = mapping[resid]
                f_out.write(f"{rank}\t{new_resid}\t{resname}\t{probability}\n")
                mapped_residues += 1
            else:
                # Skip residue that cannot be mapped
                skipped_residues.append(f"{resid}({resname})")
        else:
            f_out.write(line)

# Check if at least 50% of residues were mapped
if total_residues > 0:
    mapping_rate = mapped_residues / total_residues
    if mapping_rate < 0.5:
        print(f"ERROR: Only {mapping_rate:.1%} residues mapped (threshold: 50%)", file=sys.stderr)
        sys.exit(1)
    
    if skipped_residues:
        print(f"WARNING: Skipped {len(skipped_residues)} residues: {', '.join(skipped_residues[:10])}", file=sys.stderr)
        if len(skipped_residues) > 10:
            print(f"         ... and {len(skipped_residues) - 10} more", file=sys.stderr)

print(output_dir)
PY
}

# Generate restraints using ARCTIC3D with sequence alignment-based mapping
generate_arctic3d_restraints_with_footprint() {
    # Override log/warn to print to stderr so we don't corrupt the JSON output
    log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }
    warn() { log "WARN: $*" >&2; }

    local out_dir="$1"
    local pair_label="$2"
    local label_a="$3"
    local label_b="$4"
    local prep_dir="$5"
    local combined_a="$6" # Local cleaned PDB A
    local combined_b="$7" # Local cleaned PDB B
    
    mkdir -p "$out_dir"
    local merged_ambig="$out_dir/arctic3d_${pair_label}_ambig.tbl"
    local merged_unambig="$out_dir/arctic3d_${pair_label}_unambig.tbl"
    
    # Check if ARCTIC3D manifests exist
    local manifest_a="$prep_dir/$label_a/computational_data/$label_a/arctic3d/chain_manifest.json"
    local manifest_b="$prep_dir/$label_b/computational_data/$label_b/arctic3d/chain_manifest.json"
    
    if [[ ! -f "$manifest_a" || ! -f "$manifest_b" ]]; then
        log "   ARCTIC3D: manifests not found for $pair_label; skipping ARCTIC3D restraints"
        echo '{"ambig": "", "unambig": "", "status": "skip"}'
        return 0
    fi
    
    # Parse manifests
    # We need to know which chains have successful ARCTIC3D runs
    local chains_a_json
    chains_a_json=$(cat "$manifest_a")
    local chains_b_json
    chains_b_json=$(cat "$manifest_b")
    
    # Helper to get successful chains
    get_success_chains() {
        echo "$1" | python3 -c "import sys, json; data=json.load(sys.stdin); print(' '.join([c for c, i in data.get('chains', {}).items() if i.get('status') == 'success']))"
    }
    
    local -a chains_a=($(get_success_chains "$chains_a_json"))
    local -a chains_b=($(get_success_chains "$chains_b_json"))
    
    if (( ${#chains_a[@]} == 0 || ${#chains_b[@]} == 0 )); then
        log "   ARCTIC3D: no successful chains for $pair_label"
        echo '{"ambig": "", "unambig": "", "status": "skip"}'
        return 0
    fi
    
    local restraints_work_dir="$out_dir/arctic3d_restraints_work"
    mkdir -p "$restraints_work_dir"
    
    local -a valid_restraints=()
    
    # Iterate over all chain pairs
    for chain_a in "${chains_a[@]}"; do
        # Get ARCTIC3D dir for chain A
        local dir_a
        dir_a=$(echo "$chains_a_json" | python3 -c "import sys, json; print(json.load(sys.stdin)['chains']['$chain_a']['arctic3d_dir'])")
        
        # Find main ARCTIC3D PDB file (format: UNIPROT-pdb-chain.pdb, not _cl1/_cl2/_cl3)
        local arctic_pdb_a
        arctic_pdb_a=$(find "$dir_a" -maxdepth 1 -type f -name "*-*-*.pdb" ! -name "*_cl*.pdb" | head -n1)
        
        if [[ ! -f "$arctic_pdb_a" ]]; then
            warn "Main ARCTIC3D PDB not found for $label_a chain $chain_a"
            continue
        fi
        
        log "      Chain $chain_a: ARCTIC3D PDB = $(basename "$arctic_pdb_a")"
        
        # Extract sequences
        local arctic_seq_a="$restraints_work_dir/${label_a}_${chain_a}_arctic_seq.json"
        local local_seq_a="$restraints_work_dir/${label_a}_${chain_a}_local_seq.json"
        
        if [[ ! -f "$arctic_seq_a" ]]; then
            # ARCTIC3D PDB typically has single chain, extract first chain
            local arctic_chain_a
            arctic_chain_a=$(grep "^ATOM" "$arctic_pdb_a" | head -n1 | cut -c22)
            extract_sequence_with_numbering "$arctic_pdb_a" "$arctic_chain_a" "$arctic_seq_a" 2>&1 | grep -v "^$arctic_seq_a$" >&2 || true
        fi
        
        if [[ ! -f "$local_seq_a" ]]; then
            extract_sequence_with_numbering "$combined_a" "$chain_a" "$local_seq_a" 2>&1 | grep -v "^$local_seq_a$" >&2 || true
        fi
        
        # Align sequences and build mapping
        local map_a="$restraints_work_dir/${label_a}_${chain_a}_map.json"
        if [[ ! -f "$map_a" ]]; then
            if align_and_map_residues "$arctic_seq_a" "$local_seq_a" "$map_a" 0.8 2>&1 | grep -v "^$map_a$" >&2; then
                local identity_a
                identity_a=$(python3 -c "import json; print(f\"{json.load(open('$map_a'))['identity']:.1%}\")" 2>/dev/null || echo "N/A")
                log "         Alignment identity: $identity_a"
            else
                warn "Sequence alignment failed for $label_a chain $chain_a (identity < 80% or alignment error)"
                continue
            fi
        fi
        
        # Remap clustered_residues_probs.out
        local remapped_dir_a="$restraints_work_dir/${label_a}_${chain_a}_remapped"
        if [[ ! -d "$remapped_dir_a" ]]; then
            if remap_clustered_residues_probs "$dir_a" "$map_a" "$remapped_dir_a" 2>&1 | grep -v "^$remapped_dir_a$" >&2; then
                log "         Remapped clustered_residues_probs.out → $(basename "$remapped_dir_a")"
            else
                warn "Failed to remap clustered_residues_probs.out for $label_a chain $chain_a"
                continue
            fi
        fi
        
        for chain_b in "${chains_b[@]}"; do
            local dir_b
            dir_b=$(echo "$chains_b_json" | python3 -c "import sys, json; print(json.load(sys.stdin)['chains']['$chain_b']['arctic3d_dir'])")
            
            # Find main ARCTIC3D PDB file for chain B
            local arctic_pdb_b
            arctic_pdb_b=$(find "$dir_b" -maxdepth 1 -type f -name "*-*-*.pdb" ! -name "*_cl*.pdb" | head -n1)
            
            if [[ ! -f "$arctic_pdb_b" ]]; then
                warn "Main ARCTIC3D PDB not found for $label_b chain $chain_b"
                continue
            fi
            
            log "      Chain $chain_b: ARCTIC3D PDB = $(basename "$arctic_pdb_b")"
            
            # Extract sequences
            local arctic_seq_b="$restraints_work_dir/${label_b}_${chain_b}_arctic_seq.json"
            local local_seq_b="$restraints_work_dir/${label_b}_${chain_b}_local_seq.json"
            
            if [[ ! -f "$arctic_seq_b" ]]; then
                local arctic_chain_b
                arctic_chain_b=$(grep "^ATOM" "$arctic_pdb_b" | head -n1 | cut -c22)
                extract_sequence_with_numbering "$arctic_pdb_b" "$arctic_chain_b" "$arctic_seq_b" 2>&1 | grep -v "^$arctic_seq_b$" >&2 || true
            fi
            
            if [[ ! -f "$local_seq_b" ]]; then
                extract_sequence_with_numbering "$combined_b" "$chain_b" "$local_seq_b" 2>&1 | grep -v "^$local_seq_b$" >&2 || true
            fi
            
            # Align sequences and build mapping
            local map_b="$restraints_work_dir/${label_b}_${chain_b}_map.json"
            if [[ ! -f "$map_b" ]]; then
                if align_and_map_residues "$arctic_seq_b" "$local_seq_b" "$map_b" 0.8 2>&1 | grep -v "^$map_b$" >&2; then
                    local identity_b
                    identity_b=$(python3 -c "import json; print(f\"{json.load(open('$map_b'))['identity']:.1%}\")" 2>/dev/null || echo "N/A")
                    log "         Alignment identity: $identity_b"
                else
                    warn "Sequence alignment failed for $label_b chain $chain_b (identity < 80% or alignment error)"
                    continue
                fi
            fi
            
            # Remap clustered_residues_probs.out
            local remapped_dir_b="$restraints_work_dir/${label_b}_${chain_b}_remapped"
            if [[ ! -d "$remapped_dir_b" ]]; then
                if remap_clustered_residues_probs "$dir_b" "$map_b" "$remapped_dir_b" 2>&1 | grep -v "^$remapped_dir_b$" >&2; then
                    log "         Remapped clustered_residues_probs.out → $(basename "$remapped_dir_b")"
                else
                    warn "Failed to remap clustered_residues_probs.out for $label_b chain $chain_b"
                    continue
                fi
            fi
            
            # Generate restraints using REMAPPED directories
            local pair_name="${label_a}_${chain_a}_vs_${label_b}_${chain_b}"
            local pair_dir="$restraints_work_dir/$pair_name"
            
            log "      Generating restraints for $pair_name..."
            
            # Run arctic3d-restraints with REMAPPED directories
            local -a restraint_cmd=(arctic3d-restraints --r1 "$remapped_dir_a" --r2 "$remapped_dir_b" \
                --prob_threshold "$ARCTIC3D_PROB_THRESHOLD" \
                --ch1 "$chain_a" --ch2 "$chain_b" \
                --run_dir "$pair_dir")

            if "${restraint_cmd[@]}" > "$restraints_work_dir/${pair_name}.log" 2>&1; then
                if [[ -d "$pair_dir" ]]; then
                    # Find generated restraints (already have correct residue numbers!)
                    local ambig_tbl
                    ambig_tbl=$(find "$pair_dir" -name "*ambig.tbl" | head -n1)
                    
                    if [[ -f "$ambig_tbl" && -s "$ambig_tbl" ]]; then
                        valid_restraints+=("$ambig_tbl")
                        log "         ✓ Generated restraints with remapped residue numbers"
                    else
                        warn "No valid restraints generated for $pair_name"
                    fi
                fi
            else
                warn "arctic3d-restraints failed for $pair_name; see $restraints_work_dir/${pair_name}.log"
            fi
        done
    done
    
    if (( ${#valid_restraints[@]} == 0 )); then
        log "   ARCTIC3D: no valid restraints generated for $pair_label"
        echo '{"ambig": "", "unambig": "", "status": "no_restraints"}'
        return 0
    fi
    
    # Merge restraints
    {
        echo "! ARCTIC3D-generated ambiguous restraints (sequence alignment-based mapping)"
        echo "! Merged from ${#valid_restraints[@]} chain pair(s)"
        echo "!"
        for tbl_file in "${valid_restraints[@]}"; do
            echo "! From: $(basename "$(dirname "$tbl_file")")"
            grep -v '^!' "$tbl_file" | grep -v '^[[:space:]]*$' || true
        done
    } > "$merged_ambig"
    
    log "   ARCTIC3D: merged ${#valid_restraints[@]} restraint file(s) → $merged_ambig"
    
    # Return JSON result
    echo "{\"ambig\": \"$merged_ambig\", \"unambig\": \"\", \"status\": \"ok\"}"
}

# ---------------------------------------------------------------------------
# Standalone CLI (for testing)
# ---------------------------------------------------------------------------
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    set -euo pipefail

    usage() {
        cat >&2 <<'USAGE'
Standalone test runner for ARCTIC3D integration helpers.

Usage:
  script/arctic3d_functions.sh <command> [options]

Commands:
  check
      Checks whether `arctic3d`, `arctic3d-restraints`, and BioPython are available.

  uniprot --pdb <file.pdb> --chain <A|...|ALL> [--label <name>]
      Attempts to extract a UniProt ID from PDB headers (or via BLAST fallback).

            Arguments:
                --pdb   Path to a PDB file you want to inspect (e.g. protein_test/7KBJ.pdb).
                             This should be the ORIGINAL/unprocessed PDB if you want DBREF/COMPND.
                --chain Chain ID as it appears in the PDB ATOM records (column 22).
                             Example: if the ATOM line has chain 'G', use --chain G.
                             Use ALL to ignore chain filtering when extracting sequence for BLAST.
                --label Optional name used only for naming temporary FASTA files.

  extract-seq --pdb <file.pdb> --chain <A> --out <seq.json>
      Writes a JSON with sequence + residue_numbers + residue_names.

            Arguments:
                --pdb   Path to a PDB file.
                --chain Chain ID in that PDB (column 22 of ATOM lines).
                --out   Output JSON path to write (will be overwritten).

  align-map --arctic-seq <arctic_seq.json> --local-seq <local_seq.json> --out <map.json> [--min-identity <0-1>]
      Runs global alignment and produces a residue-number mapping.

            Arguments:
                --arctic-seq JSON produced by `extract-seq` from the ARCTIC3D PDB
                                        (the PDB downloaded/used by ARCTIC3D; residue numbers are ARCTIC3D's numbering).
                --local-seq  JSON produced by `extract-seq` from your LOCAL/cleaned PDB
                                        (the structure you will dock with; residue numbers are your local numbering).
                --out        Output mapping JSON path to write (will be overwritten).
                --min-identity Minimum allowed sequence identity for accepting the mapping (default: 0.8).

  remap --arctic3d-dir <dir> --mapping <map.json> --outdir <dir>
      Remaps clustered_residues_probs.out using the mapping.

            Arguments:
                --arctic3d-dir Directory containing ARCTIC3D outputs for ONE chain, must include:
                                            clustered_residues_probs.out
                --mapping      Mapping JSON produced by `align-map` (arctic3d_resnum -> local_resnum).
                --outdir       Output directory to write a remapped clustered_residues_probs.out.

  dry-map-remap --arctic3d-dir <dir> --local-pdb <file.pdb> --local-chain <A> [--min-identity <0-1>] [--workdir <dir>]
      Convenience: finds ARCTIC3D main PDB inside <dir>, extracts sequences,
      aligns/maps to the local PDB chain, then remaps clustered_residues_probs.out.

            Arguments:
                --arctic3d-dir Directory containing ARCTIC3D outputs for ONE chain, must include:
                                            clustered_residues_probs.out and an ARCTIC3D PDB (*.pdb).
                --local-pdb    Your local (typically cleaned) PDB you want to map INTO.
                                            Residue numbers from this file are what HADDOCK3 will use.
                --local-chain  Chain ID in the local PDB (column 22 of ATOM lines).
                --workdir      Optional scratch directory. If omitted, a temporary directory is created.

Environment variables (optional):
  ARCTIC3D_BLAST_DB   Path to BLAST database for UniProt ID fallback
  USE_ARCTIC3D        true/false (default: true for standalone)
  ARCTIC3D_AVAILABLE  unknown/yes/no (default: unknown)

Examples (WSL bash):
    # 1) Check whether arctic3d + arctic3d-restraints + BioPython are available.
    ./script/arctic3d_functions.sh check

    # 2) Extract UniProt ID from a local PDB header (DBREF/COMPND), chain = ATOM column 22.
    ./script/arctic3d_functions.sh uniprot --pdb protein_test/7KBJ.pdb --chain G --label 7KBJ

    # 3) Extract sequence+numbering JSON from a PDB chain.
    ./script/arctic3d_functions.sh extract-seq --pdb protein_test/7KBJ.pdb --chain G --out protein_test/_seq_7KBJ_G.json

    # 4) Build an ARCTIC3D->local residue mapping via global alignment.
    #    (Here, we use an existing ARCTIC3D run dir under protein_test/ and map it onto 7KBJ chain G.)
    mkdir -p protein_test/_map_demo
    ARCTIC_PDB=$(find protein_test/arctic3d-Q8BHN3 -maxdepth 1 -type f -name "*-*-*.pdb" ! -name "*_cl*.pdb" | head -n1)
    ARCTIC_CHAIN=$(grep "^ATOM" "$ARCTIC_PDB" | head -n1 | cut -c22)
    ./script/arctic3d_functions.sh extract-seq --pdb "$ARCTIC_PDB" --chain "$ARCTIC_CHAIN" --out protein_test/_map_demo/arctic_seq.json
    ./script/arctic3d_functions.sh extract-seq --pdb protein_test/7KBJ.pdb --chain C --out protein_test/_map_demo/local_seq.json
    ./script/arctic3d_functions.sh align-map --arctic-seq protein_test/_map_demo/arctic_seq.json --local-seq protein_test/_map_demo/local_seq.json --out protein_test/_map_demo/residue_mapping.json --min-identity 0.8

    # 5) Remap clustered_residues_probs.out using the mapping JSON.
    ./script/arctic3d_functions.sh remap --arctic3d-dir protein_test/arctic3d-Q8BHN3 --mapping protein_test/_map_demo/residue_mapping.json --outdir protein_test/_map_demo/remapped

    # 6) One-shot pipeline: find ARCTIC3D PDB in a run dir -> extract seqs -> align/map -> remap.
    ./script/arctic3d_functions.sh dry-map-remap --arctic3d-dir protein_test/arctic3d-Q8BHN3 --local-pdb protein_test/7KBJ.pdb --local-chain G --workdir protein_test/_dry_q8bhn3

Examples (Windows PowerShell -> WSL):
    wsl -d Ubuntu -- bash -lc "cd /home/phucdao/haddock3-znd && ./script/arctic3d_functions.sh uniprot --pdb protein_test/7KBJ.pdb --chain G --label 7KBJ"

USAGE
    }

    die() { printf 'ERROR: %s\n' "$*" >&2; exit 2; }

    need_arg() {
        local name="$1"; local value="$2"
        [[ -n "${value:-}" ]] || die "Missing required argument: $name"
    }

    # Defaults for standalone runs (do not affect sourced usage)
    : "${USE_ARCTIC3D:=true}"
    : "${ARCTIC3D_AVAILABLE:=unknown}"
    : "${ARCTIC3D_PROB_THRESHOLD:=0.5}"

    cmd="${1:-}"
    shift || true

    case "$cmd" in
        ""|"-h"|"--help"|"help")
            usage
            exit 0
            ;;

        check)
            ensure_arctic3d_available
            echo "ARCTIC3D_AVAILABLE=$ARCTIC3D_AVAILABLE"
            ;;

        uniprot)
            # --pdb: path to PDB to inspect
            # --chain: chain ID in that PDB (ATOM column 22)
            pdb=""; chain=""; label="test"
            while (( $# )); do
                case "$1" in
                    --pdb) pdb="${2:-}"; shift 2 ;;
                    --chain) chain="${2:-}"; shift 2 ;;
                    --label) label="${2:-}"; shift 2 ;;
                    -h|--help) usage; exit 0 ;;
                    *) die "Unknown option for uniprot: $1" ;;
                esac
            done
            need_arg --pdb "$pdb"
            need_arg --chain "$chain"
            [[ -f "$pdb" ]] || die "PDB file not found: $pdb"
            if ! extract_uniprot_id "$pdb" "$chain" "$label"; then
                die "UniProt ID extraction failed for pdb=$pdb chain=$chain"
            fi
            ;;

        extract-seq)
            # --pdb: path to PDB to extract sequence from
            # --chain: chain ID in that PDB (ATOM column 22)
            # --out: output JSON path
            pdb=""; chain=""; out=""
            while (( $# )); do
                case "$1" in
                    --pdb) pdb="${2:-}"; shift 2 ;;
                    --chain) chain="${2:-}"; shift 2 ;;
                    --out) out="${2:-}"; shift 2 ;;
                    -h|--help) usage; exit 0 ;;
                    *) die "Unknown option for extract-seq: $1" ;;
                esac
            done
            need_arg --pdb "$pdb"
            need_arg --chain "$chain"
            need_arg --out "$out"
            [[ -f "$pdb" ]] || die "PDB file not found: $pdb"
            if ! extract_sequence_with_numbering "$pdb" "$chain" "$out" >/dev/null; then
                die "Sequence extraction failed for pdb=$pdb chain=$chain"
            fi
            echo "$out"
            ;;

        align-map)
            # --arctic-seq: JSON from extract-seq run on ARCTIC3D PDB
            # --local-seq:  JSON from extract-seq run on local/cleaned PDB
            # --out: output mapping JSON
            arctic_seq=""; local_seq=""; out=""; min_identity="0.8"
            while (( $# )); do
                case "$1" in
                    --arctic-seq) arctic_seq="${2:-}"; shift 2 ;;
                    --local-seq) local_seq="${2:-}"; shift 2 ;;
                    --out) out="${2:-}"; shift 2 ;;
                    --min-identity) min_identity="${2:-}"; shift 2 ;;
                    -h|--help) usage; exit 0 ;;
                    *) die "Unknown option for align-map: $1" ;;
                esac
            done
            need_arg --arctic-seq "$arctic_seq"
            need_arg --local-seq "$local_seq"
            need_arg --out "$out"
            [[ -f "$arctic_seq" ]] || die "File not found: $arctic_seq"
            [[ -f "$local_seq" ]] || die "File not found: $local_seq"
            if ! align_and_map_residues "$arctic_seq" "$local_seq" "$out" "$min_identity" >/dev/null; then
                die "Alignment/mapping failed (min_identity=$min_identity)"
            fi
            echo "$out"
            ;;

        remap)
            # --arctic3d-dir: directory containing clustered_residues_probs.out
            # --mapping: mapping JSON from align-map
            # --outdir: output directory
            arctic_dir=""; mapping=""; outdir=""
            while (( $# )); do
                case "$1" in
                    --arctic3d-dir) arctic_dir="${2:-}"; shift 2 ;;
                    --mapping) mapping="${2:-}"; shift 2 ;;
                    --outdir) outdir="${2:-}"; shift 2 ;;
                    -h|--help) usage; exit 0 ;;
                    *) die "Unknown option for remap: $1" ;;
                esac
            done
            need_arg --arctic3d-dir "$arctic_dir"
            need_arg --mapping "$mapping"
            need_arg --outdir "$outdir"
            [[ -d "$arctic_dir" ]] || die "Not a directory: $arctic_dir"
            [[ -f "$mapping" ]] || die "Mapping file not found: $mapping"
            if ! remap_clustered_residues_probs "$arctic_dir" "$mapping" "$outdir" >/dev/null; then
                die "Remap failed for arctic3d_dir=$arctic_dir"
            fi
            echo "$outdir"
            ;;

        dry-map-remap)
            # Convenience: extract seqs, align/map, remap clustered_residues_probs.out
            # --arctic3d-dir: ARCTIC3D output directory containing PDB + clustered_residues_probs.out
            # --local-pdb: local/cleaned PDB you will dock
            # --local-chain: chain ID in local PDB
            arctic_dir=""; local_pdb=""; local_chain=""; min_identity="0.8"; workdir=""
            while (( $# )); do
                case "$1" in
                    --arctic3d-dir) arctic_dir="${2:-}"; shift 2 ;;
                    --local-pdb) local_pdb="${2:-}"; shift 2 ;;
                    --local-chain) local_chain="${2:-}"; shift 2 ;;
                    --min-identity) min_identity="${2:-}"; shift 2 ;;
                    --workdir) workdir="${2:-}"; shift 2 ;;
                    -h|--help) usage; exit 0 ;;
                    *) die "Unknown option for dry-map-remap: $1" ;;
                esac
            done
            need_arg --arctic3d-dir "$arctic_dir"
            need_arg --local-pdb "$local_pdb"
            need_arg --local-chain "$local_chain"

            [[ -d "$arctic_dir" ]] || die "Not a directory: $arctic_dir"
            [[ -f "$local_pdb" ]] || die "Local PDB not found: $local_pdb"

            # Locate an ARCTIC3D main PDB file in the run directory.
            arctic_pdb=$(find "$arctic_dir" -maxdepth 1 -type f -name "*-*-*.pdb" ! -name "*_cl*.pdb" | head -n1 || true)
            if [[ -z "$arctic_pdb" ]]; then
                arctic_pdb=$(find "$arctic_dir" -maxdepth 1 -type f -name "*.pdb" | head -n1 || true)
            fi
            [[ -n "$arctic_pdb" && -f "$arctic_pdb" ]] || die "No PDB found in: $arctic_dir"

            # Guess ARCTIC chain ID from first ATOM line.
            arctic_chain=$(grep "^ATOM" "$arctic_pdb" | head -n1 | cut -c22 | tr -d '[:space:]' || true)
            [[ -n "$arctic_chain" ]] || die "Could not determine chain ID from: $arctic_pdb"

            if [[ -z "$workdir" ]]; then
                workdir=$(mktemp -d)
            else
                mkdir -p "$workdir"
            fi

            log "Using ARCTIC3D PDB: $arctic_pdb (chain $arctic_chain)"
            log "Using local PDB:   $local_pdb (chain $local_chain)"
            log "Workdir:           $workdir"

            arctic_seq="$workdir/arctic_seq.json"
            local_seq="$workdir/local_seq.json"
            map_json="$workdir/residue_mapping.json"
            outdir="$workdir/remapped"

            extract_sequence_with_numbering "$arctic_pdb" "$arctic_chain" "$arctic_seq" >/dev/null
            extract_sequence_with_numbering "$local_pdb" "$local_chain" "$local_seq" >/dev/null
            if ! align_and_map_residues "$arctic_seq" "$local_seq" "$map_json" "$min_identity" >/dev/null; then
                die "Alignment/mapping failed (min_identity=$min_identity)"
            fi
            if ! remap_clustered_residues_probs "$arctic_dir" "$map_json" "$outdir" >/dev/null; then
                die "Remap failed for arctic3d_dir=$arctic_dir"
            fi

            printf 'arctic_pdb=%s\n' "$arctic_pdb"
            printf 'arctic_chain=%s\n' "$arctic_chain"
            printf 'workdir=%s\n' "$workdir"
            printf 'arctic_seq=%s\n' "$arctic_seq"
            printf 'local_seq=%s\n' "$local_seq"
            printf 'map_json=%s\n' "$map_json"
            printf 'remapped_dir=%s\n' "$outdir"
            ;;

        *)
            usage
            die "Unknown command: ${cmd:-<none>}"
            ;;
    esac
fi
