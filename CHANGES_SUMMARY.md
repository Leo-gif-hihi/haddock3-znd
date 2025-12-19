# Summary of Changes: Replacing fpocket with ARCTIC3D

## Overview

The `a.sh` script has been updated to replace fpocket with ARCTIC3D for generating protein-protein and protein-ligand docking restraints. This change addresses the limitation that arctic3d-restraints only accepts two inputs at a time, implementing a "divide and conquer" strategy for multi-chain complexes.

## Key Changes Made

### 1. Global Configuration Variables (Lines 71-76)

**Before:**
```bash
# Fpocket Settings
USE_FPOCKET=false
FPOCKET_TOP_N=3
FPOCKET_MIN_DRUG=0.5
FPOCKET_MIN_VOLUME=50
FPOCKET_AVAILABLE="unknown"
```

**After:**
```bash
# ARCTIC3D Settings
USE_ARCTIC3D=false
ARCTIC3D_PROB_THRESHOLD=0.4
ARCTIC3D_AVAILABLE="unknown"
ARCTIC3D_BLAST_DB=""
```

### 2. Command-Line Options (Lines ~112 and ~1547-1560)

**Added Options:**
- `--use-arctic3d`: Enable ARCTIC3D restraint generation
- `--arctic3d-prob <float>`: Probability threshold (default: 0.4)
- `--arctic3d-blast-db <path>`: Path to local BLAST database

**Removed Options:**
- `--use-fpocket`

### 3. Uniprot ID Extraction Function (Lines ~613-705)

**New Function: `extract_uniprot_id()`**

This function extracts uniprot IDs from PDB files through multiple methods:

1. **DBREF records**: Parses standard PDB database reference records
2. **COMPND fields**: Searches for uniprot references in compound descriptions
3. **Sequence BLAST**: Falls back to BLAST search (local or remote) if direct extraction fails

Key features:
- Supports both local and remote BLAST
- Extracts sequences from PDB ATOM records
- Converts 3-letter amino acid codes to 1-letter codes
- Handles multiple chains

### 4. ARCTIC3D Prediction Function (Lines ~705-850)

**New Function: `generate_arctic3d_predictions()`**

Replaces `generate_fpocket_predictions()` with comprehensive ARCTIC3D workflow:

**Workflow:**
1. Parses chain mapping from JSON to identify all chains
2. Extracts uniprot ID for each chain independently
3. Creates exclusion list (other chains in same partner)
4. Runs `arctic3d` for each chain:
   ```bash
   arctic3d <uniprot_id> --pdb_to_use=<pdb> --chain=<chain> --out_uniprot=<exclude_list>
   ```
5. Stores results in structured directories
6. Creates a manifest JSON file with chain-to-directory mappings

**Key improvements over fpocket:**
- Chain-aware analysis
- Automatic exclusion of intra-complex interfaces
- Structured output with full ARCTIC3D run data
- Manifest file for downstream restraint generation

### 5. ARCTIC3D Restraint Generation Function (Lines ~850-1018)

**New Function: `generate_arctic3d_restraints()`**

Implements the "divide and conquer" strategy for multi-chain restraints:

**Workflow:**
1. Loads chain manifests for both partners
2. Iterates over all chain pairs (A₁-B₁, A₁-B₂, A₂-B₁, A₂-B₂, etc.)
3. For each chain pair, runs:
   ```bash
   arctic3d-restraints --r1 <dir_A> --r2 <dir_B> \
     --prob_threshold <threshold> --ch1 <chain_A> --ch2 <chain_B>
   ```
4. Collects all generated .tbl files
5. Merges restraints into two master files:
   - `arctic3d_<pair>_ambig.tbl`: Combined ambiguous restraints
   - `arctic3d_<pair>_unambig.tbl`: Combined unambiguous restraints

**Merging Logic:**
- Identifies ambiguous vs unambiguous files by name pattern
- Preserves comments indicating source of each restraint
- Filters out comment lines during merge
- Maintains proper HADDOCK restraint format

### 6. Main Execution Flow Updates (Lines ~1680-1700)

**Partner Preparation:**
```bash
# Before
if [[ "$USE_FPOCKET" == true ]]; then
    fpocket_target="$PREP_DIR/$label/computational_data"
    mkdir -p "$fpocket_target"
    generate_fpocket_predictions "$label" "${PARTNER_COMBINED[$label]}" "$fpocket_target"
fi

# After
if [[ "$USE_ARCTIC3D" == true ]]; then
    arctic3d_target="$PREP_DIR/$label/computational_data"
    mkdir -p "$arctic3d_target"
    generate_arctic3d_predictions "$label" "${PARTNER_COMBINED[$label]}" \
        "$arctic3d_target" "${PARTNER_MAPPING[$label]}"
fi
```

### 7. Pairwise Docking Loop Updates (Lines ~1880-1935)

**Data Copying:**
```bash
# Before
if [[ "$USE_FPOCKET" == true ]]; then
    cp -r "$PREP_DIR/$lhs/computational_data/$lhs/fpocket" "$pair_comput_dir/$lhs/"
fi

# After
if [[ "$USE_ARCTIC3D" == true ]]; then
    cp -r "$PREP_DIR/$lhs/computational_data/$lhs/arctic3d" "$pair_comput_dir/$lhs/"
fi
```

**Restraint Generation (NEW):**
```bash
# Generate ARCTIC3D restraints if enabled
arctic3d_json=""
if [[ "$USE_ARCTIC3D" == true ]]; then
    log "   Generating ARCTIC3D restraints for $lhs vs $rhs..."
    arctic3d_json=$(generate_arctic3d_restraints "$pair_rest_dir" "$pair_label" \
        "$lhs" "$rhs" "$PREP_DIR" \
        "$mapping_lhs_dest" "$mapping_rhs_dest")
    
    # Parse results and add to restraint files
    if [[ -n "$arctic3d_json" ]]; then
        arctic3d_ambig=$(extract ambig from JSON)
        arctic3d_unambig=$(extract unambig from JSON)
        append_unique pair_ambig_files "$arctic3d_ambig"
        append_unique pair_unambig_files "$arctic3d_unambig"
        arctic3d_restraints_count=$((arctic3d_restraints_count + 1))
    fi
fi
```

### 8. Summary Statistics (Lines ~2100-2110)

**Added Counter:**
```bash
arctic3d_restraints_count=0
```

**Summary Output:**
```bash
if [[ "$USE_ARCTIC3D" == true ]]; then
    log "ARCTIC3D restraints: $arctic3d_restraints_count"
fi
```

## Critical Design Decisions

### 1. Chain-Level Processing

**Problem:** ARCTIC3D only accepts two partners at a time, but complexes have multiple chains.

**Solution:** Process each chain independently:
- Extract uniprot ID per chain
- Run ARCTIC3D per chain with exclusions
- Generate pairwise restraints for all chain combinations
- Merge all restraints into master files

### 2. Exclusion Strategy

**Problem:** When docking complex A+B to ligand C, we don't want A-B interface restraints.

**Solution:** Use `--out_uniprot` to exclude other chains:
```bash
# For chain A, exclude chain B's uniprot
arctic3d Q8BHN3 --out_uniprot=O08795 ...
```

### 3. Restraint Mapping Preservation

**Problem:** Chain renaming could break restraints.

**Solution:** 
- Store original→new chain mapping in JSON
- ARCTIC3D restraints use the **new** chain IDs
- Mapping is applied during `arctic3d-restraints` call via `--ch1` and `--ch2`

### 4. Restraint Merging

**Problem:** Multiple `.tbl` files from different chain pairs need to be combined.

**Solution:**
- Classify files as ambiguous/unambiguous by name pattern
- Merge by concatenating non-comment lines
- Add source comments for traceability

## Benefits Over fpocket

1. **Interface-Specific**: ARCTIC3D finds actual protein-protein interfaces from database knowledge
2. **Multi-Chain Aware**: Handles complex topologies automatically
3. **Exclusion Support**: Prevents unwanted intra-complex restraints
4. **Comprehensive Output**: Full ARCTIC3D analysis available for inspection
5. **Restraint Quality**: Based on experimental interface data, not just pockets

## Backward Compatibility

- All existing options remain functional
- Old workflow (without `--use-arctic3d`) unchanged
- Can still use manual restraints via `--ambig` and `--unambig`
- Version 1/2/3 workflows still supported

## Testing Recommendations

1. **Simple Two-Partner**: Test basic A+B docking
2. **Multi-Chain Complex**: Test A+B vs C docking
3. **Without Uniprot IDs**: Test BLAST fallback
4. **Local BLAST**: Test with `--arctic3d-blast-db`
5. **Multiple Pairs**: Test `--pair A,B --pair A,C` with ARCTIC3D
6. **Edge Cases**: Empty chains, missing data, failed ARCTIC3D runs

## Documentation Created

1. **ARCTIC3D_USAGE.md**: Comprehensive user guide with examples
2. **This file**: Technical summary of implementation changes

## Future Enhancements

Possible improvements:

1. **Parallel ARCTIC3D**: Run chains in parallel for speed
2. **Cache Uniprot IDs**: Store extracted IDs to avoid re-BLASTing
3. **Restraint Statistics**: Report number of restraints per source
4. **Interface Visualization**: Generate plots of predicted interfaces
5. **Confidence Filtering**: Allow per-restraint confidence thresholds
