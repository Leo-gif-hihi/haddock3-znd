# ARCTIC3D Functions Refactoring

## Overview
This document describes the refactoring of ARCTIC3D-related code in the HADDOCK3 workflow script for improved code organization and maintainability.

## Changes Made

### 1. New File Created: `arctic3d_functions.sh`
All ARCTIC3D-related functions have been extracted into a separate file located at:
```
script/arctic3d_functions.sh
```

### 2. Functions Moved
The following functions were moved from `a.sh` to `arctic3d_functions.sh`:

- **`ensure_arctic3d_available()`** - Checks if ARCTIC3D commands are available in PATH
- **`extract_uniprot_id()`** - Extracts UniProt IDs from PDB files using various methods (DBREF, COMPND, BLAST)
- **`extract_pdb_metadata()`** - Extracts PDB ID from file headers or filenames
- **`run_arctic3d_analysis()`** - Runs ARCTIC3D analysis on partner structures
- **`build_residue_footprint()`** - Creates spatial footprint of residues for structure matching
- **`match_footprint_to_local_pdb()`** - Matches ARCTIC3D output to local PDB structures using footprints
- **`generate_arctic3d_restraints_with_footprint()`** - Generates interface restraints from ARCTIC3D predictions

### 3. Main Script Changes
In `a.sh`, the following changes were made:

#### Added Sourcing Logic (around line 498-510)
```bash
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
```

#### Removed Function Definitions
All ARCTIC3D function definitions (approximately 600 lines) were removed from the main script.

### 4. Function Calls Remain Unchanged
All calls to ARCTIC3D functions in the main script remain the same:
- `run_arctic3d_analysis()` is called during partner preparation (line ~788)
- `generate_arctic3d_restraints_with_footprint()` is called during restraint generation (line ~1726)

## Benefits

### 1. **Improved Readability**
- Main script is now ~600 lines shorter
- ARCTIC3D logic is isolated in its own file
- Easier to understand the main workflow without ARCTIC3D details

### 2. **Better Maintainability**
- ARCTIC3D code can be updated independently
- Easier to test ARCTIC3D functions in isolation
- Clear separation of concerns

### 3. **Modular Design**
- ARCTIC3D functionality is now a pluggable module
- Other features could be similarly modularized
- Easier to disable/enable ARCTIC3D support

### 4. **Error Handling**
- Script gracefully handles missing arctic3d_functions.sh
- Warns user if ARCTIC3D is requested but functions file is not found
- Falls back to non-ARCTIC3D mode automatically

## Usage

### No Changes Required for End Users
The script works exactly as before. All command-line options remain the same:

```bash
./a.sh --partner A=proteinA.pdb --partner B=proteinB.pdb \
       --use-arctic3d --arctic3d-prob 0.4
```

### For Developers
When modifying ARCTIC3D functionality:
1. Edit `script/arctic3d_functions.sh` instead of the main script
2. Keep the function signatures unchanged to maintain compatibility
3. The functions have access to all global variables from the main script
4. Use `log()` and `warn()` helper functions for consistent output

## File Structure

```
script/
├── a.sh                      # Main workflow script (sources arctic3d_functions.sh)
└── arctic3d_functions.sh     # ARCTIC3D integration functions (new file)
```

## Testing

After this refactoring, test the following scenarios:

1. **With ARCTIC3D enabled**: Verify ARCTIC3D functions work correctly
   ```bash
   ./a.sh --partner A=A.pdb --partner B=B.pdb --use-arctic3d
   ```

2. **Without ARCTIC3D**: Verify script works without --use-arctic3d flag
   ```bash
   ./a.sh --partner A=A.pdb --partner B=B.pdb
   ```

3. **Missing functions file**: Temporarily rename arctic3d_functions.sh and verify graceful fallback
   ```bash
   mv arctic3d_functions.sh arctic3d_functions.sh.bak
   ./a.sh --partner A=A.pdb --partner B=B.pdb --use-arctic3d
   # Should warn and continue without ARCTIC3D
   ```

## Global Variables Used by ARCTIC3D Functions

The following global variables from the main script are accessed by ARCTIC3D functions:
- `USE_ARCTIC3D` - Flag to enable/disable ARCTIC3D
- `ARCTIC3D_AVAILABLE` - Cache of ARCTIC3D availability check
- `ARCTIC3D_BLAST_DB` - Path to BLAST database
- `ARCTIC3D_PROB_THRESHOLD` - Probability threshold for restraints
- `ARCTIC3D_POSITION_TOLERANCE` - Tolerance for position matching

These variables are set in the main script and used by the sourced functions.

## Future Improvements

Potential future enhancements:
1. Move ARCTIC3D global variable initialization to arctic3d_functions.sh
2. Create similar modular files for other features (e.g., restraint generation)
3. Add unit tests for ARCTIC3D functions
4. Document ARCTIC3D workflow in separate documentation
