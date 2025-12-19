#!/bin/bash

extract_uniprot_id() {
    local pdb_file="$1"
    local chain="$2"
    local label="$3"
    
    [[ -f "$pdb_file" ]] || return 1
    
    # Try to extract from DBREF records first
    local uniprot_id
    uniprot_id=$(grep "^DBREF" "$pdb_file" | grep "UNP" | head -n1 | awk '{for(i=1;i<=NF;i++) if($i=="UNP") print $(i+1)}')
    
    echo "Extracted: '$uniprot_id'"
}

cat <<PDB > test.pdb
HEADER    TEST PDB 1A00
DBREF  1A00 A    1   100  UNP    P12345   P12345_HUMAN     1    100
ATOM      1  CA  ALA A 101      10.000  10.000  10.000  1.00  0.00           C
END
PDB

extract_uniprot_id test.pdb A LabelA
