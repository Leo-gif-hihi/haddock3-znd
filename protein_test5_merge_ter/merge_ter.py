def merge_pdb_chains(input_filename, output_filename, new_chain_id='B'):
    with open(input_filename, 'r') as f_in, open(output_filename, 'w') as f_out:
        
        atom_serial = 0
        current_res_seq = 0
        # specific tracking to handle residue numbering changes across chains
        prev_res_id = None # (chain, res_seq, ins_code)
        
        # New variable to track the chain from the previous atom
        last_seen_chain = None

        for line in f_in:
            if line.startswith(('ATOM', 'HETATM')):
                # Parse current residue identity and chain
                old_chain = line[21]
                old_res_seq = line[22:26]
                old_ins_code = line[26]
                current_res_id = (old_chain, old_res_seq, old_ins_code)
                
                # --- CHANGE START: Insert TER if chain changes ---
                # If we have seen a chain before, and the current one is different
                if last_seen_chain is not None and old_chain != last_seen_chain:
                    f_out.write("TER\n")
                
                # Update the last seen chain to the current one
                last_seen_chain = old_chain
                # --- CHANGE END ---

                # Check if we moved to a new residue
                if current_res_id != prev_res_id:
                    current_res_seq += 1
                    prev_res_id = current_res_id
                
                atom_serial += 1
                
                # Reconstruct the line
                new_line = (
                    line[:6] +                  # Record name
                    f"{atom_serial:5d}" +       # New Atom Serial
                    line[11:21] +               # Atom Name, Res Name
                    new_chain_id +              # New Chain ID
                    f"{current_res_seq:4d}" +   # New Res Seq
                    " " +                       # Insertion Code (cleared)
                    line[27:]                   # Coordinates and rest
                )
                f_out.write(new_line)
                
            elif line.startswith('END'):
                pass
            elif line.startswith('TER'):
                pass
            else:
                if line.startswith(('CRYST1', 'SCALE', 'REMARK')):
                    f_out.write(line)

        # Write final TER and END
        f_out.write("TER\n")
        f_out.write("END\n")

    print(f"Successfully created {output_filename} with {atom_serial} atoms.")

# Usage
input_pdb = 'P7KBJ_0_P7KBJ_combined.pdb'
output_pdb = 'P7KBJ_merged.pdb'
merge_pdb_chains(input_pdb, output_pdb)