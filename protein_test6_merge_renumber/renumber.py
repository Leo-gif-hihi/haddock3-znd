import sys

def renumber_pdb_continuously(input_pdb, output_pdb):
    """
    Renumbers residues in a PDB file continuously across all chains.
    Handles short TER lines to prevent index errors.
    """
    global_res_num = 0
    prev_res_id = None # Stores (chain, res_seq, i_code)
    
    try:
        with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
            for line in f_in:
                # Process ATOM and HETATM records
                if line.startswith(('ATOM', 'HETATM')):
                    # Ensure line is long enough to parse
                    if len(line) < 27:
                        f_out.write(line)
                        continue

                    # Parse identity: Chain(21), ResNum(22-26), InsCode(26)
                    chain_id = line[21]
                    res_seq = line[22:26]
                    i_code = line[26]
                    
                    current_res_id = (chain_id, res_seq, i_code)
                    
                    # Increment counter only if we hit a new residue ID
                    if current_res_id != prev_res_id:
                        global_res_num += 1
                        prev_res_id = current_res_id
                    
                    # Format new number (4 chars, right aligned)
                    new_res_str = "{:>4}".format(global_res_num)
                    
                    # Reconstruct line: chars 0-22 + new_num + chars 26-end
                    new_line = line[:22] + new_res_str + line[26:]
                    f_out.write(new_line)

                # Process TER records
                elif line.startswith('TER'):
                    # If TER has residue info (long line), update it to match previous atom
                    if len(line) >= 27:
                        new_res_str = "{:>4}".format(global_res_num)
                        new_line = line[:22] + new_res_str + line[26:]
                        f_out.write(new_line)
                    else:
                        # Short TER line (just "TER"), write as is
                        f_out.write(line)
                
                # Write all other lines (HEADER, REMARK, END, etc.) as is
                else:
                    f_out.write(line)
            
        print(f"Successfully created '{output_pdb}' (Last residue number: {global_res_num})")

    except FileNotFoundError:
        print(f"Error: The file '{input_pdb}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# --- Execute ---
# Make sure this filename matches your uploaded file exactly
input_filename = 'P7KBJ_0_P7KBJ_combined.pdb'
output_filename = 'P7KBJ_combined_renumbered.pdb'

renumber_pdb_continuously(input_filename, output_filename)