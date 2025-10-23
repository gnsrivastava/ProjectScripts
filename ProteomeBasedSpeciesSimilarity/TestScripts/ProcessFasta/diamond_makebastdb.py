import os
import glob
import subprocess
import pandas as pd
# replace location to installation with your diamond path

def create_diamond_db(input_dir: str, output_dir: str = 'NewBVBRC_db'):
    """
    Creates a Diamond protein database from all .faa (protein FASTA) 
    files in the specified directory.
    
    Args:
        input_dir (str): Path to the folder containing the strain .faa files.
    """
    
    # Find all protein FASTA files (.faa)
    protein_files = glob.glob(os.path.join(input_dir, "*.faa"))
        
    if not protein_files:
        print(f"Error: No protein FASTA files (.faa) found in {input_dir}")
        return
    created_dbs = []
    for fasta_path in protein_files:
        # Determine the base name for the database (e.g., 'StrainA' from 'StrainA.faa')
        fasta_basename = os.path.basename(fasta_path)
        db_name = os.path.splitext(fasta_basename)[0] 

        print(f"\t-->Creating database for: {fasta_basename}")
        try:
            command = [
                "/location to installation/diamond", "makedb", 
                "--in", fasta_path, 
                "-d", os.path.join(output_dir, db_name)
            ] 
            # Execute the command
            # Note: subprocess.run is typically used, but the original logic is preserved.
            subprocess.run(command, check=True, capture_output=True, text=True)
            print(f"  -> Success: {db_name}.dmnd created.")
            created_dbs.append(db_name)
        except subprocess.CalledProcessError as e:
            print(f"  -> Error running Diamond makedb for {fasta_basename}: {e.stderr}")
        except FileNotFoundError:
            print("Error: Diamond command not found. Make sure Diamond is in your PATH.")
            return []
    return created_dbs

if __name__=='__main__':
    dbs = create_diamond_db('./NewBVBRC')
