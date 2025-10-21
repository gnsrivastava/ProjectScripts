import pandas as pd
import numpy as np 
from scipy.optimize import linear_sum_assignment 
from Bio import SeqIO 
import os
import sys
from mpi4py import MPI
import subprocess
from itertools import product
from typing import List, Tuple, Any, Dict

# Initialize MPI
COMM = MPI.COMM_WORLD
RANK = COMM.Get_rank()
SIZE = COMM.Get_size()

# Global variable to store the mapping from Genome ID to Taxon ID (BVBRC)
# Type hint for the global variable
BVBRC_GENOME_TO_TAXID: Dict[str, str] = {}


# --- Helper Functions (Unchanged from previous versions) ---

def list_files_low_memory(folder_path: str) -> str:
    """Yields file names (base name without extension) from a folder using low memory."""
    try:
        if not os.path.isdir(folder_path):
            if RANK == 0:
                print(f"[Warning] Folder not found at {folder_path}. Proceeding.")
            return
        for entry in os.scandir(folder_path):
            if entry.is_file():
                yield os.path.splitext(entry.name)[0]
    except Exception as e:
        if RANK == 0:
            print(f"[Error] An error occurred in list_files_low_memory: {e}")
        sys.exit(1)

def Species(data: pd.DataFrame, ColumnID: str) -> pd.Series:
    """Extracts species name (first two words) from a column."""
    if ColumnID in data.columns:
        return (data[f'{ColumnID}'].astype(str).str.split(' ').apply(lambda x: " ".join(x[:2])))
    return pd.Series([None] * len(data))
    
def Subset(data: pd.DataFrame, stringdf: pd.DataFrame, genomesLST: List[str]) -> pd.DataFrame:
    """Filters BVBRC data based on species in STRING and available FASTA files."""
    data['Genome ID'] = data['Genome ID'].astype('str')
    data = data[data.Species.isin(stringdf.Species.unique())]
    data = data[data['Genome ID'].isin(genomesLST)]
    return data

# --- Diamond Execution Function (Processes a single pair) ---

# NOTE: The BVBRC_GENOME_TO_TAXID map is now accessed directly within the function,
# having been correctly populated and broadcast to all ranks via the global variable.
def run_diamond_single_pair(pair: Tuple[str, str], QPATH: str, TPATH: str, OUTPATH: str) -> str:
    """Run DIAMOND blastp for a single query-target pair with hierarchical output."""
    global RANK, BVBRC_GENOME_TO_TAXID # Access the global map
    
    try:
        query_genome_id, target_tax_id = pair[0], str(pair[1])
        
        # 1. Determine the BVBRC Taxon ID for the query genome
        bvbrc_taxon_id = BVBRC_GENOME_TO_TAXID.get(query_genome_id)
        if not bvbrc_taxon_id:
            return f"[Rank {RANK}] ERROR: BVBRC Taxon ID not found for {query_genome_id} in the map. Skipping."
            
        # 2. Define the hierarchical output path
        # Output: OUTPATH/BVBRC_TAXON_ID/QUERY_ID_TARGET_ID.tsv
        species_out_dir = os.path.join(OUTPATH, bvbrc_taxon_id)
        out_file = f'{query_genome_id}_{target_tax_id}.tsv'
        full_out_path = os.path.join(species_out_dir, out_file)
        
        # Create the species-specific output directory
        # Since this runs on all ranks, only the ranks assigned tasks for a given
        # BVBRC Taxon ID will create that specific folder. This is safe with os.makedirs(exist_ok=True)
        os.makedirs(species_out_dir, exist_ok=True)

        # 3. Define DIAMOND command
        cmd = [
            "/project/gsriva2/diamond", "blastp",
            "-q", f'{QPATH}/{query_genome_id}_nonredundant.faa',
            "-d", f'{TPATH}/{target_tax_id}', 
            "-o", full_out_path,
            "-f", "6", "-k", "1", "--evalue", "1e6", "--header"
        ]
        
        print(f"[Rank {RANK}] Submitting: {' '.join(cmd)}")
        
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        
        return full_out_path 
        
    except subprocess.CalledProcessError as e:
        return f"[Rank {RANK}] ERROR: DIAMOND failed for pair {pair}. Exit code: {e.returncode}"
    except Exception as e:
        return f"[Rank {RANK}] ERROR: An unexpected error occurred for pair {pair}: {e}"


# --- Main MPI Workflow ---
def run_main_workflow_parallel_diamond(QPATH: str, TPATH: str, OUTPATH: str):
    
    global BVBRC_GENOME_TO_TAXID # Declare intent to modify the global variable
    
    all_diamond_tasks = []
    
    # ----------------------------------------------------
    # PHASE 1: Data Loading & Task Generation (Only on Rank 0)
    # ----------------------------------------------------
    
    if RANK == 0:
        print("[Rank 0] Starting data loading and task generation...")
        
        # 1. Load BVBRC instances
        df = pd.read_csv('LabTested_BVBRC_Instances.csv', low_memory=False)
        df['Genome ID'] = df['Genome ID'].astype('str')
        df['Species'] = Species(df, 'Genome Name')
        
        # Populate BVBRC_GENOME_TO_TAXID
        if 'Taxon ID' not in df.columns:
             print("[Rank 0] ERROR: 'Taxon ID' column not found in BVBRC file. Cannot create hierarchical output.")
             sys.exit(1)
             
        # Populate the global map using 'Genome ID' and 'Taxon ID' columns
        BVBRC_GENOME_TO_TAXID = dict(zip(df['Genome ID'], df['Taxon ID'].astype(str)))
        print(f"[Rank 0] Populated BVBRC map with {len(BVBRC_GENOME_TO_TAXID)} entries.")
        
        # 2. Load String bacteria dataset
        string = pd.read_csv('StringBacteria.tsv', sep='\t')
        string['Species'] = Species(string, 'STRING_name_compact')
        
        # 3. Filter data
        genomesLST = list(list_files_low_memory('../Fasta/fastaBVBRC_org'))
        df_filtered = Subset(df, string, genomesLST)
        # Ensure 'Taxon ID' is kept for filtering/lookup sanity, though the global map holds the core info
        df_filtered = df_filtered[['Genome ID', 'Genome Name', 'Species', 'Taxon ID']].drop_duplicates()
        
        # 4. Generate ALL DIAMOND tasks (Cartesian Product per species)
        for species in df_filtered['Species'].unique():
            df_ids = df_filtered[df_filtered.Species==species]['Genome ID'].unique()
            string_ids = string[string.Species==species]['#taxon_id'].unique()
            all_diamond_tasks.extend(list(product(df_ids, string_ids)))

        print(f"[Rank 0] Generated {len(all_diamond_tasks)} total DIAMOND tasks.")
        
        # 5. Calculate workload distribution
        N_diamond = len(all_diamond_tasks)
        avg_d, rem_d = divmod(N_diamond, SIZE)
        counts_d = [avg_d + 1 if i < rem_d else avg_d for i in range(SIZE)]
        displs_d = [sum(counts_d[:i]) for i in range(SIZE)]
        
        task_slices_d = []
        for i in range(SIZE):
            start = displs_d[i]
            end = displs_d[i] + counts_d[i]
            task_slices_d.append(all_diamond_tasks[start:end])
    else:
        task_slices_d = None
        
    # 6. Scatter the DIAMOND tasks
    local_diamond_pairs = COMM.scatter(task_slices_d, root=0)
    
    # 7. Broadcast the BVBRC Genome ID -> Taxon ID map
    # This sends the populated dictionary from Rank 0 to all other ranks, 
    # updating their global BVBRC_GENOME_TO_TAXID variable.
    BVBRC_GENOME_TO_TAXID = COMM.bcast(BVBRC_GENOME_TO_TAXID, root=0)

    # ----------------------------------------------------
    # PHASE 2: Parallel DIAMOND Runs (All ranks)
    # ----------------------------------------------------
    
    if local_diamond_pairs:
        for pair in local_diamond_pairs:
            # The run_diamond_single_pair function now relies on the global map
            result = run_diamond_single_pair(pair, QPATH, TPATH, OUTPATH)
            print(result) 

    # Wait for all DIAMOND runs to complete before exiting
    COMM.barrier() 

    if RANK == 0:
        print("\n[Rank 0] === Parallel DIAMOND Runs Completed. ===")


# --- Execution Block ---
if __name__ == '__main__':
    # Define paths
    QPATH = '../Fasta/removed_redundant'
    TPATH = './StringSeqsDiamond'      # Path to the DIAMOND DBs
    OUTPATH = './OUTPUT'               # Base directory for hierarchical DIAMOND TSVs

    # Ensure the base output directory exists on Rank 0
    if RANK == 0:
        os.makedirs(OUTPATH, exist_ok=True)
    
    # Synchronize all processes before starting the main work
    COMM.barrier() 
    
    run_main_workflow_parallel_diamond(QPATH, TPATH, OUTPATH)
