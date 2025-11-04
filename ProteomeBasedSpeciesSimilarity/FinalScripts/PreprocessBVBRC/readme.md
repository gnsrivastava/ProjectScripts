# Fasta cleaning pipeline before Diamond BLASTP

## Prerequisites
1. Conda
2. Biopython
3. Python v13
4. CD-HIT

**Step 1: Remove peptides from the original BV-BRC sequences** 

```
sh Remove_peptides_parallel.sh
```
- It will create slrum_scripts and chunk folders with required SLRUM scripts and input files
- Use command below to run SLRUM script

```
for file in slrum_scripts/*.sh
  do
    sbatch $file; 
    sleep 5; # Sleep is to ensure that jobs do not get blocked by the servers
  done; 
```
After each step we will get set of slrum_scripts and chunks folder and thus we need to ensure that we ```delete``` the ```slrum_scripts``` and ```chunks``` folder after each step. 
And for each step job submission step will remain same as ```Step1```. 

**Step 2: Remove hypothetical proteins from the Step 1 BV-BRC sequences**
Here we do not remove **undefined** as we are only using STRING protein IDs for PPIs 
```
sh Remove_Hypothetical_proteins.sh
``` 
**Step 3: Run CD-HIT**
- On the files with at least 1 sequence, run `CDHIT` with ```-c 0.95``` to ensure only proteins with extremely high similarity get removed.
- Use following command to run ```CD-HIT```.
```
# Before running the script ensure that you have a cdhit conda environment in your server or system
sh Remove_redundant_proteins.sh
```
Step 4: Run diamond on these non-redundant sequences against STRING (Script: 
Step 5: Run Hungarian Assignment (Script:
