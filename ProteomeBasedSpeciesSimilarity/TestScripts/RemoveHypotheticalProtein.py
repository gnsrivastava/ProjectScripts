##Author: Gopal Srivastava
## The script uses a fasta file as input and removes the hypothetical or undefined proteins from the fasta file

from Bio import SeqIO
import re
import os, sys, argparse

def parseFasta(inputfile, type):
  sequences = SeqIO.parse(f'{inputfile}', 'fasta')
  filtered = [seq for seq in sequences if f'{type}' not in seq.description]

  with open(f'{inputfile}', 'wt') as output:
    SeqIO.write(filtered, output, 'fasta')

def main():
  ap = argparse.ArgumentParser(description="Fetch protein FASTA from BV-BRC genome_feature with filters.")
  ap.add_argument("--file", "-f", help="File with one ID per line")
  ap.add_argument("--type", "-t", help="Either use hypothetical or undefined as input.")
  args = ap.parse_args()

  if args.type:
    type=args.type
    
  if args.file:
    inputfile = args.file
    parseFasta(inputfile, type)

if __name__ == "__main__":
  main()
