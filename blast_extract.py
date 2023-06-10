# -*- coding: utf-8 -*-
"""
Created on Wed May 17 13:49:29 2023

@author: 44741
"""

import click 
import pandas as pd 
from Bio import SeqIO 

@click.command()
@click.option("-f", "--fasta", required=True, help="path to fasta file")
@click.option("-b", "--blast_file", required=True, help="path to blast file, must be in standard outputfmt6")
@click.option("-o", "--output_file", required=True, help="output file name")

def __main__(fasta, blast_file, output_file):
    blastdf = pd.read_csv(blast_file, sep="\t", names = ["qseqid", "sseqid", "pident", "length", "mismatch","gapopen","qstart","qend","sstart","send","evalue", "bitscore"])
    resorted_df = blastdf.sort_values(by=["evalue","bitscore", "pident", "length"], ascending=[True, False, False, False])
    no_dups = resorted_df.drop_duplicates(subset = "sseqid", keep = "first", inplace = False)
    blast_results = no_dups.reset_index() #loads all of your blast results, sorts them so the best hits are at the top and picks one hit per sequence 
    
    records = list(SeqIO.parse(fasta,"fasta"))
    with open(output_file, "a") as f:
        for r in records:
            for i in range(len(blast_results)):
                if r.id == blast_results.loc[i, "sseqid"]:
                    SeqIO.write(r, f, "fasta")
                    
__main__()