import pandas as pd 
from Bio import SeqIO
from Bio import SeqUtils 
import click

@click.command()
@click.option("--input_file", "-i", required=True)
@click.option("--output_file", "-o", default="GCcontent.txt", required=False)
@click.option("--window_size", "-w", defualt=0, required=False)

def gcContent(input_file, output_file, window_size):
    records = list(SeqIO.parse(input_file,"fasta"))
    data = pd.DataFrame(columns = ["Scaffold", "Start", "End", "GC"])
    if window_size == 0:
        for r in records: 
            x = len(data)
            seq = r.seq
            GC = SeqUtils.GC(seq)
            data.loc[x,"Scaffold"] = r.id
            data.loc[x, "Start"] = 1
            data.loc[x, "End"] = len(seq)
            data.loc[x, "GC"] = GC
    else:
        for r in records: 
            seq = r.seq
            for i in range(0, len(seq), window_size):
                x = len(data)
                chunk = seq[i:i+window_size]
                GC = SeqUtils.GC(chunk)
                data.loc[x,"Scaffold"] = r.id
                data.loc[x, "Start"] = i
                data.loc[x, "End"] = i+window_size
                data.loc[x, "GC"] = GC
    data.to_csv(output_file, sep="\t", index=False)

gcContent()
