#! /usr/bin/env python

"""
Input:
* Folder with aligned (and trimed) fasta files (.algn extension)

concatenates all alignments and creates a nexus file for iqtree
Expects headers to be in the format
>[BUSCO id]_[assembly acc.] [optional: description]
"""

import sys
import os
import io
import argparse
from collections import defaultdict
from pathlib import Path


def parameters_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputfolder", help="Folder with aligned sequences", \
        required=True, type=Path)
    parser.add_argument("-n", "--name", help="Base name for the output", \
        type=str, required=True)
    return parser.parse_args()


def read_fasta(fasta: Path, data: dict) -> str:
    """
    fasta: Path object to each fasta file
    data: dictionary where key: acc, value: a list of sequences

    All headers should have the format
    >[BUSCO id]_[assembly acc.] [optional: description]
    """

    gene_id = ""

    sequence_lengths = set() # use to detect any possible variation
    with open(fasta) as f:
        header = ""
        acc = ""
        sequence = list()
        
        for line in f:
            if line.strip() == "":
                continue
            
            if line[0] == ">":
                # First time; we have no accession yet
                if header != "":
                    # there shouldn't be headers without sequences; 
                    # even missing genes should have deletion characters
                    assert(sequence)

                    data[acc].append("".join(sequence))
                    sequence_lengths.add(len(data[acc]))
                    sequence = list()
                
                header = line[1:].split(" ")[0].strip() # get rid of description
                if gene_id == "":
                    gene_id = header.split("_")[0]
                elif gene_id != header.split("_")[0]:
                        sys.exit("Error! found a different gene id in {}: {}, {}".format(fasta, gene_id, header.split("_")[0]))
                acc = "_".join(header.split("_")[1:])
                    
            else:
                sequence.append(line.strip())

        # record the last sequence
        assert(sequence)
        data[acc].append("".join(sequence))
        sequence_lengths.add(len(data[acc]))


    if len(sequence_lengths) != 1:
        sys.exit("Error! found {} different sequence lengths in {}".format(len(sequence_lengths), fasta))

    return gene_id


def sequence80(seq):
    """
    Re-formats a sequence to limit to 80 columns
    """
    length = len(seq)
    part_one = "\n".join([seq[row*80:(row+1)*80] for row in \
            range(length // 80)])
        
    part_two = ""
    remainder = length % 80
    if remainder > 0:
        part_two = "{}\n".format(seq[-remainder:])
    return "{}\n{}".format(part_one, part_two)


if __name__ == "__main__":
    args = parameters_parser()
    
    i = args.inputfolder
    n = args.name
    
    if not i.is_dir():
        sys.exit("Error: given input folder is not a valid folder")
        
    # three output files concatenated seqs. nexus
    data = defaultdict(list)
    geneids = list()
    for fasta in sorted(i.glob("*.algn")):
        geneids.append(read_fasta(fasta, data))


    with open(f"{args.name}.fasta", "w") as f, open(f"{args.name}.nex", "w") as n:
        concatenation_lengths = set() # to verify all accs. have all positions
        for acc, seqs in data.items():
            f.write(f">{acc}\n")
            concatenation = "".join(seqs)
            f.write(sequence80(concatenation))
            concatenation_lengths.add(len(concatenation))

        n.write("#nexus\n")
        n.write("begin sets;\n")
        current_position = 1
        for i in range(len(data[acc])):
            gene_id = geneids[i]
            seq = data[acc][i]
            n.write(f"\tcharset\t{gene_id} = {current_position}-{current_position+len(seq)-1};\n")
            current_position += len(seq)
        n.write("end;")