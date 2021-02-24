#! /usr/bin/env python

"""
Input:
* Folder with Target Gene fasta files

Uses MAFFT to make multiple sequence alignments of all Target Genes
It also concatenates all alignments and creates the RAxML parts file for IQ-TREE
"""

import sys
import os
import argparse
from pathlib import Path
from multiprocessing import Pool
from subprocess import run, PIPE


def parameters_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputfolder", help="Folder with unaligned \
        .fasta files", required=True, type=Path)
    parser.add_argument("-o", "--outputfolder", help="Folder with aligned \
        sequence file.", required=True, type=Path)
    parser.add_argument("-p", "--processes", help="Number of MAFFT instances. \
        Default: 1", type=int, default=1)
    parser.add_argument("-t", "--threads", help="--threads parameter for each \
        MAFFT process. Default: 2", type=int, default=1)
    return parser.parse_args()


def do_mafft(fasta, algn, thread):
    cmd = ["mafft", "--auto", "--quiet", "--thread", str(thread), str(fasta)]
    #cmd = ["mafft", "--auto", "--thread", str(thread), str(fasta)]
    print(" ".join(cmd))
    with open(algn, "w") as f:
        proc_mafft = run(cmd, encoding="utf-8", check=True, stdout=f)
    return


if __name__ == "__main__":
    pars = parameters_parser()
    
    i = pars.inputfolder
    if not i.is_dir():
        sys.exit("Error, {} not a folder".format(i))

    o = pars.outputfolder
    if not o.is_dir():
        os.makedirs(o, exist_ok=True)
        
    with Pool(pars.processes) as pool:
    #if True: # comment above and uncomment this for serialized processing
        for fasta in i.glob("*.fasta"):
            aligned_file = o / "{}.algn".format(fasta.stem)
            pool.apply_async(do_mafft, args=(fasta, aligned_file, pars.threads, ))
            #do_mafft(fasta, aligned_file, pars.threads) # comment above and uncomment this for serialized processing
        pool.close()
        pool.join()

