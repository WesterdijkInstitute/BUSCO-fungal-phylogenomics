#! /usr/bin/env python

"""
Input:
* Folder with Target Gene fasta files
Output:
* A folder with aligned fasta files
* A folder with trimmed alignment files

Uses MAFFT to make multiple sequence alignments of all Target Genes
It also curates the alignments using trimal "gappyout" strategy
"""

import sys
import os
import argparse
from pathlib import Path
from multiprocessing import Pool
from subprocess import run, PIPE, DEVNULL


def parameters_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputfolder", help="Folder with unaligned \
        .fasta files", required=True, type=Path)
    parser.add_argument("-a", "--alignedfolder", help="Folder for aligned \
        sequence files.", required=True, type=Path)
    parser.add_argument("-t", "--trimmedfolder", help="Folder for trimmed \
        aligned sequence files (using command 'trimal -gappyout')", 
        required=True, type=Path)
    parser.add_argument("-p", "--processes", help="Number of MAFFT instances. \
        Default: 1", type=int, default=1)
    parser.add_argument("--threads", help="--threads parameter for each \
        MAFFT process. Default: 2", type=int, default=1)
    return parser.parse_args()


def do_mafft(fasta, algn, trimmed_algn, thread):
    cmd = ["mafft", "--auto", "--quiet", "--thread", str(thread), str(fasta)]
    print(" ".join(cmd))
    
    try:
        with open(algn, "w") as f:
            proc_mafft = run(cmd, encoding="utf-8", check=True, stdout=f)
    except CalledProcessError:
        print("Error running '{}'".format(" ".join(cmd)))
    else:
        # launch trimal
        # Keep empty sequences, keep original headers
        cmd = ["trimal", "-keepseqs", "-keepheader", "-in", str(algn), "-out", 
               str(trimmed_algn), "-gappyout"]
        print("\t{}".format(" ".join(cmd)))
        # NOTE: use capture_output=True but don't use it; to suppress warnings
        # about empty sequences
        proc_trimal = run(cmd, check=True, capture_output=True)
    return


if __name__ == "__main__":
    pars = parameters_parser()
    
    i = pars.inputfolder
    if not i.is_dir():
        sys.exit("Error, {} not a folder".format(i))

    o = pars.alignedfolder
    if not o.is_dir():
        os.makedirs(o, exist_ok=True)
        
    t = pars.trimmedfolder
    if not t.is_dir():
        os.makedirs(t, exist_ok=True)
        
    with Pool(pars.processes) as pool:
    #if True: # comment above and uncomment this for serialized processing
        for fasta in i.glob("*.fasta"):
            aligned_file = o / "{}.algn".format(fasta.stem)
            trimmed_algn = t / "{}.trimal.algn".format(fasta.stem)
            pool.apply_async(do_mafft, args=(fasta, aligned_file, trimmed_algn, pars.threads, ))
            #do_mafft(fasta, aligned_file, trimmed_algn, pars.threads) # comment above and uncomment this for serialized processing
        pool.close()
        pool.join()

