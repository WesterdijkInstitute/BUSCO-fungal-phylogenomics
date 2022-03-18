#! /usr/bin/env python

"""
Input:
* Folder with BUSCO results
* List of Filtered Assemblies
* List of Target Genes

Extracts sequences of each Target Gene, for each Filtered Assemblies
"""

import sys
import os
import io
import argparse
from pathlib import Path
from zipfile import ZipFile, BadZipFile


def parameters_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--results", help="Base path with BUSCO results", \
        required=True, type=Path)
    parser.add_argument("-a", "--assemblies", help="List of assemblies from \
        which to extract the sequences (e.g. 'matrix_analysis_Top_0.XX_Assemblies.tsv')", \
        required=True, type=Path)
    parser.add_argument("-t", "--targetgenes", help="List of Target Genes (e.g. 'matrix_analysis_S_genes_in_0.XX_of_Top_0.XX_assemblies.tsv')", \
        required=True, type=Path)
    parser.add_argument("-o", "--outputfolder", help="Folder with unaligned sequence file. \
        Default: ./output/Target_Genes_unaligned", \
        type=Path, default="./output/Target_Genes_unaligned")
    parser.add_argument("--aa", help="Extract protein sequences instead of DNA",
        default=False, action="store_true")
    return parser.parse_args() 


def read_list(filepath):
    rset = set()
    
    try:
        with open(filepath) as f:
            for line in f:
                if line[0] == "#" or line.strip() == "" or line.startswith("Assembly"):
                    continue
                rset.add(line.strip().split("\t")[0])
    except IOError:
        sys.exit("Error: cannot open {}".format(filepath))
    else:
        return rset


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
    
    base_folder = args.results
    if not base_folder.is_dir():
        sys.exit("Error: {} not a folder".format(base_folder))
        
    TargetGenes = read_list(args.targetgenes)
    FilteredAssemblies = read_list(args.assemblies)
    
    o = args.outputfolder
    if not o.is_dir():
        os.makedirs(o, exist_ok=True)
        
    # Choose DNA or AA output
    suffix = "fna"
    file_type = "dna"
    if args.aa:
        suffix = "faa"
        file_type = "aa"
    
    # Pre-assemble all zip file names, and find out if any of them is missing
    target_zips = dict()
    not_found = set()
    for asm in FilteredAssemblies:
        asm_folder = base_folder / asm
        rfs = list(asm_folder.glob("run_*"))
        try:
            assert(len(rfs) == 1)
        except AssertionError:
            print(f"Can't find the run folder for assembly {asm}")
            continue
        run_folder = rfs[0]
        
        # get either the folder or its compressed version
        if (run_folder / "busco_sequences.zip").is_file():
            target_zips[asm] = run_folder / "busco_sequences.zip"
        elif (run_folder / "busco_sequences/single_copy_busco_sequences").is_dir():
            target_zips[asm] = run_folder / "busco_sequences/single_copy_busco_sequences/"
        else:
            not_found.add(asm)
            
    if not_found:
        print("Not found: {} BUSCO 'busco_sequences' zip".format(len(not_found)))
        sys.exit(", ".join(not_found))
    
    for gene in TargetGenes:
        with open(o / "{}.{}.fasta".format(gene, file_type), "w") as f:
            for asm in sorted(FilteredAssemblies - not_found):
                target = target_zips[asm]
                
                if target.is_file():
                    with ZipFile(target_zip) as z:
                        try:
                            with io.TextIOWrapper(z.open("busco_sequences/single_copy_busco_sequences/{}.{}".format(gene, suffix)), encoding="utf-8") as fasta:
                                old_header = str(fasta.readline())[1:].split(" ")[0]
                                new_header = ">{}_{} {}".format(gene, asm, old_header)
                                seq = "".join([str(l).strip() for l in fasta.readlines()])
                                
                                f.write("{}\n".format(new_header))
                                f.write(sequence80(seq))
                        # this assembly doesn't have a copy of this (S) gene
                        except KeyError:
                            f.write(">{}_{}\n".format(gene, asm))
                else:
                    try:
                        with open(target / f"{gene}.{suffix}") as fasta:
                            old_header = str(fasta.readline())[1:].split(" ")[0]
                            new_header = ">{}_{} {}".format(gene, asm, old_header.strip())
                            seq = "".join([str(l).strip() for l in fasta.readlines()])
                            
                            f.write("{}\n".format(new_header))
                            f.write(sequence80(seq))
                    except IOError:
                        f.write(">{}_{}\n".format(gene, asm))
