#! /usr/bin/env python

"""
This script creates a "summary of summaries" of all busco results' short summaries.

It also compares the number of "Complete and single-copy BUSCOs" reported in the 
short summary against the actual number of files to detect possible truncated results
Counted result files are dna sequence files (fna)
"""

import sys
import os
import argparse
from pathlib import Path
#import zipfile as zip
from zipfile import ZipFile, BadZipFile

def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--buscofolders", help="Path to folder with busco\
        results (each result is a subfolder)", required=True, type=Path)
    parser.add_argument("-m", "--metadata", help="A tab-separated file in which \
        first column is the assembly name, third column is the species name, and\
        fourth column is the strain name. Will be used to annotate summary file \
        (Optional).", type=Path)
    return parser.parse_args()


def read_metadata_file(filepath):
    name_dictionary = dict()
    
    # user didn't use parameter
    if not filepath:
        return name_dictionary
    
    # user used parameter but argument doesn't point to valid file
    if not filepath.is_file():
        return name_dictionary
    
    try:
        with open(filepath) as f:
            for line in f:
                if line[0] == "#" or line.strip() == "":
                    continue
                x = line.split("\t")
                
                if len(x) != 4:
                    print(x)
                    sys.exit()
                
                name = x[2]
                if x[3].strip() != "":
                    name = "{} {}".format(name, x[3].strip())
                name_dictionary[x[0]] = name
    except IOError:
        return name_dictionary
    else:
        return name_dictionary

if __name__ == "__main__":
    args = arg_parser()
    
    i = args.buscofolders
    if not i.is_dir():
        sys.exit("Error: given input folder not a folder")
        
    name_dictionary = read_metadata_file(args.metadata)
    
    # Create a pandas data frame to store all the single copy busco hits
    discrepancies = dict() # key: assembly. value: tuple of reported complete+single BUSCOs (S), number of files
    summary = dict() # value: a list of [C, S, D, F, M] BUSCO results
    for num, folder in enumerate(sorted(i.glob("*"))):
        # only traverse folders
        if not folder.is_dir():
            continue
        
        assembly = folder.parts[-1]
        
        # place with summary
        summary_file = folder / "run_ascomycota_odb10/short_summary.txt"
        if not summary_file.is_file():
            print("Warning: Can't find summary file for assembly {}".format(assembly))
            continue
        with open(summary_file) as f:
            lines = f.readlines()
            C = int(lines[8].strip().split("\t")[0])
            S = int(lines[9].strip().split("\t")[0])
            D = int(lines[10].strip().split("\t")[0])
            F = int(lines[11].strip().split("\t")[0])
            M = int(lines[12].strip().split("\t")[0])
            summary[assembly] = [C, S, D, F, M]
        
        # place with fasta files
        target_folder = folder / "run_ascomycota_odb10/busco_sequences/single_copy_busco_sequences/"
        if not target_folder.is_dir():
            # try to see if results where zipped            
            try:
                target_zip = folder / "run_ascomycota_odb10/busco_sequences.zip"
                fnas = set()
                faas = set()
                with ZipFile(target_zip) as z:
                    for x in z.namelist():
                        if x.startswith("busco_sequences/single_copy_busco_sequences/"):
                            if x[-3:] == "fna":
                                fnas.add(x)
                            elif x[-3:] == "faa":
                                faas.add(x)
                            else:
                                print("Unknown type {}".format(x))
                #sys.exit()
                    
            except:
                sys.exit("Error: Can't find results for assembly {}".format(target_folder))
                
        else:
            fnas = set([fasta_file.stem for fasta_file in target_folder.glob("*.fna")])
            faas = set([fasta_file.stem for fasta_file in target_folder.glob("*.faa")])
            
        if len(fnas) != len(faas):
            print("Warning! BUSCO results for {} have different number of fna and faa files ({}, {})".format(assembly, len(fnas), len(faas)))
            
        # Detect discrepancy between reported number of complete and single-copy hits vs actual files
        if len(fnas) != S:
            discrepancies[assembly] = (S, len(fnas))
            
    
    
    # Finalize. 
    print("Checked {} result folders".format(num+1))
    
    with open("./busco_set_results_summary.tsv", "w") as f:
        # write header
        f.write("Assembly\t[C]omplete BUSCOs\tComplete and [S]ingle-copy BUSCOs\tComplete and [D]uplicated BUSCOs\t[F]ragmented BUSCOs\t[M]issing BUSCOs\tName\n")
        for assembly, numbers in summary.items():
            try:
                name = name_dictionary[assembly]
            except KeyError:
                name = ""
            f.write("{}\t{}\t{}\n".format(assembly, "\t".join([str(n) for n in numbers]), name))
    
    if len(discrepancies) > 0:
        print("Found {} assemblies with discrepancies".format(len(discrepancies)))
        with open("./busco_set_results_discrepancies.tsv", "w") as f:
            f.write("Assembly\tSinge copy reported\tSingle copy found\n")
            for assembly in discrepancies:
                S, got = discrepancies[assembly]
                f.write("{}\t{}\t{}\n".format(assembly, S, got))
    else:
        print("No discrepancies were found!")
    
    
    
 
