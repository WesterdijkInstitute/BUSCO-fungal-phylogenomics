#! /usr/bin/env python

"""
Once busco has run on a set of assemblies, run this program to get a summary 
of everything
"""

import sys
import os
import argparse
from pathlib import Path
import pandas as pd
from zipfile import ZipFile

def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputfolder", help="Path to folder with busco\
        results (each result is a subfolder)", required=True, type=Path)
    parser.add_argument("-f", "--filter_list", help="Optional. File with list \
        of assemblies. Only the assemblies from the input folder that are in \
        this list will be processed (it can be a tab-separated file)", type=Path)
    return parser.parse_args()

if __name__ == "__main__":
    args = arg_parser()
    
    i = args.inputfolder
    if not i.is_dir():
        sys.exit("Error: given input folder not a folder")
        
    filter_list = set()
    if args.filter_list:
        try:
            with open(args.filter_list) as f:
                for line in f:
                    filter_list.add(line.split("\t")[0].strip())
            print("Option --filter_list used. Got {} assembly accessions".format(len(filter_list)))
        except IOError:
            sys.exit("Error: --filter_list used, but cannot open file")
    
    # Create a pandas data frame to store all the single copy busco hits
    all_busco_hits = set()
    data = dict()
    for num, folder in enumerate(sorted(i.glob("*"))):
        # only traverse folders
        if not folder.is_dir():
            continue
        
        assembly = folder.parts[-1]
        
        if filter_list:
            if assembly not in filter_list:
                continue
        
        # place with fasta files
        rfs = list(folder.glob("run_*"))
        assert(len(rfs) == 1)
        run_folder = rfs[0]

        target_folder = run_folder / "busco_sequences/single_copy_busco_sequences/"
        
        if not target_folder.is_dir():
            # try to see if results have been zipped
            try:
                scbs_set = set()
                target_zip = folder / "run_ascomycota_odb10/busco_sequences.zip"
                with ZipFile(target_zip) as z:
                    for x in z.namelist():
                        xpath = Path(x)
                        # only get complete single_copy_busco_sequences
                        if len(xpath.parts) < 2:
                            continue
                        if xpath.parts[1] != "single_copy_busco_sequences":
                            continue
                        if xpath.suffix == ".fna":
                            scbs_set.add(xpath.stem)
            except:
                sys.exit("Error with zip file {}".format(target_zip))
        else:
            scbs_set = set([fasta_file.stem for fasta_file in target_folder.glob("*.faa")])
            
        data[assembly] = scbs_set
        all_busco_hits.update(scbs_set)
            
        
    # got all data, now create dataframe...
    busco_table = pd.DataFrame(data=False, dtype=bool, index=sorted(data.keys()), columns=sorted(all_busco_hits))
    # ...and fill with values
    for assembly, hits_set in data.items():
        for hit in hits_set:
            busco_table.loc[assembly, hit] = True
    #print(busco_table.head())    
    
    # Finalize. 
    print("Checked {} result folders".format(num+1))
    busco_table.to_csv("./busco_a-p_matrix.tsv", sep="\t")
    
    
    
    
