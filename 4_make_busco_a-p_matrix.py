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
    return parser.parse_args()

if __name__ == "__main__":
    args = arg_parser()
    
    i = args.inputfolder
    if not i.is_dir():
        sys.exit("Error: given input folder not a folder")
    
    # Create a pandas data frame to store all the single copy busco hits
    all_busco_hits = set()
    data = dict()
    for num, folder in enumerate(sorted(i.glob("*"))):
        # only traverse folders
        if not folder.is_dir():
            continue
        
        assembly = folder.parts[-1]
        
        # place with fasta files
        target_folder = folder / "run_ascomycota_odb10/busco_sequences/single_copy_busco_sequences/"
        
        if not target_folder.is_dir():
            # try to see if results have been zipped
            try:
                scbs_set = set()
                target_zip = folder / "run_ascomycota_odb10/busco_sequences.zip"
                with ZipFile(target_zip) as z:
                    for x in z.namelist():
                        if x[-3:] == "fna":
                            scbs_set.add(Path(x).stem)
            except:
                sys.exit("Error with target folder {}".format(target_folder))
        else:
            scbs_set = set([fasta_file.stem for fasta_file in target_folder.glob("*.fna")])
            
        data[assembly] = scbs_set
        all_busco_hits.update(scbs_set)
            
        
    # got all data, now create dataframe...
    busco_table = pd.DataFrame(data=False, dtype=bool, index=sorted(data.keys()), columns=sorted(all_busco_hits))
    # ...and fill with values
    for assembly, hits_set in data.items():
        for hit in hits_set:
            busco_table.loc[assembly, hit] = True
    print(busco_table.head())    
    
    # Finalize. 
    print("Checked {} result folders".format(num+1))
    busco_table.to_csv("./busco_a-p_matrix.tsv", sep="\t")
    
    
    
    
