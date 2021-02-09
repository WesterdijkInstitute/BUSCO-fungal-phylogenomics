#! /usr/bin/env python

"""
STEP 1
given a json file with assembly data downloaded with NCBI datasets, download 
and update assembly data.

Get the JSON file with (using version 10.21.0 of ncbi-datasets)
ncbi-datasets summary genome taxon 147537 | python -m json.tool > <name>.json

The script is designed to be able to update the folder with downloaded assemblies when using newer json files
"""

import sys
import os
import json
import argparse
from pathlib import Path
import subprocess
from subprocess import STDOUT
import zipfile as zip

def command_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--json", help="JSON file downloaded with NCBI \
                        datasets with assembly data", required=True, 
                        type=Path)
    parser.add_argument("-t", "--tries", help="Try this many times to get each\
                        accession file. If it still fails, a warning will be \
                        issued. Default: 3 times", type=int, default=3)
    parser.add_argument("-o", "--outputfolder", help="Folder where the \
                        assemblies will be stored (default: assemblies/)", 
                        type=Path, default=(Path(__file__).parent/"assemblies"))
    parser.add_argument("-n", help="Optional: Number of assemblies to download (in the\
                        order of the JSON file)", type=int)
    
    return parser.parse_args()


def download_accession(acc, tries, outputfolder):
    zipfilename = outputfolder / (acc + ".zip")
    
    cmd = ["ncbi-datasets", "download", "genome", "accession"]
    #cmd.append("ncbi-datasets")
    #cmd.append("download")
    #cmd.append("assembly")
    cmd.append(acc)
    cmd.append("--filename")
    cmd.append(str(zipfilename))
    #print(" ".join(cmd))
    
    for n in range(tries):
        proc = subprocess.run(cmd, stderr=STDOUT, encoding="utf-8")
        try:
            proc.check_returncode()
        except subprocess.CalledProcessError:
            print("Error downloading {} (try {}/{})".format(acc, n, tries))
            print(proc.stderr)
        else:
            zipfilename = outputfolder / (acc + ".zip")
            # check here if file already exists (it should)
            if zipfilename.is_file():
                # and whether it can be unzipped
                try:
                    z = zip.ZipFile(zipfilename, "r")
                except zip.BadZipFile:
                    print("Error with zip file for {}".format(acc))
                else:
                    z.close()
                    return True
    print("Could not get {}. Try again (and possibly increase number of tries)".format(acc))
    return False


if __name__ == "__main__":
    options = command_parser()
    
    json_file = options.json
    if not json_file.is_file():
        sys.exit("JSON file is not a valid file")
        
    o = options.outputfolder
    if not o.is_dir():
        os.makedirs(o, exist_ok=True)
    
    accession_metadata_summary = list()
    with open(json_file) as f:
        j = json.load(f)
        
        for n, item in enumerate(j["assemblies"]):
            if options.n:
                if n+1 > options.n:
                    break
                
            asm = item["assembly"]
            
            asm_ac = asm["assembly_accession"]
            org = asm["org"]
            sci_name = org["sci_name"]
            strain = org.get("strain", "")
            tax_id = org.get("tax_id", "")
            
            zipfilename = o / (asm_ac + ".zip")
            # check here if file already exists
            if zipfilename.is_file():
                # and whether it can be unzipped
                try:
                    z = zip.ZipFile(zipfilename, "r")
                except zip.BadZipFile:
                    print(" Warning: Zip error for {}, re-downloading")
                else:
                    z.close()
                    accession_metadata_summary.append((asm_ac, tax_id, sci_name, strain))
                    continue
            
            if download_accession(asm_ac, options.tries, o):
                accession_metadata_summary.append((asm_ac, tax_id, sci_name, strain))
            
            
    with open(o/"metadata.tsv", "w") as f:
        f.write("\n".join(["\t".join(x) for x in accession_metadata_summary]))
            
            
            
