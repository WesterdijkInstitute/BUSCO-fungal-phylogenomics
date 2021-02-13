#! /usr/bin/env python

"""
STEP 2
given a folder with assemblies (as downloaded by using NCBI's datasets),
read the genome file to make busco analysis

Needs a busco database

Changes in version 2: now reads a file with a list of GCAs (from possibly
failed analysis) and doesn't skip those when launching busco

change in version 3: add a filter list. Only GCAs from inputfolder in the filter
list will be analized. When busco is done, zip some folders within the results to 
avoid millions of tiny files

"""

import sys
import os
import argparse
from pathlib import Path
import subprocess
from subprocess import STDOUT
from zipfile import ZipFile, BadZipFile, is_zipfile, ZIP_DEFLATED
from multiprocessing import Pool, cpu_count
import tempfile
import io
from shutil import rmtree

def command_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputfolder", help="Folder with zipped\
        assemblies", required=True, type=Path)
    parser.add_argument("-o", "--outputfolder", help="Folder with BUSCO\
        results (Default: 'Busco_results')", type=Path,
        default=(Path(__file__).parent/"Busco_results"))
    parser.add_argument("-d", "--dbfolder", help="Folder with a BUSCO \
        database", required=True, type=Path)
    parser.add_argument("--re_analyze_file", help="A file containing 1-assembly\
        accession per line. Those assemblies will be re-analized", type=Path, \
        default=None)
    parser.add_argument("--filter_list", type=Path, help="Read a txt file with \
        one-GCA per line. Only analyze GCAs from inputfolder that appear in \
        this list")
    parser.add_argument("-c", "--cpus", help="Number of cpus to pass to BUSCO\
        (default: all available)", type=int, default=cpu_count())
    parser.add_argument("-p", "--processes", help="Number of BUSCO processes to \
        launch simultaneously. Default: 2", default=2, type=int)
    return parser.parse_args()


def read_re_analyze(file_path):
    re_analyze_gca = set()
    
    if not file_path:
        return re_analyze_gca
    
    try:
        with open(file_path) as f:
            for line in f:
                if line[0] == "#" or line.strip() == "":
                    continue
                re_analyze_gca.add(line.strip())
    except IOError:
        pass
    
    return re_analyze_gca


def zipfile_ok(zipfile_path):
    if not zipfile_path.exists():
        return False
    
    if not is_zipfile(zipfile_path):
        return False
    
    try:
        ZipFile(zipfile_path)
    except BadZipFile:
        return False
    else:
        return True
    
    
def compress_folder(folder, name):
    """
    folder: path
    name: zip file
    """
    
    if not folder.is_dir():
        return
    
    try:
        with ZipFile(name, "w", compression=ZIP_DEFLATED, compresslevel=9) as z:
            for target in folder.glob("**/*"):
                #print(target)
                if target.is_dir():
                    continue
                
                # NOTE: the path needs to be relative to the end point of 'folder'
                # "relative_to" doesn't work...
                starting_folder = folder.parts[-1] # either 'augustus_output' or 'hmmer_output'
                target_starting_folder_idx = target.parent.parts.index(starting_folder)
                target_relative_path = "/".join(target.parent.parts[target_starting_folder_idx:])
                #print(target_relative_path)
                arcname_ = "{}/{}".format(target_relative_path, target.name)
                #print(arcname_)
                z.write(target, arcname=arcname_)
    except IOError:
        print("Error! compressing folder {} didn't work...".format(folder))
        return False
    else:
        rmtree(folder)
        return True
        

def busco(cpus, o, gca, db, zipfile, fna_filenames):
    # in parameters, use "delete=False" to inspect /tmp/*.fna files
    with tempfile.NamedTemporaryFile(prefix=gca, suffix=".fna") as fasta_file, \
            ZipFile(zipfile) as gcazip:
        # read the zipped genome and put it in the temporary file
        for zipped_fna in fna_filenames:
            fasta_file.write(gcazip.read(zipped_fna))
        
        cmd = []
        cmd.append("busco")
        cmd.extend(["--mode", "genome"])
        cmd.extend(["--cpu", str(cpus)])
        cmd.extend(["--out_path", str(o)])
        cmd.extend(["--out", gca])
        cmd.extend(["--lineage_dataset", str(db.name)])
        # TODO: choose something else here?
        cmd.extend(["--augustus_species", "saccharomyces_cerevisiae_S288C"])
        #cmd.append("--long") # I wonder how bad this can be
        cmd.extend(["--in", fasta_file.name])
        print(" ".join(cmd))
        
        try:
            proc = subprocess.run(cmd, stderr=STDOUT, encoding="utf-8")
        except FileNotFoundError as e:
            print("Error running busco command:")
            print(e)
            return False
        else:
            base_target_folder = o / gca / "run_ascomycota_odb10"
            
            augustus_folder = base_target_folder / "augustus_output/"
            hmmer_output_folder = base_target_folder / "hmmer_output/"
            busco_seq_folder = base_target_folder / "busco_sequences"
            
            augustus_zip = base_target_folder / "augustus_output.zip"
            hmmer_output_zip = base_target_folder / "hmmer_output.zip"
            busco_seq_zip = base_target_folder / "busco_sequences.zip"
            
            # Augustus
            # Check first as BUSCO 5 doesn't necessarily use Augustus
            if augustus_folder.is_dir():
                print("\tCompressing augustus output")
                compress_folder(augustus_folder, augustus_zip)
                # check if it worked
                if not zipfile_ok(augustus_zip):
                    print("Error zipping file {}".format(augustus_zip))
                
            # hmmer
            print("\tCompressing hmmer output")
            compress_folder(hmmer_output_folder, hmmer_output_zip)
            # check if it worked
            if not zipfile_ok(hmmer_output_zip):
                print("Error zipping file {}".format(hmmer_output_zip))
                
            # busco_sequences
            print("\tCompressing busco_seq output")
            compress_folder(busco_seq_folder, busco_seq_zip)
            if not zipfile_ok(busco_seq_zip):
                print("Error zipping file {}".format(busco_seq_zip))
                
            return True


if __name__ == "__main__":
    options = command_parser()
    
    # handle parameters
    cpus = options.cpus
    i = options.inputfolder
    if not i.is_dir():
        sys.exit("Error (--inputfolder). {} does not seem a valid folder".format(i))
    o = options.outputfolder
    if not o.is_dir():
        os.makedirs(o, exist_ok=True)
    db = options.dbfolder
    if not db.is_dir():
        sys.exit("Error (--dbfolder). {} does not seem a valid folder".format(db))
    re_analyze_gca = read_re_analyze(options.re_analyze_file)
    
    gca_filter = set()
    if options.filter_list:
        with open(options.filter_list) as f:
            for line in f:
                if line.strip() == "":
                    continue
                else:
                    gca_filter.add(line.strip)
    
    # traverse zip files
    with Pool(processes=options.processes) as pool:
        for zipfile in i.glob("*.zip"):
            gca = zipfile.stem
            
            # only calculate new stuff
            if re_analyze_gca:
                if gca not in re_analyze_gca:
                    continue
            
            if gca_filter:
                if gca not in gca_filter:
                    continue
            
            # Check if results folder exist already and we don't need to re-analyze
            if (o / gca).is_dir():
                # if no re-analyze file is given, this set is empty
                if gca not in re_analyze_gca:
                    continue
        
            try:
                with ZipFile(zipfile) as gcazip:
                    fna_filenames = set()
                    # traverse all names inside the zip. Genome may be split into
                    # more than one (chromosome-level) file
                    for item in gcazip.namelist():
                        if item.startswith("ncbi_dataset/data/"+gca) and item[-3:] == "fna":
                            fna_filenames.add(item)
                
                    if not fna_filenames:
                        print("Warning: could not find any .fna file for {}".format(gca))
                        continue
                    
                    # This doesn't work with apply_sync because of gcazip. 
                    # Someone on stackoverflow (questions/37907350) suggests that
                    # all parameters need to be pickle-able...
                    #busco(cpus, o, gca, db, gcazip, fna_filename, )
                    
                # so passing the zipfile location and opening on each children
                # process actually does the trick
                pool.apply_async(busco, args=(cpus, o, gca, db, zipfile, fna_filenames, ))
                    
            except BadZipFile:
                print("Warning: Cannot open {}".format(zipfile))
                continue
            
        pool.close()
        pool.join()
            
