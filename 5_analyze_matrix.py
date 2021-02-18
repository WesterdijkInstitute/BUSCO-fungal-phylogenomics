#! /usr/bin/env python

"""
Reads back the matrix created in previous script and filters assemblies and BUSCOs
"""

import sys
import os
import argparse
from pathlib import Path
import pandas as pd


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--matrix", help="Path to tsv busco a/p matrix",
        required=True, type=Path)
    parser.add_argument("-l", "--links", help="Path to 'links_to_ODB10.txt' \
        file, which contains information about the BUSCO genes. It will be \
        used for the gene report. Optional", type=Path)
    parser.add_argument("-s", "--summary", help="Path to the summary .tsv file, \
        which contains information about the assembly set. It will be used \
        for the assembly reports. Optional", type=Path)    
    parser.add_argument("-t", "--threshold", help="A number between 0 and 1 \
        representing the percentage of Busco completeness (relative to the \
        number of total Busco hits in the set) necessary to pass to downstream \
        analysis. Assemblies below this number will also be reported. Default: \
        0.7", type=float, default=0.7)
    return parser.parse_args()


if __name__ == "__main__":
    args = arg_parser()
    
    i = args.matrix
    if not i.is_file():
        sys.exit("Error: given input folder not a file")
        
    t = args.threshold
    if t < 0.0 or t > 1.0:
        sys.exit("Error: --threshold argument must be in the range [0.0, 1.0]")
    
    gene_info = dict()
    if args.links:
        with open(args.links) as f:
            for line in f:
                x = line.strip().split("\t")
                gene_info[x[0]] = "{}\t{}".format(x[1], x[2])
                
    asm_info = dict()
    # expects the following header
    report_header = "Assembly\t[C]omplete BUSCOs\tComplete and [S]ingle-copy BUSCOs\tComplete and [D]uplicated BUSCOs\t[F]ragmented BUSCOs\t[M]issing BUSCOs\tName\n"
    if args.summary:
        with open(args.summary) as f:
            for line in f:
                x = line.strip().split("\t")
                asm_info[x[0]] = "\t".join(x[1:])
    
    # read Busco Hits dataframe
    bh = pd.read_csv(i, sep="\t", header=0, index_col=0)
    print("Got a dataframe of {} assemblies and {} busco genes.".format(len(bh.index), len(bh.columns)))
    
    # make a new dataframe: completeness
    bh.loc[:, "Completeness"] = bh.sum(axis=1)
    print("\nSample of the absence/presence matrix, including completeness:")
    print(bh.head())
    print()
    
    
    # Filter data frame based on requested completeness threshold
    # i.e. find "good" assemblies
    minimum_hits = int(len(bh.columns) * t) # calculate minimum busco hits based on threshold
    x = bh["Completeness"].ge(minimum_hits) # alternative way of selecting last column
    filtered_bh = bh[x]
    print("{}/{} assemblies pass the requirement of having at least {} ({}) BUSCO hits".format(len(filtered_bh.index), \
        len(bh.index), minimum_hits, t))
    
    # Report underperforming assemblies. Can be used with the launch_busco script to try to re-analyze
    discarded_bh = bh[~x]
    with open("matrix_analysis_Bottom_{:04.2f}_Assemblies.tsv".format(1-t), "w") as f:
        f.write(report_header)
        for asm in discarded_bh.index:
            f.write("{}\t{}\n".format(asm, asm_info.get(asm, "\t\t\t\t\t")))
        #f.write("\n".join(discarded_bh.index))
    
    
    # proceed with filtered data frame. Which genes are present in all assemblies?
    busco_sharedness = filtered_bh.iloc[:,:-1].all(axis=0) # creates a True/False series. Leave out "Completeness" column
    shared_by_all = filtered_bh.iloc[:, :-1].columns[busco_sharedness] # data series of busco ids shared by all assemblies
    # instead of the following, report genes for various presence (enrichment?) levels
    with open("matrix_analysis_Top_{:04.2f}_Assemblies.tsv".format(t), "w") as f:
        for asm in filtered_bh.index:
            f.write(report_header)
            f.write("{}\t{}\n".format(asm, asm_info.get(asm, "\t\t\t\t\t")))
    print("{}/{} busco genes are present in all remaining assemblies".format(len(shared_by_all), \
        len(bh.columns)-1)) # "-1" because we have extra "Completeness" extra column
    
    # Output genes for various levels of presence in assemblies
    # i.e. analyze columns to find "good" genes
    print("\nBUSCO hits analysis (on {} filtered assemblies)".format(len(filtered_bh)))
    print("Target asm. enrichment %\tTarget asm. enrichment #\tBUSCOs found in target num. of asms.")
    for asm_perc in [1.0, 0.95, 0.9]:
        # basically, filter all gene labels with 'asms' or higher presence in assemblies
        asms = int(asm_perc * len(filtered_bh.index)) # min. num of assemblies the gene must be found (S)
        
        # make the filter
        hits_presence_bool = filtered_bh.iloc[:,:-1].sum(axis=0).ge(asms, axis=0)
        
        # from filtered assembly df, select everything except for the 
        # completeness column, then filter using hits_presence_bool
        genes = filtered_bh.columns[:-1][hits_presence_bool]
        print("{}\t{}\t{}".format(asm_perc, asms, len(genes)))
        
        with open("matrix_analysis_S_genes_in_{:04.2f}_of_Top_{:04.2f}_assemblies.tsv".format(asm_perc, t), "w") as f:
            for g in genes:
                f.write("{}\t{}\n".format(g, gene_info.get(g, "\t")))
        
