# Saccharomycotina taxonomy

# Preamble

We'll use a conda environment with the BUSCO software and its dependencies. First, [install conda](https://docs.conda.io/en/latest/miniconda.html). 

The current guide to install BUSCO using conda are [here](https://busco.ezlab.org/busco_userguide.html#conda-package) but they are for the new version 5.0.

Using conda:
```
conda create -n busco406 -c bioconda -c conda-forge busco=4.0.6 biopython=1.77 python>=3.7 pandas
```
This also installs [Pandas](https://pandas.pydata.org/), which will be used later on.

Another (easier) option is to use a `.yml` file included in this repository by doing 
```
conda env create --file busco406.yml
```

We'll use that environment to run the rest of the scripts. Activate with
```
conda activate busco406
```

We will also need a specific [BUSCO database](https://busco.ezlab.org/frames/fungi.htm). For this project, I used the [Ascomycota Odb10 set](https://busco-data.ezlab.org/v4/data/lineages/ascomycota_odb10.2019-11-20.tar.gz). Download and decompress.

For the latter part of the project, another environment will be used, containing `MAFFT`, `trimal` and `IQ-Tree` (these could in theory be included in the `busco406` environment but it's already difficult for conda to solve all those dependencies):
```
conda create -n phylogeny -c bioconda mafft trimal iqtree python=3
```

# Obtain information about available assemblies

The assemblies were obtained using NCBI's Datasets tool. [Download it](https://www.ncbi.nlm.nih.gov/datasets/docs/command-line-start/) and make sure it's available in your path. I have renamed the executable from `datasets` to `ncbi-datasets`. I am using version `10.21.0`.

Obtain a `json file` with the information of all assemblies. Taxon `147537` corresponds to [*Saccharomycotina*](https://www.ncbi.nlm.nih.gov/taxonomy/147537). The command after the pipe character will re-format the file to make it human-readable:
```
ncbi-datasets summary genome taxon 147537 | python -m json.tool > saccharomycotina_2021-02-08.json
```

The `json file` now looks like this:
```JSON
{
    "assemblies": [
        {
            "assembly": {
                "annotation_metadata": {},
                "assembly_accession": "GCA_001600815.1",
                "assembly_category": "representative genome",
                "assembly_level": "Scaffold",
                "bioproject_lineages": [
                    {
                        "bioprojects": [
                            {
                                "accession": "PRJDB3734",
                                "parent_accessions": [
                                    "PRJDB3556"
                                ],
                                "title": "NBRP: Genome sequencing of Alloascoidea hylecoeti JCM 7604"
                            },
                            {
                                "accession": "PRJDB3556",
                                "title": "NBRP: Genome sequencing of eukaryotic microorganisms"
                            }
                        ]
                    }
                ],
                "chromosomes": [
                    "Un"
                ],
                "contig_n50": 224558,
                "display_name": "JCM_7604_assembly_v001",
                "estimated_size": "7801514",
                "org": {
                    "assembly_counts": {
                        "node": 1,
                        "subtree": 1
                    },
                    "key": "54196",
                    "parent_tax_id": "1540230",
                    "rank": "SPECIES",
                    "sci_name": "Alloascoidea hylecoeti",
                    "strain": "JCM 7604",
                    "tax_id": "54196",
                    "title": "Alloascoidea hylecoeti"
                },
                "seq_length": "24815695",
                "submission_date": "2016-03-01"
            }
        },
...
```


# Download assemblies

Now we'll download each assembly listed in that file.

* Script: `1_get_assemblies_from_json.py`
* Input: the `json file`
* Output: 
  - a set of zip files with the assemblies
  - a `metadata.tsv` file.
* Usage:
```
usage: 1_get_assemblies_from_json.py [-h] -j JSON [-t TRIES] [-o OUTPUTFOLDER]
                                     [-n N]

optional arguments:
  -h, --help            show this help message and exit
  -j JSON, --json JSON  JSON file downloaded with NCBI datasets with assembly
                        data
  -t TRIES, --tries TRIES
                        Try this many times to get each accession file. If it
                        still fails, a warning will be issued. Default: 3
                        times
  -o OUTPUTFOLDER, --outputfolder OUTPUTFOLDER
                        Folder where the assemblies will be stored (default:
                        assemblies/)
  -n N                  Optional: Number of assemblies to download (in the
                        order of the JSON file)
```

If pointing to a previous output folder, the script will verify whether each file already exists (and can be opened). This allows easy updating of the assembly files.

The `metadata.tsv` file contains formatted information from the `json file`: assembly accession, NCBI tax ID, species name and strain:
```
GCA_001600815.1	54196	Alloascoidea hylecoeti	JCM 7604
GCA_001600695.1	1301101	Ascoidea asiatica	JCM 7603
GCF_001661345.1	1344418	Ascoidea rubescens DSM 1968	DSM 1968
...
```

The `metadata` file for the latest results is [here](./files/metadata_2021-01.tsv).


# Launch BUSCO on each assembly

Next step will be to launch BUSCO on the assemblies contained on each zipped file. If a results folder is found (`[output folder]/[accession]`), the BUSCO analysis will be skipped for the corresponding assembly. 

:warning: A BUSCO results folder could have been created but the run may have actually failed (e.g. user cancelled, lack of space, etc.). So be careful with this simply check when restarting the analysis!

Internally, the script tries to read the assembly zip file and traverses its internal structure to find fasta files with the assembly's sequences. With the list of files, it launches a process that will join all the sequence files into one temporary file and use it as input for BUSCO. When done, it will compress the contents of three output folders within `[outputfolder]/[accession]/run_ascomycota_odb10`: `augustus_output`, `hmmer_output` and `busco_sequences`. The reason for this is that each of these subfolders contain thousands of small output files, which can create fragmentation issues for hard drives.

:warning: Originally, I created a separate script to compress all results (when I started to have fragmentation issues). When I integrated some of this code into `2_launch_busco` to compress the results after the BUSCO run, it turned out that the zipping part was non-funcional due to differences in the Python versions used. I created a new environment keeping the same BUSCO version, but a newer Python version. I have tested this new environment on a random assembly and it works.

* Script: `2_launch_busco.py`
* Input: folder with zipped assemblies, path to folder with BUSCO database
* Optional input: 
  - A file with assembly accessions to be re-analyzed in case they have a results folder already
  - A file with assembly accessions to filter input. If used, only the accessions on the input folder that match these accessions will be analyzed.
* Output: 
  - A folder with BUSCO results
* Parameters for BUSCO command: `--mode genome --lineage_dataset [path to ascomycota_odb10] --augustus_species saccharomyces_cerevisiae_S288C`
* Usage:
```
usage: 2_launch_busco.py [-h] -i INPUTFOLDER [-o OUTPUTFOLDER] -d DBFOLDER
                         [--re_analyze_file RE_ANALYZE_FILE]
                         [--filter_list FILTER_LIST] [-c CPUS] [-p PROCESSES]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFOLDER, --inputfolder INPUTFOLDER
                        Folder with zipped assemblies
  -o OUTPUTFOLDER, --outputfolder OUTPUTFOLDER
                        Folder with BUSCO results (Default: 'Busco_results')
  -d DBFOLDER, --dbfolder DBFOLDER
                        Folder with a BUSCO database
  --re_analyze_file RE_ANALYZE_FILE
                        A file containing 1-assembly accession per line. Those
                        assemblies will be re-analized
  --filter_list FILTER_LIST
                        Read a txt file with one-GCA per line. Only analyze
                        GCAs from inputfolder that appear in this list
  -c CPUS, --cpus CPUS  Number of cpus to pass to BUSCO (default: all
                        available)
  -p PROCESSES, --processes PROCESSES
                        Number of BUSCO processes to launch simultaneously.
                        Default: 2
```

Each BUSCO result folder will contain a subfolder with data specific to the database used (in this case, `ascomycota_odb10`). Inside this foler, a small file contains a summary of the results (`[outputfolder]/[accession]/run_ascomycota_odb10/short_summary.txt`). For example, for assembly `GCA_001600695.1`, the `short_summary` file includes de following:
```
	C:81.7%[S:79.2%,D:2.5%],F:0.6%,M:17.7%,n:1706
	1395	Complete BUSCOs (C)
	1352	Complete and single-copy BUSCOs (S)
	43	Complete and duplicated BUSCOs (D)
	10	Fragmented BUSCOs (F)
	301	Missing BUSCOs (M)
	1706	Total BUSCO groups searched
```

:warning::warning::warning: When writing this guide, I tried versions 4.1.4 (when dealing with the compression issue) and 5 of BUSCO. BUSCO v5.0.0 uses another default method to find genes ([MetaEuk](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00808-x) instead of Augustus), so I expected differences. What I didn't expect were differences with v4.1.4. With the same test assembly as above, these are the `Missing BUSCOs (M)`:

```
BUSCO 4.1.4:
	751	Missing BUSCOs (M)
	
BUSCO 5.0.0:
	243	Missing BUSCOs (M)
```
As can be seen, v4.1.4 misses a lot of BUSCO hits while v5.0.0 looks better. I will continue to work with the v4.0.6 results but will try to make a full v5.0.0 run in the new server. BUSCO's [changelog](https://gitlab.com/ezlab/busco/-/blob/master/CHANGELOG) doesn't mention any breaking changes (only bug fixes), so it's difficult to say why this version seems to perform much worse than v4.0.6.

Three assemblies had issues and could not be processed
- GCA_009666835.1
- GCA_015345625.1
- GCA_015345745.1


# Verify BUSCO results

The next script reads the `short_summary` files and compares the number of singe-copy BUSCOs (S) reported there with the actual number of files. It also produces a report of all summaries, [`busco_set_results_summary`](./files/busco_set_results_summary.tsv).

* Script: `3_verify_busco_results.py`
* Input:
  - Folder with all BUSCO results
  - Tab-separated `metadata file`
* Output:
  - `busco_set_results_summary.tsv`.
* Usage
```
usage: 3_verify_busco_results.py [-h] -b BUSCOFOLDERS [-m METADATA]

optional arguments:
  -h, --help            show this help message and exit
  -b BUSCOFOLDERS, --buscofolders BUSCOFOLDERS
                        Path to folder with busco results (each result is a
                        subfolder)
  -m METADATA, --metadata METADATA
                        A tab-separated file in which first column is the
                        assembly name, third column is the species name, and
                        fourth column is the strain name. Will be used to
                        annotate summary file (Optional).
```

Example of the `busco_set_results_summary` file:
```
Assembly	[C]omplete BUSCOs	Complete and [S]ingle-copy BUSCOs	Complete and [D]uplicated BUSCOs	[F]ragmented BUSCOs	[M]issing BUSCOs	Name
GCA_000002515.1	1610	1607	3	6	90	Kluyveromyces lactis NRRL Y-1140
GCA_000002525.1	1583	1575	8	5	118	Yarrowia lipolytica CLIB122 CLIB122
GCA_000002545.2	1601	1563	38	4	101	[Candida] glabrata CBS 138
```


# Make an absence/presence matrix of all BUSCO results

With the BUSCO results per assembly, we want to evaluate which of these genes are present in all (or most) assemblies. For this we will build a presence/absence matrix of all BUSCO results (only complete and single copy hits).

* Script: `4_make_busco_a-p_matrix.py`
* Input: a path to the folder with all BUSCO results
* Output: an `a/p matrix` (rows: assemblies; columns: BUSCOs) in tsv format
* Usage:
```
usage: 4_make_busco_a-p_matrix.py [-h] -i INPUTFOLDER [-f FILTER_LIST]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFOLDER, --inputfolder INPUTFOLDER
                        Path to folder with busco results (each result is a
                        subfolder)
  -f FILTER_LIST, --filter_list FILTER_LIST
                        Optional. File with list of assemblies. Only the
                        assemblies from the input folder that are in this list
                        will be processed (it can be a tab-separated file)

```


# Analyze the absence/presence matrix

As a simple strategy to analyze the whole BUSCO results, we'll analyze the absence/presence matrix using the "completeness" of each assembly.

In the first stage, the assemblies will be analyzed. Only the assemblies that contain at least a desired number of BUSCO hits in the `Complete and single-copy BUSCOs (S)` category will be used downstream. This completeness threshold is `0.7` by default (1,194/1,706 BUSCO hits for the `ascomycota_odb10` database).

A second stage analyzes the BUSCOs in that filtered set of assemblies. Currently, the script reports lists of genes that are found in 100%, 95% and 90% of (filtered) assemblies.


* Script: `5_analyze_matrix.py`
* Input: 
  - The `a/p matrix`
  - Optional: the `links_to_ODB10.txt` file, part of the `ascomycota_odb10` contents.
  - Optional: the `busco_set_results_summary.tsv` file from script 3.
  - Threshold of completeness to filter assemblies.
* Output:
  - The list of assemblies that have the requested number of BUSCOs ("Top Assemblies"). They will be annotated if the `busco_set_results_summary.tsv` file was used.
  - The list of assemblies that did not contain the requested number of (S) BUSCOs ("Bottom Assemblies")
  - Lists of genes for different levels of presence in the filtered assemblies. They will be annotated if the `links_to_ODB10.txt` file was used.
* Usage:
```
usage: 5_analyze_matrix.py [-h] -m MATRIX [-l LINKS] [-s SUMMARY]
                           [-t THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  -m MATRIX, --matrix MATRIX
                        Path to tsv busco a/p matrix
  -l LINKS, --links LINKS
                        Path to 'links_to_ODB10.txt' file, which contains
                        information about the BUSCO genes. It will be used for
                        the gene report. Optional
  -s SUMMARY, --summary SUMMARY
                        Path to the summary .tsv file, which contains
                        information about the assembly set. It will be used
                        for the assembly reports. Optional
  -t THRESHOLD, --threshold THRESHOLD
                        A number between 0 and 1 representing the percentage
                        of Busco completeness (relative to the number of total
                        Busco hits in the set) necessary to pass to downstream
                        analysis. Assemblies below this number will also be
                        reported. Default: 0.7
```

For example, for the default top `0.7` assemblies, 14 BUSCOs were found in all those assemblies.


# Obtain the sequence of all BUSCOs

With the list of `Top Assemblies` and enriched BUSCOs, the following script will read all the (zipped) BUSCO results to extract the (DNA or AA) sequences into unaligned fasta files for each BUSCO.

If a sequence is not found in any results folder, the script currently inserts an empty sequence in the file (i.e. only its header)

* Script: `6_assemble_unaligned_TargetGenes.py.py`
* Input: 
  - The base BUSCO results folder
  - A list of assemblies. For example, `matrix_analysis_Top_0.70_Assemblies.tsv`
  - A list of BUSCOs. For example, `matrix_analysis_S_genes_in_1.00_of_Top_0.70_assemblies.tsv`
* Output:
  - Fasta files for each BUSCO
* Usage:
```
usage: 6_assemble_unaligned_TargetGenes.py [-h] -r RESULTS -a ASSEMBLIES -t
                                           TARGETGENES [-o OUTPUTFOLDER]
                                           [--aa]

optional arguments:
  -h, --help            show this help message and exit
  -r RESULTS, --results RESULTS
                        Base path with BUSCO results
  -a ASSEMBLIES, --assemblies ASSEMBLIES
                        List of assemblies from which to extract the sequences
  -t TARGETGENES, --targetgenes TARGETGENES
                        List of Target Genes
  -o OUTPUTFOLDER, --outputfolder OUTPUTFOLDER
                        Folder with unaligned sequence file. Default:
                        ./output/Target_Genes_unaligned
  --aa                  Extract protein sequences instead of DNA
```


# Multiple sequence alignment and curation

With the set of sequence files, the following script will launch `MAFFT` on each file. `Trimal` will be used to trim each alignment file (using option `-gappyout`). The `phylogeny` environment (see first section) will be used here:

* Script: `7_align_Target_Genes.py`
* Input: a folder with fasta files (can be nucloetide or amino acid)
* Output: 
  - A folder with the aligned versions of the fasta files in the input
  - A folder with the trimmed versions of each aligned file
* Usage:
```
usage: 7_align_Target_Genes.py [-h] -i INPUTFOLDER -a ALIGNEDFOLDER -t TRIMMEDFOLDER [-p PROCESSES] [--threads THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFOLDER, --inputfolder INPUTFOLDER
                        Folder with unaligned .fasta files
  -a ALIGNEDFOLDER, --alignedfolder ALIGNEDFOLDER
                        Folder for aligned sequence files.
  -t TRIMMEDFOLDER, --trimmedfolder TRIMMEDFOLDER
                        Folder for trimmed aligned sequence files (using command 'trimal -gappyout')
  -p PROCESSES, --processes PROCESSES
                        Number of MAFFT instances. Default: 1
  --threads THREADS     --threads parameter for each MAFFT process. Default: 2
```


# Concatenate alignments

The last script before launching IQ-Tree reads the set of aligned and curated files from the previous steps and concatenates their sequences, producing also the required partition file.

* Script: `8_concatenate_alignments.py`
* Input: a folder with all curated alignments
* Output:
  - A concatenated sequence file. 
  - A nexus partition file
```
usage: 8_concatenate_alignments.py [-h] -i INPUTFOLDER -n NAME

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFOLDER, --inputfolder INPUTFOLDER
                        Folder with aligned sequences
  -n NAME, --name NAME  Base name for the output
```


# Use IQ-Tree

Use the previous files for the phylogenomic analysis. There are many options for IQ-Tree. E.g.:
```
iqtree -s [concatenated].fasta -p [partition].nex -T AUTO -B 1000 --msub nuclear -m MFP
```
Instead of the ultrafase bootstrap, we can substitute `-B 1000` with `--fast`. The `-m MFP` will use the "new ModelFinder" which will try to find the most appropriate substitution model for each partition, out of the models designed for sequences localized in the nuclear genome (option `--msub nuclear`)
