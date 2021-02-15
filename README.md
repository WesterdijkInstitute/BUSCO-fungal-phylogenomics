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
conda activate busco412
```

We will also need a specific [BUSCO database](https://busco.ezlab.org/frames/fungi.htm). For this project, I used the [Ascomycota Odb10 set](https://busco-data.ezlab.org/v4/data/lineages/ascomycota_odb10.2019-11-20.tar.gz). Download and decompress.


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


# Launch BUSCO on each assembly

Next step will be to launch BUSCO on the assemblies contained on each zipped file. If a results folder is found (`[output folder]/[accession]`), the BUSCO analysis will be skipped for the corresponding assembly. 

:warning: A BUSCO results folder could have been created but the run may have actually failed (e.g. user cancelled, lack of space, etc.). So be careful with this simply check when restarting the analysis!

Internally, the script tries to read the assembly zip file and traverses its internal structure to find fasta files with the assembly's sequences. With the list of files, it launches a process that will join all the sequence files into one temporary file and use it as input for BUSCO. When done, it will compress the contents of three output folders within `[outputfolder]/[accession]/run_ascomycota_odb10`: `augustus_output`, `hmmer_output` and `busco_sequences`. The reason for this is that each of these subfolders contain thousands of small output files, which can create fragmentation issues for hard drives.

:warning: Originally, I created a separate script to compress all results (when I started to have fragmentation issues). When I integrated some of this code into `2_launch_busco` to compress the results after the BUSCO run, it turned out that the zipping part was non-funcional due to differences in the Python versions used. I created a new environment with BUSCO 4.1.2, which includes a newer version of Python where zipping works. The BUSCO version I used originally was **4.0.6** but as far as I can see, no significant changes were made to the algorithm (just [bugfixes](https://gitlab.com/ezlab/busco/-/blob/master/CHANGELOG)).

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



# Check busco results

* Script: `check_busco_results.py`. Performs a quick check comparing the numbers reported in each BUSCO results folder 'short_summary' file against the actual number of files with sequences.
* Input:
  - Folder with all BUSCO results
  - Tab-separated `metadata file`
* Output:
  - `busco_set_results_summary.tsv`. Example:

```
Assembly	[C]omplete BUSCOs	Complete and [S]ingle-copy BUSCOs	Complete and [D]uplicated BUSCOs	[F]ragmented BUSCOs	[M]issing BUSCOs	Name
GCA_000002515.1	1610	1607	3	6	90	Kluyveromyces lactis NRRL Y-1140
GCA_000002525.1	1583	1575	8	5	118	Yarrowia lipolytica CLIB122 CLIB122
GCA_000002545.2	1601	1563	38	4	101	[Candida] glabrata CBS 138
```

