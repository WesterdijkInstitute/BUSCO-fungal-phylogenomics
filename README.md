# Saccharomycotina taxonomy

# Preamble: conda environment

We'll use a conda environment with the BUSCO software and its dependencies. First, [install conda](https://docs.conda.io/en/latest/miniconda.html). 

The current guide to install BUSCO using conda are [here](https://busco.ezlab.org/busco_userguide.html#conda-package) but in this repository there is a file to install version 4 (which is what I used). To install this version, do:

```
conda env create --file busco4env.yml
```

# Obtain data

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

