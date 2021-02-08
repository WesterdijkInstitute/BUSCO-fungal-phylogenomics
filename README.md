# Saccharomycotina taxonomy

# Obtain data

The assemblies were obtained using NCBI's Datasets tool. [Download it](https://www.ncbi.nlm.nih.gov/datasets/docs/command-line-start/) and make sure it's available in your path. I have renamed the executable from `datasets` to `ncbi-datasets`.

Obtain a `json file` with the information of all assemblies. Taxon `147537` corresponds to *Saccharomycotina*. The command after the pipe character will re-format the file to make it human-readable:
```
ncbi-datasets summary genome taxon 147537 | python -m json.tool > saccharomycotina_2021-02-08.json
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

