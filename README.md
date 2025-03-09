# NOTE
This pipeline analyzes RNA crosslinking and proximity ligation data using DuplexDiscovereR and the [Snakemake workflow management system](https://snakemake.readthedocs.io/en/stable/) 🐍. 

🌟 Reference 🌟<br/>
https://bioconductor.org/packages/release/bioc/manuals/DuplexDiscovereR/man/DuplexDiscovereR.pdf

## Installing Snakemake
Please follow the [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install the Snakemake workflow management tool. We recommend using `Conda/Mamba` to install Snakemake.

This Snakemake workflow has been tested with `v7.32.4` and Python `v3.11.11`.

## Python Virtual Environment Requirement
A virtual environment is required for the last step of this pipeline when merging duplex groups obtained from various samples. Once installed, please add the path to `config/config.yaml`.

## Downloading FASTQ files from GEO
This Snakemake workflow includes the following datasets. The datasets/samples to be included in the workflow can be modified in `config/config.yaml`.

**PARIS**<br/>
[Lu et al. 2016 Cell](https://www.cell.com/fulltext/S0092867416304226#secsectitle0085)
| Accession  | Sample Name  |
| -----------| ------------ |
| SRR2814761 | HeLa         |
| SRR2814762 | HeLa         |
| SRR2814763 | HEK293T      |
| SRR2814764 | HEK293T      |
| SRR2814765 | HEK293T      |

**PARIS2**<br/>
[Zhang et al. 2021 Nature Communications](https://www.nature.com/articles/s41467-021-22552-y)<br/>
[Zhang et al. 2023 PNAS](https://pubmed.ncbi.nlm.nih.gov/37792516/)
| Accession   | Sample Name                          |
| ----------- | ------------------------------------ |
| SRR11624581 | HEK293T - AMT                        |
| SRR11624582 | HEK293T - AMT                        |
| SRR11624583 | HEK293T - AMT                        |
| SRR11624584 | HEK293T - AMT                        |
| SRR11624585 | HEK293T - Amoto                      |
| SRR11624586 | HEK293T - Amoto                      |
| SRR11624587 | HEK293T - Amoto                      |
| SRR11624588 | HEK293T - Amoto                      |
| SRR11624589 | HEK293T - mRNA                       |
| SRR11951629 | HEK293T - snoRNA                     |
| SRR24883901 | Neuron - chromatin associated RNA    |
| SRR24883902 | Astrocyte - chromatin associated RNA |
| SRR24883903 | NPC - chromatin associated RNA       |
| SRR24883904 | iPS - chromatin associated RNA       |
| SRR24883905 | ECs - chromatin associated RNA       |
| SRR24883906 | HEK293T - chromatin associated RNA   |
| SRR24883907 | SHSY5Y - snoRNA                      |
| SRR24883908 | HEK293T - snoRNA                     |
| SRR24883909 | HEK293T - Total RNA                  |

**SPLASH**<br/>
[Aw et al. 2016 Molecular Cell](https://www.sciencedirect.com/science/article/pii/S1097276516301046#:~:text=SPLASH%20Uncovers%20New%20rRNA%2DrRNA,intra%2D%20and%20intermolecular%20RNA%20interactions.)
| Accession  | Sample Name                |
| ---------- | -------------------------- |
| SRR3404924 | Lymphoblastoid - Total RNA |
| SRR3404925 | Lymphoblastoid - Total RNA |
| SRR3404926 | hES - PolyA                |
| SRR3404927 | RA - PolyA                 |
| SRR3404928 | RA - PolyA                 |
| SRR3404936 | Lymphoblastoid - Total RNA |
| SRR3404937 | Lymphoblastoid - Total RNA |
| SRR3404938 | Lymphoblastoid - snoRNA    |
| SRR3404939 | Lymphoblastoid - PolyA     |
| SRR3404940 | Lymphoblastoid - PolyA     |
| SRR3404941 | Lymphoblastoid - PolyA     |
| SRR3404942 | Lymphoblastoid - PolyA     |
| SRR3404943 | hES - PolyA                |

**LIGR-seq**<br/>
[Sharma et al. 2016 Molecular Cell](https://www.sciencedirect.com/science/article/pii/S109727651630106X?via%3Dihub)
| Accession  | Sample Name |
| ---------- | ----------- |
| SRR3361013 | HEK293T     |
| SRR3361017 | HEK293T     |

The FASTQ files can be downloaded using:
- `download_wget.smk`: Both the `SRRXXXXXXX` run number and the download URL (Available on [SRA Explorer](https://sra-explorer.info/#)) are required as input. Add the run number to `config/config.yaml` and the URL to `resources/sra_ids.txt`.

We recommend verifying the download of all required datasets before moving on to the next steps of the pipeline.

## Running the Snakemake workflow
For a dry-run of this Snakemake workflow, simply run the following code from `DuplexDiscovereR-Snakemake/`.
```
snakemake -n
```
To run this Snakemake workflow, simply run the following code from `DuplexDiscovereR-Snakemake/`.
```
snakemake --profile profile_slurm

OR

snakemake --profile profile_local
```

If you have any specific questions regarding this Snakemake workflow, please contact [Kristina Sungeun Song](mailto:kristina.song@usherbrooke.ca). Questions on the technicalities of DuplexDiscovereR should be addressed to the original authors.
