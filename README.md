# Constraints on Bacterial Genome Size  

HST.508 Final Project (Fall 2025)  

**Authors:** Aidan Pavao, Terry Cho  



## Project Description  

Bacterial genomes vary nearly an order of magnitude in size. Interestingly, bacteria from similar environments across diverse taxa tend to have similar genome sizes (e.g., human skin commensals ~2.5 Mb). It has been proposed that this convergence reflects metabolic complexity of the environment, while an alternative hypothesis suggests that genome size is constrained by nutrient limitation (e.g., nitrogen conservation).  



The goal of this project is to model how environmental factors constrain bacterial genome size. We will analyze genome size, metabolic coding capacity (via KEGG), and other genomic signatures of nutrient limitation (e.g., amino acid usage, pathway composition) across multiple well-defined environments.



## Data  

Filtered reference genomes were obtained from NCBI:  

[NCBI Reference Genome Dataset](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2&reference_only=true&annotated_only=true&refseq_annotation=true&genbank_annotation=true&typical_only=true&exclude_mags=true&exclude_multi_isolates=true&assembly_level=3:3)  

**Criteria:**
- Reference genomes only
- Complete assembly level
- Released since 2010
- Annotated (RefSeq and GenBank)
- Excluded atypical, MAGs, and multi-isolate assemblies

**Total genomes:** 5,986 (downloaded and organized)

**Data location:** `data/ncbi/`

**Directory structure:**
```
data/ncbi/
├── assemblies/          # 5,986 genome directories (126 GB)
│   └── GCF_*/          # One directory per assembly accession
│       ├── *_genomic.fna    # Genome sequences
│       ├── protein.faa      # Protein sequences
│       ├── cds_from_genomic.fna  # Coding sequences
│       ├── genomic.gff      # Gene annotations
│       └── genomic.gbff     # GenBank format
└── metadata/            # Summary files and manifests (29 MB)
    ├── assembly_manifest.txt              # List of all accessions
    ├── metadata_table.tsv                 # Basic metadata (6 columns)
    ├── assembly_data_report_extracted.tsv # Detailed metadata (73 columns)
    └── dataset_catalog_extracted.tsv      # File catalog
```

**For detailed documentation:** See `scripts/0_data_download/README.md` for complete directory structure, file descriptions, usage patterns, and script workflow.



## Preliminary Analysis  

Initial exploration includes:  

- Histogram and scatter plots of genome size distribution  

- Assessment of genome size variation across environments/ecology  

- Verification that data are appropriate for downstream KEGG and nutrient limitation analyses  

**Analysis scripts:** `scripts/1_exploratory_analyses/01_preliminary_analysis.ipynb`  



## Planned Analysis  

- KEGG annotation of protein sequences to estimate metabolic capacity  

- Correlation of genome size with environment, metabolic coding capacity, and nutrient limitation signatures  

- Model building (linear regression or machine learning) incorporating phylogeny and environment  



## Contributions  

- **Aidan Pavao:** Primary data analysis, visualization, modeling  

- **Terry Cho:** Data setup, repository organization, documentation  



## Repository Structure  

```
bac_genome_constraint/
├── data/                    # Input datasets
│   ├── ncbi/               # NCBI genome data (5,986 assemblies)
│   └── metadata/           # Metadata files
├── scripts/                # Analysis and plotting scripts
│   └── 0_data_download/    # NCBI data download scripts and documentation
├── results/                # Output figures and summary tables
└── README.md               # Project overview
```

**Data download documentation:** `scripts/0_data_download/README.md`



## Timeline  

- **Data criteria due:** Nov 7, 2025  

- **Preliminary analysis due:** Nov 17, 2025  

- **Presentation:** Dec 8–10, 2025  

- **Final paper due:** Dec 12, 2025  

