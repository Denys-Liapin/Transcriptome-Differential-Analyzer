## RNA-Seq Differential Expression & Enrichment Analyzer

A R pipeline for differential expression analysis (DEA) using DESeq2 and functional enrichment analysis (GSEA & ORA) via clusterProfiler. Optimized for 2x2 factorial designs.

## How to Run

### Installation of necessary libraries
```bash
sudo Rscript install.R
```

### Start of analysis
```bash
Rscript Analyzer.R data.xlsx
```

# Pipeline Features & Logic
Duplicate Gene Resolution: Automatically handles duplicate gene identifiers using make.unique() to ensure mathematical compatibility with DESeq2 engine requirements.

Dynamic Downstream Protection: Safely skips visualization steps if a specific contrast drops below the absolute statistical threshold of <2 significant genes, preventing runtime script crashes.

Interaction Evaluation: Runs complete statistical interaction designs (Genotype:Diet) to assess whether independent variables compound synergistically or act sequentially as standard downstream physiological phenotypic triggers.

<details> <summary>

MA plot example

</summary>

<img alt="image" src="https://github.com/Denys-Liapin/Transcriptome-Differential-Analyzer/blob/main/output%20example/MA_Disease_100.png" />

</details>

<details> <summary>

Heatmap plot example

</summary>

<img alt="image" src="https://github.com/Denys-Liapin/Transcriptome-Differential-Analyzer/blob/main/output%20example/Heatmap_Disease_50.png" />

</details>
