# miRNA_CAD_Thesis

**Author:** Viktorija Daukšaitė (@Viktorija0719)  
**Thesis Title:** *Elucidating the Expression Profile of Circulating Micro-RNAs in Early-Onset and Late-Onset Coronary Artery Disease Patients*

## Overview

This repository contains the complete codebase, datasets, and documentation for my master's thesis, which investigates the role of circulating microRNAs (miRNAs) as potential biomarkers in the diagnosis and prognosis of coronary artery disease (CAD). The study uses computational methods to analyze miRNA expression profiles and identify significant patterns and associations related to early- and late-onset CAD.

## Getting Started

### Prerequisites

- Python 3.8 or higher
- R 4.0 or higher

### Installation

1. Clone the repository:

```
git clone https://github.com/Viktorija0719/miRNA_CAD_Thesis.git
cd miRNA_CAD_Thesis
```

2. (Optional) Create and activate a virtual environment:

```
python -m venv venv
source venv/bin/activate        # On Windows: venv\Scripts\activate
```

3. Install Python dependencies:

```
pip install -r requirements.txt
```

4. Install R dependencies by opening an R session and running:

```
install.packages("BiocManager")
BiocManager::install(c("tidyverse", "limma", "edgeR", "caret", "EnhancedVolcano",
                       "ggplot2", "patchwork", "cowplot", "glmnet", "pROC", 
                       "dplyr", "stringr", "tibble", "tidyr", "optparse", "psych"))
install.packages("randomForest")
```

## Usage

- Use the scripts in the `scripts/` directory for data preprocessing and model training.
- To run the full pipeline:

```
./scripts/run_all.sh
```

## Synthetic Data Generation

As part of this project, synthetic datasets were generated to supplement the original miRNA expression data. The following files contain synthetic data:

- `synthetic_clinical_data.csv`
- `countData.txt`

These datasets were created using the [CTGANSynthesizer](https://docs.sdv.dev/sdv/single-table-data/modeling/synthesizers/ctgansynthesizer) from the [Synthetic Data Vault (SDV)](https://docs.sdv.dev/sdv) library.

**Reference:**
Xu, L., Skoularidou, M., Cuesta-Infante, A., & Veeramachaneni, K. (2019). Modeling Tabular data using Conditional GAN. *Advances in Neural Information Processing Systems*, 32. [Link](https://papers.nips.cc/paper/2019/hash/254ed7d2de3b23ab10936528e9c54a9f-Abstract.html)

The synthetic data were generated to enhance the dataset for analysis while ensuring data privacy and compliance with ethical standards.



