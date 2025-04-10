# ckIGeS-BS

This repository contains IGeS-BS workflow and relevant analysis.

---

## Repository Structure

```bash
├── IGeS_BS_tool/             # IGeS-BS workflow
	├──  IGeS_Profiling   
	├──  IGeS_BS
	├──  Data	  
	├──  dependency
├── analysis/                 # Scripts used for data preprocessing and constructing IGeS-BS and result analysis
└── README.md                 # Repository documentation
```

---

## Requirements

The scripts and tools in this repository were developed and tested using the following software environment:

* R ≥ 4.0.2 ()the development of IGeS-BS and relevant analysis version)
* openjdk version "1.8.0_372" for running GSEA step

A full list of R library dependencies for environment setup can be found in `environment.yml`.

## IGeS Usage

### 1. Clone the repository

clone IGeS-BS local folder:

```bash
git clone https://github.com/LiHongCSBLab/IGeS_BS.git
cd IGeS_BS_tool
```

### 2. Set up the environment

Using Conda (recommended):

```bash
conda env create -f environment.yml
conda activate IGeS_BS_env
```

For checking minimum required dependency, please take a look at: `dependency_check.R`

### 3. Run IGeS-profiling

#### Data dependency

Current profling support two raw data input from LINCS:IGeS-BS

**GSE70138:**

gene expression profile - GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx

    *Please download from link below and place under ./IGeS_BS_tool/Data/LINCS_data*

meta info: GSE70138_Broad_LINCS_inst_info.txt (provided in `./Data` folder)

**GSE92742:**

gene expression profile - GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx

    *Please download from link below and place under ./IGeS_BS_tool/Data/LINCS_data*

meta info: GSE92742_Broad_LINCS_inst_info.txt (provided in `./Data` folder)

Additional annotation files:

gene annotation file: GSE92742_Broad_LINCS_gene_info.txt (provided in `./Data` folder)

### 4. Run IGeS-BS

After profiling for each Compounds, boosting score can be calculated with command-line script.

**command-line options:**

```bash
filepath=/Path/to/IGeS_Profiling_output/
workpath=/Path/to/toolFolder/
savepath=/Path/to/IGeS_BS_output/

Rscript --vanilla plot_for_drug_IGeS.R -a 70138 -c LIHC -t Primary -u TUMERIC -f FALSE \
-d "$filepath}" \
-w "${workpath}" \
-s "${savepath}"> IGeS-BS.log 2>&1&

```

## Data Availability

Due to size constraints, the datasets used in this study are hosted externally. They can be accessed from the following links:

| Dataset Description             | Link                                                                                   |
| ------------------------------- | -------------------------------------------------------------------------------------- |
| GSE70138 data                   | [Download from GEO-GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138) |
| GSE92742 dataset                | [Download from GEO-GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) |
| in-house mouse experiment data | [Download from NODE [National Omics Data Encyclopedia](https://www.biosino.org/node/project/detail/OEP0000602) |

> **Note:** Please cite the corresponding data repository when using these datasets in your work.

---

## Citation

If you use this code, please cite the following:

> [Full citation of the manuscript]

BibTeX:

```bibtex
@article{YourCitationKey,
  title     = {Title of Manuscript},
  author    = {Author Names},
  journal   = {Journal Name},
  year      = {Year},
  doi       = {DOI}
}
```

---

## License

This project is licensed under the [MIT License / GPL-3.0 / Other]. See the `LICENSE` file for details.

---

## Contact

For questions or feedback, please contact:

Li Hong

Email: lihong01@sinh.ac.cn

Affiliation: Shanghai Institute of Nutrition and Health
