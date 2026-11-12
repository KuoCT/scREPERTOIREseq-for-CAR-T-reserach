# Single-Cell Repertoire Sequencing Analysis for CAR-T Research

## Overview
This repository contains the R scripts used for single-cell gene expression (GEX) and repertoire (REP) sequencing analysis in CAR-T cell research. The workflow processes Cell Ranger outputs to perform downstream analyses and generate datasets for AI model training.

The repository is provided to support the reproducibility of the computational analyses presented in the associated study.

## Usage
The analysis workflow is controlled by `main.sh`, which executes the R scripts in the required order.

Run the following command from the repository root directory:

```bash
bash main.sh
```

## Requirements
The analysis requires R and the packages used by the GEX and repertoire analysis scripts.

The R session information files provide details of the environments used for the original analyses:

- `sessionInfo_RNA.txt`
- `sessionInfo_REP.txt`

## Citation
If you use this repository, please cite the associated publication and the archived software version.

- **Publication**: To be updated upon publication.
- **Software archive**: Zenodo DOI (pending release).
- **Raw sequencing data:** PRJNA1492306