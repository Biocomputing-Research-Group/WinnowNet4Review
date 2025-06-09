# WinnowNet: A Deep Learning-Based Filtering Framework for Metaproteomics

This repository contains the datasets, scripts, and instructions required to reproduce the experiments presented in our paper, particularly the sections on **Training WinnowNet**, **Benchmark Dataset Descriptions**, **Performance Comparison to State-of-the-art Filtering Algorithms**, and **Performance Evaluation of WinnowNet-Integrated Protein Identification Pipelines**.

WinnowNet is a deep-learning-based PSM filtering algorithm for metaproteomics that demonstrates superior performance over existing statistical and machine learning tools. This repository is organized to facilitate easy access and reproducibility of results.

## Repository Structure

```
WinnowNet-Reproduction/
│
├── Training WinnowNet/
│   └── Scripts and configuration files to train both the self-attention-based and CNN-based WinnowNet models.
│
├── Benchmark Dataset Descriptions/
│   └── Descriptions of the twelve benchmark datasets (Synthetic, P1–P3, Marine 1–3, Soil 1–3, Human Gut, and Human Gut TimsTOF).
│       Includes guidelines for data preprocessing and database construction with entrapment proteins.
│
├── Performance Comparison to State-of-the-art Filtering Algorithms/
│   └── Commands to benchmark WinnowNet against Percolator, Q-ranker, PeptideProphet, iProphet, MS²Rescore, and DeepFilter.
│       Includes results and plotting scripts for comparisons on PSM, peptide, and protein levels.
│
├── Performance Evaluation of WinnowNet-Integrated Protein Identification Pipelines/
│   └── Instructions for integrating WinnowNet with Sipros-Ensemble, FragPipe, Peaks, and AlphaPept pipelines.
│       Includes evaluation scripts and performance metrics (e.g., FDR control and identification gain).
│
└── README.md
```

## Reproducibility Overview

Each subdirectory contains:
- **Step-by-step instructions** for reproducing the results reported in the paper.
- **Dependencies and environment setup** details (e.g., specific software versions like ProteoWizard 3.0.11841).
- **Preprocessed or synthetic datasets**, where applicable.
- **Output files** and **scripts for evaluation**, including FDR calculations.

## Highlights

- **Entrapment-based FDR estimation** methods (including paired and combined strategies).
- **Self-contained training scripts** for WinnowNet, enabling retraining on new datasets.
- **Pipeline integration and benchmarking** across diverse platforms and sample complexities.

## Citation

If you find this repository helpful, please cite our paper:

> Citation coming soon.

## License

> MIT License.

## Contact

For any questions or issues, please [submit here](https://github.com/Biocomputing-Research-Group/WinnowNet4Review/issues)
