# Metaproteomics Protein Identification Pipeline

This folder contains the computational pipeline used for protein identification in our MS-based metaproteomics analysis, as described in our publication. The pipeline includes multiple tools (Sipros Ensemble, FragPipe, Peaks, and AlphaPept) and is organized to ensure reproducibility and compliance with the FAIR principles.

## Overview

The pipeline supports protein identification through database searching, filtering, and protein assembly using four major tools:

- [Sipros Ensemble](#sipros-ensemble)
- [FragPipe](#fragpipe)
- [Peaks](#peaks)
- [AlphaPept](#alphapept)

---

## FAIR Principles

- **Findable**: Each tool and script used in the pipeline is documented and referenced by name. Software tools are standard in the proteomics field.
- **Accessible**: We provide all command-line steps and screenshots required for GUI tools. If scripts are needed, users can find them in this folder or request them.
- **Interoperable**: Configuration files (e.g., `SiprosConfig.cfg`, `config.yaml`) follow community standards and can be reused.
- **Reusable**: This pipeline can be easily adapted to different metaproteomic datasets. Parameter settings are provided or referenced via screenshots.

---

## Tool-Specific Pipelines

### Sipros Ensemble

#### Database Searching
```bash
Sipros_Openmp -f filename.ms2 -c work_Dir/SiprosConfig.cfg -o work_Dir
```

#### Filtering and Protein Assembly
```bash
./runSiprosFiltering.sh -in work_Dir -c work_Dir/SiprosConfig.cfg -o work_Dir
```

- *Note*: Adjust `FDR_Threshold` in `SiprosConfig.cfg` to control filtering stringency.

---

### FragPipe

#### Database Searching
- Load parameters using GUI:
  - Screenshot: `FP_DB_search.png`
  - Click "Load" to use `fragger.params`

#### Filtering and Protein Assembly
- Use GUI settings shown in:
  - Screenshot: `FP_filtering_protein_assembly.png`

---

### Peaks

Peaks provides a one-click workflow.

- GUI screenshot: `Peaks.png`

**Parameter adjustments:**
- Precursor mass tolerance: 0.09 Da
- Fragment ion tolerance: 0.01 Da
- Enzyme: Trypsin
- PTM: Oxidation
- Peptide length: 7 to 60

---

### AlphaPept

One-click command for the full workflow:
```bash
alphapept workflow config.yaml
```

- Configuration file: `config.yaml` (user must provide or customize)

---

## Additional Notes

- Screenshots referenced in this README must be placed in the same folder.
- Please include any used scripts that are not publicly available or request them if needed.

