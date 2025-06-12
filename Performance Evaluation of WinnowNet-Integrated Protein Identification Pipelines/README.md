# Performance Evaluation of WinnowNet-Integrated Protein Identification Pipelines

This folder contains the computational pipeline used for protein identification in our MS-based metaproteomics analysis, as described in our publication. The pipeline includes four tools (Sipros Ensemble, FragPipe, Peaks, and AlphaPept) and is organized to ensure reproducibility and compliance with the FAIR principles.

## Overview

The pipeline supports protein identification through database searching, filtering, and protein assembly using four major tools:

- [Sipros Ensemble](#sipros-ensemble)
- [FragPipe](#fragpipe)
- [Peaks](#peaks)
- [AlphaPept](#alphapept)
- [Sipros Ensemble (with WinnowNet)](#sipros-ensemble-with-winnownet)
- [FragPipe (with WinnowNet)](#fragpipe-with-winnownet)
- [Peaks (with WinnowNet)](#peaks-with-winnownet)
- [AlphaPept (with WinnowNet)](#alphapept-with-winnownet)

*Note*: All FDRs are controlled at 1% in these benchmarking experiments.
---

## FAIR Principles

- **Findable**: Each tool and script used in the pipeline is documented and referenced by name. Software tools are standard in the proteomics field.
- **Accessible**: We provide all command-line steps and screenshots required for GUI tools. If scripts are needed, users can find them in this folder or request them.
- **Interoperable**: Configuration files (e.g., `SiprosConfig.cfg`, `config.yaml`) follow community standards and can be reused.
- **Reusable**: This pipeline can be easily adapted to different metaproteomic datasets. Parameter settings are provided or referenced via screenshots.

---

### Software Download Links

- **Sipros Ensemble**  
  https://github.com/guo-xuan/Sipros-Ensemble/archive/refs/heads/master.zip

- **FragPipe (version 21.0)** and **Philosopher**
  https://github.com/Nesvilab/FragPipe/releases/download/21.0/FragPipe-21.0.zip

- **Peaks (version 12.5)**  
  https://www.bioinfor.com/peaks-studio-trial/ (Need to request a registration code for free trial)
  
- **AlphaPept (version 0.5)**  
  https://github.com/MannLabs/alphapept/archive/refs/tags/v0.5.0.zip or install the package using the command `pip install alphapept.`

---

## Original Pipelines

### Sipros Ensemble

#### Database Searching
```bash
Sipros_Openmp -f filename.ms2 -c work_Dir/SiprosConfig.cfg -o work_Dir
```

- `Sipros_Openmp` is in `Sipros-Ensemble-master/ReleaseOpenMP` folder.
- Raw MS data could be converted to `.ms2` by using MSConvert as in [Benchmark Dataset Descriptions](../Benchmark%20Dataset%20Descriptions/).
- Protein database could be updated in `SiprosConfig.cfg`. Reverse proteins could be added by executing this command `sipros_prepare_protein_database.py -i ./protein_database.fasta -o protein_database_rev.fasta -c ./SiprosConfig.cfg`

#### Filtering and Protein Assembly
```bash
./runSiprosFiltering.sh -in work_Dir -c work_Dir/SiprosConfig.cfg -o work_Dir
```

- *Note*: Adjust `FDR_Threshold` in [SiprosConfig.cfg](./SiprosConfig.cfg) to control filtering stringency.

---

### FragPipe

#### Database Searching
- Load parameters using GUI:
  - Screenshot: [FP_DB_search.png](./FP_DB_search.png)
  - Click "Load" to use [fragger.params](./fragger.params)

#### Filtering and Protein Assembly
- Use GUI settings shown in:
  - Screenshot: [FP_filtering_protein_assembly.png](./FP_filtering_protein_assembly.png)

---

### Peaks

Peaks provides a one-click workflow.

- GUI screenshot: [Peaks.png](./Peaks.png)

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

- Configuration file: [config.yaml](./config.yaml) (user must provide or customize)

---

## WinnowNet-Integrated Pipelines

The following alternative pipelines integrate WinnowNet for PSM filtering and scoring. They replace the original filtering steps with WinnowNet-based re-scoring and FDR control.

### Sipros Ensemble (with WinnowNet)
#### Database Searching
```bash
Sipros_Openmp -f filename.ms2 -c work_Dir/SiprosConfig.cfg -o work_Dir
```

#### Filtering with WinnowNet
```bash
# Step 1: Convert search results to TSV
python SE2win.py -i filename.Spe2Pep.txt -o filename.tsv

# Step 2: Re-score using WinnowNet
python SpectraFeatures.py -i filename.tsv -s filename.ms2 -o spectra.pkl -t 48 -f att
python Prediction.py -i spectra.pkl -o rescore.out.txt -m att_pytorch.pt

# Step 3: Apply FDR control at PSM/Peptide levels
python filtering_shuffle.py -i rescore.out.txt -p filename.tsv -o filtered -f 0.01
```

#### Protein Assembly
```bash
python sipros_peptides_assembling.py
```
*Note: Adjust `-f` in the previous step to ensure protein-level FDR is 1% and rerun if necessary.*

### FragPipe (with WinnowNet)
#### Database Searching
GUI-based setup (same as original): [FP_DB_search.png](./FP_DB_search.png)

#### Filtering with WinnowNet
```bash
# Step 1: Convert to TSV
python XML2win.py -i filename.pepXML -o filename.tsv

# Step 2: Re-score
python SpectraFeatures.py -i filename.tsv -s filename.ms2 -o spectra.pkl -t 48 -f att
python Prediction.py -i spectra.pkl -o rescore.out.txt -m att_pytorch.pt

# Step 3: Convert to ProteinProphet input
python win2prophet.py -i filename.pepXML -r rescore.out.txt -o filename.pep.xml
```

#### Protein Assembly
```bash
philosopher workspace --clean
philosopher workspace --init
philosopher proteinprophet --maxppmdiff 2000000 --minprob 0.5 --output combined filelist_proteinprophet.txt
philosopher.exe database --annotate protein_db.fasta --prefix shuffle_
philosopher filter --picked --prot 0.01 --minPepLen 7 --tag shuffle_ --pepxml combined.prot.xml --razor
philosopher report
```

### Peaks (with WinnowNet)
#### Database Searching
GUI-based setup: [Peaks.png](./Peaks.png)

#### Filtering with WinnowNet
```bash
# Step 1: Convert to TSV
python XML2win.py -i filename.pepXML -o filename.tsv

# Step 2: Re-score
python SpectraFeatures.py -i filename.tsv -s filename.ms2 -o spectra.pkl -t 48 -f att
python Prediction.py -i spectra.pkl -o rescore.out.txt -m att_pytorch.pt

# Step 3: Convert to ProteinProphet input
python win2prophet.py -i filename.pepXML -r rescore.out.txt -o filename.pep.xml
```

#### Protein Assembly
```bash
philosopher workspace --clean
philosopher workspace --init
philosopher proteinprophet --maxppmdiff 2000000 --minprob 0.5 --output combined filelist_proteinprophet.txt
philosopher.exe database --annotate protein_db.fasta --prefix shuffle_
philosopher filter --picked --prot 0.01 --minPepLen 7 --tag shuffle_ --pepxml combined.prot.xml --razor
philosopher report
```

### AlphaPept (with WinnowNet)
#### Database Searching
```bash
alphapept workflow config.yaml
```

#### Filtering with WinnowNet
```bash
# Step 1: Convert AlphaPept result
python alpha2win.py -i filename.csv -o filename.tsv

# Step 2: Re-score
python SpectraFeatures.py -i filename.tsv -s filename.ms2 -o spectra.pkl -t 48 -f att
python Prediction.py -i spectra.pkl -o rescore.out.txt -m att_pytorch.pt

# Step 3: Apply WinnowNet filtering
python alphapept_filtering.py -i filename.csv -r rescore.out.txt -o filename_output.csv
```

- `filename.csv` is the PSM result file ended by `_ids.csv`.

#### Protein Assembly
```bash
python AlphaPept_protein_assembly.py -i filename_output.csv -o filename_assembly_output.csv
```