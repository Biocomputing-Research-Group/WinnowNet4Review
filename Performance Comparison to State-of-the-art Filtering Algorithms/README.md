# Benchmarking Metaproteomics State-of-the-Art Filtering Algorithms

This directory provides all necessary commands, configuration files, and conversion scripts to reproduce our benchmarking of **WinnowNet** against six leading post-processing algorithms:

- **Percolator**
- **Q-ranker**
- **PeptideProphet**
- **iProphet**
- **MS2Rescore**
- **DeepFilter**

We used peptide-spectrum matches (PSMs) generated from three protein database search engines: **Comet**, **Myrimatch**, and **MS-GF+**, applied across multiple metaproteomics datasets. More details about these datasets are in [Data/](./Data/). All FDRs are controlled at 1% in these benchmarking experiments.

---

## Contents

- [1. Protein Database Searching](#1-protein-database-searching)
- [2. PSM Filtering and Protein Inference with FDR Control](#2-psm-filtering-and-protein-inference-with-fdr-control)
  - [Percolator](#21-percolator)
  - [Q-ranker](#22-q-ranker)
  - [PeptideProphet](#23-peptideprophet)
  - [iProphet](#24-iprophet)
  - [MS2Rescore](#25-ms2rescore)
  - [DeepFilter](#26-deepfilter)
  - [CNN-based WinnowNet](#27-cnn-based-winnownet)
  - [Self-attention-based WinnowNet](#28-self-attention-based-winnownet)

---

## 1. Protein Database Searching

### Software Download Links

- **Comet (version 2018012)**  
  https://sourceforge.net/projects/comet-ms/files/comet_2018012.zip

- **Myrimatch (version 2.2, released in 2012)**  
  https://figshare.com/articles/dataset/myrimatch/29144549?file=54806417

- **MS-GF+ (version 2019.04.18)**  
  https://github.com/MSGFPlus/msgfplus

---

### Running Comet

1. Place `comet.params` in the same directory as `comet.exe`.  
2. Update the `database_name` field in `comet.params` to point to the appropriate FASTA file.
3. Comet will include decoy search if using the provided configuration file.

**Command:**
```
export OMP_NUM_THREADS=16
./comet.exe workDir/*.ms2
```

**Relevant snippet in `comet.params`:**
```
database_name = ./Marine_shuffled.fasta           # Marine
database_name = ./Soil_shuffled.fasta             # Soil
database_name = ./Mock_Comm_RefDB_V3_shuffled.fasta  # Mock Community
database_name = ./stool_nrnew_shuffled.fasta         # Human Gut
```

---

### Running Myrimatch

1. Place `myrimatch.cfg` in your working directory.  
2. Modify the `ProteinDatabase` field as needed.
3. Myrimatch will include decoy search if using the provided configuration file.

**Command:**
```
export LC_ALL=C
./myrimatch -cfg /path/to/myrimatch.cfg -workdir /path/to/workDir/ -cpus 16 /path/to/workDir/*.ms2
```

**Relevant snippet in `myrimatch.cfg`:**
```
ProteinDatabase = ./Marine_shuffled.fasta           # Marine
ProteinDatabase = ./Soil_shuffled.fasta             # Soil
ProteinDatabase = ./Mock_Comm_RefDB_V3_shuffled.fasta  # Mock Community
ProteinDatabase = ./stool_nrnew_shuffled.fasta         # Human Gut
```

---

### Running MS-GF+

1. Place `MSGFPlus_Mods1.txt` and the corresponding `.fasta` database in the `workDir`.  
2. Ensure Java 8 or newer is installed.
3. MS-GF+ will include decoy search if using the provided parameter option `-tda 1`.

**Command:**
```
cd workDir
for name in *.ms2; do
  java -Xmx250000M -jar /path/to/MSGFPlus.jar \
    -s "$name" \
    -d "./DATABASE.fasta" \
    -inst 1 -t 0.09Da -ti -1,3 -ntt 2 -e 1 \
    -thread 16 -m 3 -minLength 7 -maxLength 60 -n 5 \
    -addFeatures 1 -maxMissedCleavages 3 -numMods 1 \
    -decoy Rev -tda 1 -mod MSGFPlus_Mods1.txt \
    -o "${name%.*}.mzid"
done
```

**Replace `DATABASE.fasta` with:**
```
- Marine_shuffled.fasta
- Soil_shuffled.fasta
- Mock_Comm_RefDB_V3_shuffled.fasta
- stool_nrnew_shuffled.fasta
```

---

## 2. PSM Filtering and Protein Inference with FDR Control

In this section, we detail how to post-process PSMs and infer proteins using different filtering algorithms. Each subsection explains:

- **Download link** for the tool
- **Input conversion** from search engine output to the tool-specific format
- **Command** to run the tool
- **FDR filtering** at PSM/peptide level (1% FDR in these benchmark experiments, `sipros_post_module.py` and `parseconfig.py` are needed.)
- **Protein inference** step to assemble peptides into proteins and control protein-level FDR

---

### 2.1. Percolator

**Download URL:**  
https://github.com/percolator/percolator (Version 3.2)

#### Input Conversion to `.pin` Format

- **Comet**: Produces `.pin` files natively; no conversion script needed.  
- **Myrimatch**:  
  ```bash
  python myrimatch2pin.py -i filename.pepXML -o filename.pin
  ```  
  (`-i`: input Myrimatch `.pepXML`; `-o`: output `.pin`)  
- **MS-GF+**:  
  ```bash
  python msgf2pin.py -i filename.mzid -o filename.pin
  ```  
  (`-i`: input MS-GF+ `.mzid`; `-o`: output `.pin`)

#### Running Percolator

To combine search engine output files for Percolator, concatenate all .pin files while retaining only the first file's header line. Use the following command:

```bash
head -n 1 file1.pin > filename.pin && tail -n +2 -q *.pin >> filename.pin
```

Then run Percolator:

```bash
percolator filename.pin -m filename.target.tsv -M filename.decoy.tsv
```
- `-m`: output target PSMs TSV  
- `-M`: output decoy PSMs TSV

#### FDR Filtering (PSM/Peptide Level)

Apply 1% FDR:
```bash
python filtering_benchmark_shuffle.py -t filename.target.tsv -d filename.decoy.tsv -m percolator -o filename.filtered -f 0.01
```
- `-t`: Percolator target TSV  
- `-d`: Percolator decoy TSV  
- `-m`: method name `percolator`  
- `-o`: output prefix  
- `-f`: FDR threshold (e.g., `0.01`)

**Output files:**  
- `filename.filtered.psm.txt`  
- `filename.filtered.pep.txt`

#### Protein-Level Assembly and FDR Control

```bash
python sipros_peptides_assembling.py
```
> Note: To achieve 1% protein-level FDR, adjust the `-f` parameter in the previous filtering step and rerun until the desired protein FDR is met.

---

### 2.2. Q-ranker

**Download URL:**  
https://figshare.com/articles/dataset/crux/29206184?file=55020527

#### Input Conversion to Q-ranker-Compatible Format

- **Comet**:  
  ```bash
  python comet2txt.py -i filename.pep.xml -o filename.txt
  ```  
- **Myrimatch**:  
  ```bash
  python myriamtch2txt.py -i filename.pepXML -o filename.txt
  ```  
- **MS-GF+**:  
  ```bash
  python msgf2txt.py -i filename.mzid -o filename.txt
  ```

#### Running Q-ranker

To combine search engine output files for Q-ranker, concatenate all .txt files while retaining only the first file's header line. Use the following command:

```bash
head -n 1 file1.txt > filename.txt && tail -n +2 -q *.txt >> filename.txt
```

Then run Q-ranker:

```bash
./crux q-ranker --decoy-prefix Rev_ --output-dir output_directory filename.ms2 filename.txt
```

#### FDR Filtering (PSM/Peptide Level)

```bash
python filtering_benchmark_shuffle.py -t output_directory/q-ranker.target.psms.txt -d output_directory/q-ranker.decoy.psms.txt -m qranker -o filename.filtered -f 0.01
```
- Generates: `filename.filtered.psm.txt`, `filename.filtered.pep.txt`

#### Protein-Level Assembly

```bash
python sipros_peptides_assembling.py
```
> Note: Adjust `-f` threshold iteratively to meet 1% protein FDR.

---

### 2.3. PeptideProphet

**Download URL:**  
http://tools.proteomecenter.org/wiki/index.php?title=TPP_5.2.0:_Installing_on_Ubuntu_18.04_LTS

#### Input Format

- **Comet** & **Myrimatch**: Generate `.pep.xml` files directly.  
- **MS-GF+**: Convert using TPP’s `idconvert`:
  ```bash
  idconvert filename.mzid --pepXML
  ```

#### Running PeptideProphet

```bash
InteractParser combined_filename.pep.xml *.pep.xml -Dprotein.fasta -Tfasta -Estricttrypsin -a /ms2_file_directory/
PeptideProphetParser combined_filename.pep.xml ZERO DECOY=Rev_ DECOYPROBS
```

> This produces a combined `.pep.xml` with probability scores for each PSM.

#### FDR Filtering (PSM/Peptide Level)

```bash
python filtering_benchmarkXML_shuffle.py -i combined.pep.xml -m peptideprophet -o filename.filtered -f 0.01
```

Output:
- `filename.filtered.psm.txt`: PSMs filtered at 1% FDR
- `filename.filtered.pep.txt`: Peptides filtered at 1% FDR

#### Protein-Level Assembly

```bash
python sipros_peptides_assembling.py
```
This script assembles peptides into proteins and reports the false discovery portion.

> To control protein-level FDR at 1%, adjust the `-f` value in the filtering step to a lower threshold and rerun this script until the desired FDR is achieved.

---

### 2.4. iProphet

**Download URL:**  
http://tools.proteomecenter.org/wiki/index.php?title=TPP_5.2.0:_Installing_on_Ubuntu_18.04_LTS

#### Input

Uses the `.pep.xml` output file generated by PeptideProphet.

#### Running iProphet

```bash
InterProphetParser THREADS=8 DECOY=Rev_ NONSI NONSM NONSP NONRS combined_filename.pep.xml iProphet_output.pep.xml
```
> iProphet refines probability estimates by integrating multiple search engine results.

#### FDR Filtering (PSM/Peptide Level)

```bash
python filtering_benchmarkXML_shuffle.py -i iProphet_output.pep.xml -m iprophet -o filename.filtered -f 0.01
```

Generates:
- `filename.filtered.psm.txt`
- `filename.filtered.pep.txt`

#### Protein-Level Assembly

```bash
python sipros_peptides_assembling.py
```

Adjust the `-f` value and re-run to achieve 1% protein-level FDR.

---

### 2.5. MS2Rescore

**Documentation & Installation:**  
https://ms2rescore.readthedocs.io/en/stable/installation/

MS2Rescore enhances PSM confidence using additional features like retention time and MS2 spectral intensity.

#### Input Conversion to `.tsv` Format

- **Comet**:  
  ```bash
  python comet2tsv.py -i filename.pep.xml -o filename.tsv
  ```  
- **Myrimatch**:  
  ```bash
  python myrimatch2tsv.py -i filename.pep.xml -o filename.tsv
  ```  
- **MS-GF+**:  
  ```bash
  python msgf2tsv.py -i filename.mzid -o filename.tsv
  ```

#### Running MS2Rescore

To combine search engine output files for MS2Rescore, concatenate all .tsv files while retaining only the first file's header line. Use the following command:

```bash
head -n 1 file1.tsv > filename.tsv && tail -n +2 -q *.tsv >> filename.tsv
```

Then run Q-ranker:

```bash
ms2rescore -c config.json -s filename.mgf -p filename.tsv -n 20
```
- `-c`: path to `config.json` which specifies model parameters
- `-s`: path to the mass spectrometry data in mgf format, the conversion configuration for MSConvertGUI is shown in [here](./tomgf.png)
- `-p`: input PSM `.tsv`  
- `-n`: number of CPU cores

#### FDR Filtering (PSM/Peptide Level)

Use output from MS2Rescore (e.g., Mokapot PSMs):

```bash
python filtering_benchmark_shuffle.py -t filename.ms2rescore.mokapot.psms.txt -d filename.ms2rescore.mokapot.decoy.psms.txt -m ms2rescore -o filename.filtered -f 0.01
```
- Produces: `filename.filtered.psm.txt`, `filename.filtered.pep.txt`

#### Protein-Level Assembly and FDR Control

```bash
python sipros_peptides_assembling.py
```
> Adjust `-f` to reach 1% protein FDR.

---

### 2.6. DeepFilter

**Documentation & Installation:**  
https://github.com/Biocomputing-Research-Group/DeepFilter

DeepFilter applies a deep learning model to rescore PSMs.

#### Input Conversion to `.pin` Format

- **Comet**: Native `.pin` output; no conversion needed.  
- **Myrimatch**:  
  ```bash
  python myrimatch2df.py -i filename.pepXML -o filename.pin
  ```  
- **MS-GF+**:  
  ```bash
  python msgf2df.py -i filename.mzid -o filename.pin
  ```

#### Running DeepFilter

```bash
./inference.sh -in filename.ms2 -s filename.pin -m benchmark.pt -o rescore.out.txt
```

#### FDR Filtering (PSM/Peptide Level)

```bash
python filtering.py rescore.out.txt filename.pin filename.psm.txt filename.pep.txt
```
- `filename.psm.txt`: filtered PSM results  
- `filename.pep.txt`: filtered peptide results

#### Protein-Level Assembly and FDR Control

```bash
python sipros_peptides_assembling.py
```
> Iterate `-f` to ensure 1% protein FDR.

---

### 2.7. CNN-based WinnowNet

**Download URL:**  
https://github.com/Biocomputing-Research-Group/WinnowNet/tree/main

A CNN-based deep learning approach for PSM rescoring.

#### Input Conversion to `.tsv` Format

- **Comet**:  
  ```bash
  python comet2win.py -i filename.pep.xml -o filename.tsv
  ```  
- **Myrimatch**:  
  ```bash
  python myrimatch2win.py -i filename.pepXML -o filename.tsv
  ```  
- **MS-GF+**:  
  ```bash
  python msgf2win.py -i filename.mzid -o filename.tsv
  ```

#### Running CNN-based WinnowNet

1. Generate spectral features:  
   ```bash
   python SpectraFeatures.py -i filename.tsv -s filename.ms2 -o spectra.pkl -t 48 -f cnn
   ```
2. Predict PSM scores:  
   ```bash
   python Prediction_CNN.py -i spectra.pkl -o rescore.out.txt -m cnn_pytorch.pt
   ```

#### FDR Filtering (PSM/Peptide Level)

```bash
python filtering_shuffle.py -i rescore.out.txt -p filename.tsv -o filtered
```

#### Protein-Level Assembly and FDR Control

```bash
python sipros_peptides_assembling.py
```
> Adjust `-f` for 1% protein FDR.

---

### 2.8. Self-attention-based WinnowNet

**Download URL:**  
https://github.com/Biocomputing-Research-Group/WinnowNet/tree/main

A transformer-inspired self-attention method for PSM rescoring.

#### Input Conversion to `.tsv` Format

- **Comet**:  
  ```bash
  python comet2win.py -i filename.pep.xml -o filename.tsv
  ```  
- **Myrimatch**:  
  ```bash
  python myrimatch2win.py -i filename.pepXML -o filename.tsv
  ```  
- **MS-GF+**:  
  ```bash
  python msgf2win.py -i filename.mzid -o filename.tsv
  ```

#### Running Self-attention-based WinnowNet

1. Generate spectral features:  
   ```bash
   python SpectraFeatures.py -i filename.tsv -s filename.ms2 -o spectra.pkl -t 48 -f att
   ```
2. Predict PSM scores:  
   ```bash
   python Prediction.py -i spectra.pkl -o rescore.out.txt -m marine_att.pt
   ```

#### FDR Filtering (PSM/Peptide Level)

```bash
python filtering_shuffle.py -i rescore.out.txt -p filename.tsv -o filtered
```

#### Protein-Level Assembly and FDR Control

```bash
python sipros_peptides_assembling.py
```
> Adjust `-f` until protein FDR is ≤1%.

---

## 3. General Notes

- All configuration files (`comet.params`, `myrimatch.cfg`, `MSGFPlus_Mods1.txt`, etc.) and FASTA databases should be placed in the `experiments/` directory.  
- Use absolute paths for reproducibility.  
- Most commands were run using 16 CPU threads.  
- Java (version 8 or newer) is required for MS-GF+.  
- Ensure FASTA databases are the **shuffled decoy versions** listed in the paper.

---
