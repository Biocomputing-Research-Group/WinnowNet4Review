# Benchmarking Metaproteomics State-of-the-Art Filtering Algorithms

This directory provides all necessary commands, configuration files, and conversion scripts to reproduce our benchmarking of **WinnowNet** against six leading post-processing algorithms:

- **Percolator**
- **Q-ranker**
- **PeptideProphet**
- **iProphet**
- **MS2Rescore**
- **DeepFilter**

We used peptide-spectrum matches (PSMs) generated from three protein database search engines: **Comet**, **Myrimatch**, and **MS-GF+**, applied across multiple metaproteomics datasets. More details about these datasets are in [Data/](../Data/)

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

**Command:**
```
export LC_ALL=C
./myrimatch -cfg myrimatch.cfg -workdir workDir/ -cpus 16 workDir/*.ms2
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
    -decoy Rev -mod MSGFPlus_Mods1.txt \
    -o "${name%.*}.mzid"
done
```

**Replace `DATABASE.fasta` with:**

- `Marine_shuffled.fasta`
- `Soil_shuffled.fasta`
- `Mock_Comm_RefDB_V3_shuffled.fasta`
- `stool_nrnew_shuffled.fasta`

---

## 2. PSM Filtering and Protein Inference with FDR Control

### 2.1. Percolator

Percolator is a machine learning-based tool for post-processing PSMs. We benchmarked version 3.2.

**Download URL:**  
https://github.com/percolator/percolator

---

#### Input Conversion to `.pin` Format

##### Comet

- **Script Needed:** None  
- **Note:** Comet directly outputs `.pin` files during execution.

##### Myrimatch

- **Script:** `myrimatch2pin.py`
- **Command:**
  ```
  python myrimatch2pin.py -i filename.pepXML -o filename.pin
  ```
  - `-i`: Myrimatch output (`.pepXML`)  
  - `-o`: Output `.pin` for Percolator

##### MS-GF+

- **Script:** `msgf2pin.py`
- **Command:**
  ```
  python msgf2pin.py -i filename.mzid -o filename.pin
  ```
  - `-i`: MS-GF+ output (`.mzid`)  
  - `-o`: Output `.pin` for Percolator

---

#### Running Percolator

Run Percolator on a `.pin` file:

```
percolator filename.pin -m filename.target.tsv -M filename.decoy.tsv
```

- `-m`: Output file for target PSMs  
- `-M`: Output file for decoy PSMs

---

#### FDR Filtering (PSM/Peptide Level)

Apply 1% FDR control at the PSM and peptide level:

```
python filtering_benchmark_shuffle.py -t filename.target.tsv -d filename.decoy.tsv -m percolator -o filename.filtered -f 0.01
```

- `-t`: Target TSV file from Percolator  
- `-d`: Decoy TSV file from Percolator  
- `-m`: Method (set to `percolator`)  
- `-o`: Output prefix  
- `-f`: FDR threshold (e.g., `0.01` for 1%)

**Output files:**

- `filename.filtered.psm.txt`  
- `filename.filtered.pep.txt`

---

#### Protein-Level Assembly and FDR Control

Assemble peptides into proteins and assess FDR:

```
python sipros_peptides_assembling.py
```

> **Note:** To control FDR at the protein level to 1%, iteratively reduce the `-f` parameter in the filtering step above and re-run until the desired protein-level FDR is achieved.

---

### 2.2. MS2Rescore

[MS2Rescore documentation and installation](https://ms2rescore.readthedocs.io/en/stable/installation/)

MS2Rescore is a rescoring tool designed to improve peptide-spectrum match classification using additional features. Below are the instructions for converting search engine outputs, running MS2Rescore, and controlling FDR at multiple levels.

---

#### Input Conversion to MS2Rescore-Compatible `.tsv` Format

##### Comet

- **Script:** `comet2tsv.py`
- **Command:**
  ```
  python comet2tsv.py -i filename.pep.xml -o filename.tsv
  ```

##### Myrimatch

- **Script:** `myrimatch2tsv.py`
- **Command:**
  ```
  python myrimatch2tsv.py -i filename.pep.xml -o filename.tsv
  ```

##### MS-GF+

- **Script:** `msgf2tsv.py`
- **Command:**
  ```
  python msgf2tsv.py -i filename.pep.xml -o filename.tsv
  ```

---

#### Running MS2Rescore

Once `.tsv` files are prepared and you have a `config.json` file ready:

```
ms2rescore -c config.json -p filename.tsv -n 20
```

- `-c`: Configuration JSON file
- `-p`: Input `.tsv` file generated by conversion scripts
- `-n`: Number of CPU cores to use

---

#### FDR Filtering (PSM/Peptide Level at 1%)

You can apply 1% FDR filtering using:

```
python filtering_benchmark_shuffle.py -t output_directory/q-ranker.target.psms.txt -d filename.decoy.tsv -m ms2rescore -o filename.filtered -f 0.01
```

**Alternative command (e.g., when using MokaPot outputs):**

```
python filtering_benchmark_shuffle.py -t filename.ms2rescore.mokapot.psms.txt -d filename.ms2rescore.mokapot.decoy.psms.txt -m ms2rescore -o filename.filtered -f 0.01
```

This generates:
- `filename.filtered.psm.txt` (PSM-level results)
- `filename.filtered.pep.txt` (Peptide-level results)

---

#### Protein-Level Assembly and FDR Control

Use the following command to assemble peptides and assess protein-level FDR:

```
python sipros_peptides_assembling.py
```

> **Note:** To control protein-level FDR to 1%, iteratively reduce the `-f` value in the filtering command and re-run until the resulting protein FDR reaches 1%.


## 3. General Notes

- All configuration files (`comet.params`, `myrimatch.cfg`, `MSGFPlus_Mods1.txt`) and FASTA databases should be placed in the `experiments/` directory.
- Use absolute paths where necessary for reproducibility.
- All commands were executed using 16 CPU threads unless otherwise specified.
- Java (version 8 or later) is required for MS-GF+.
- Make sure the `.fasta` files used are the **shuffled decoy versions** referenced in the paper.
