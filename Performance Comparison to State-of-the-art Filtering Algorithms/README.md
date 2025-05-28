# Benchmarking Metaproteomics Search Engines

This directory contains the commands and configuration files used to run protein database searches with Comet, Myrimatch, and MS-GF+ for multiple metaproteomics datasets in our study.

## Software Download URLs

- [Comet (version 2018012)](https://sourceforge.net/projects/comet-ms/files/comet_2018012.zip)  
- [Myrimatch](https://figshare.com/articles/dataset/myrimatch/29144549?file=54806417)  
- [MS-GF+ (version 2019.04.18)](https://github.com/MSGFPlus/msgfplus)  

---

## Running Comet

Place `comet.params` in the same directory as `comet.exe`. Update the `database_name` parameter in the file for each dataset as shown below.

**Common command:**
```
export OMP_NUM_THREADS=16
./comet.exe workDir/*.ms2
```

**Configuration file:** `comet.params`

```
database_name = ./Marine_shuffled.fasta # for Marine Metaproteomes
database_name = ./Soil_shuffled.fasta # for Soil Metaproteomes
database_name = ./Mock_Comm_RefDB_V3_shuffled.fasta # for Mock Community Metaproteomes
database_name = ./stool_nrnew_shuffled.fasta # for Human Gut Metaproteome
```

---

## Running Myrimatch

Place `myrimatch.cfg` in your working directory. Update the `ProteinDatabase` field for each dataset.

**Common command:**
```
export LC_ALL=C
./myrimatch -cfg myrimatch.cfg -workdir workDir/ -cpus 16 workDir/*.ms2
```

**Configuration file:** `myrimatch.cfg`
```
ProteinDatabase = ./Marine_shuffled.fasta # Marine
ProteinDatabase = ./Soil_shuffled.fasta # Soil
ProteinDatabase = ./Mock_Comm_RefDB_V3_shuffled.fasta # Mock Community
ProteinDatabase = ./stool_nrnew_shuffled.fasta # Human Gut
```

---

## Running MS-GF+

Place `MSGFPlus_Mods1.txt` and the appropriate FASTA database in the same `workDir`.

**Common command:**
```
cd workDir
for name in *.ms2; do
java -Xmx250000M -jar MSGFPlus.jar
-s "${name}"
-d "./DATABASE.fasta"
-inst 1 -t 0.09Da -ti -1,3 -ntt 2 -e 1
-thread 16 -m 3 -minLength 7 -maxLength 60 -n 5
-addFeatures 1 -maxMissedCleavages 3 -numMods 1
-decoy Rev -mod MSGFPlus_Mods1.txt
-o "${name%.}.mzid"
done
```

**Configuration file:** `MSGFPlus_Mods1.txt`  
Replace `DATABASE.fasta` with the appropriate file:  
- Marine: `Marine_shuffled.fasta`  
- Soil: `Soil_shuffled.fasta`  
- Mock Community: `Mock_Comm_RefDB_V3_shuffled.fasta`  
- Human Gut: `stool_nrnew_shuffled.fasta`  

---

## Notes

- All configuration files (`comet.params`, `myrimatch.cfg`, `MSGFPlus_Mods1.txt`) and FASTA databases should be placed in the `experiments/` directory.  
- Ensure Java 8 or newer is installed for MS-GF+.  
- All commands were executed using 16 CPU threads.  