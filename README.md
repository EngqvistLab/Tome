# Tome: Temperature optima for microorgianisms and enzymes
A machine learning model for the prediction of optimal growth temperature of microorganisms<br/>
This method has two submodules: (1) model for the prediction of optimal growth temperature; (2) find the homologues for a given seqeuence with the same EC nubmer. 

## Installation
Pyhton 2.7 or Python 3 should both work.
## Depedences
* ncbi-blast-2.7.1+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* pandas
* Biopython
* numpy
* collections
* sklearn


## Usage:
### 1. Prediction of optimal growth temperature for one microorganism
```linux
./tome predOGT -fasta test/proteomes/95_pyrococcus_horikoshii_archaea.fasta
```
Then you will get following results:<br/>
```linux
FileName	predOGT (C)
95_pyrococcus_horikoshii_archaea.fasta	94.0
```

### 2. Predict optimal growth temperatures for a list of microorganisms
```linux
./tome predOGT -indir test/proteomes/ -o test/proteomes/predicted_ogt.tsv
```
Then you will get an tab-sperated output file predicted_ogt.tsv with following contents:<br/>
```linux
FileName	predOGT (C)
38_succinivibrio_dextrinosolvens_bacteria.fasta	38.27
95_pyrococcus_horikoshii_archaea.fasta	94.0
69_caldanaerobacter_subterraneus_bacteria.fasta	70.0
```
### 3. Get enzymes for a given ec number. One can specify the temperature range for enzyme temperature optima.
```linux
tome getEnzymes -ec [ec number] -temp_range 0,100 -outdir outdir
```
### 4. Get enzymes for an given ec number. One can specify the temperature range for enzyme temperature optima. In addition, only homologues of the given enzyme sequence would be exported.
```linux
tome getEnzymes -seq seq.fasta -ec [ec number] -temp_range 0,100 -outdir outdir
```

Gang Li<br/>
2018-10-30
