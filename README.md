# Tome: Temperature optima for microorgianisms and enzymes
A machine learning model for the prediction of optimal growth temperature of microorganisms<br/>
This method has two submodules: (1) model for the prediction of optimal growth temperature; (2) find the homologues for a given seqeuence with the same EC nubmer. 

## Installation

## Depedences
* ncbi-blast-2.7.1+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

## Usage:
### 1. Prediction of optimal growth temperature of a proteome
```python
tome predOGT -fasta proteome.fasta -o outfile
```
### 2. Predict optimal growth temperatures for a list of microorganisms
```python
tome predOGT -indir [dir_to_proteomes] -o outfile
```
### 3. Get enzymes for an given ec number. One can specify the temperature range for enzyme temperature optima
```python
tome getEC -ec [ec number] -temp_range 0,100 -outdir outdir
```
### 4. Get enzymes for an given ec number. One can specify the temperature range for enzyme temperature optima. In addition, only homologues of the given enzyme sequence would be selected.
```python
tome getHomo -seq seq.fasta -ec [ec number] -temp_range 0,100 -outdir outdir
```

Gang Li<br/>
2018-10-30
