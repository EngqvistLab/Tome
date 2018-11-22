# Tome: Temperature optima for microorgianisms and enzymes
A machine learning model for the prediction of optimal growth temperature of microorganisms<br/>
This method has two submodules: (1) model for the prediction of optimal growth temperature; (2) find the homologues for a given seqeuence with the same EC nubmer. 

## Installation

## Depedences
* ncbi-blast-2.7.1+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)


## Prediction of optimal growth temperature
```python
tome predOGT -fasta proteome.fasta -o outfile
tome predOGT -indir [dir_to_proteomes] -o outfile
```

## Get homologues for a given enzyme sequence.
```python
tome getEC -ec [ec number] -trg 0,100 -outdir outdir
tome getHomo -seq seq.fasta -ec [ec number] -temprange 0,100 -outdir outdir
```

Gang Li<br/>
2018-10-30
