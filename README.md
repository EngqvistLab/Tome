# light
A machine learning model for the prediction of optimal growth temperature of microorganisms<br/>
This method has two submodules: (1) model for the prediction of optimal growth temperature; (2) find the homologues for a given seqeuence with the same EC nubmer. 

## Installation

## Depedences


## Prediction of optimal growth temperature
```python
light predOGT -infile proteome.fasta -out outfile
light predOGT -indir [dir_to_proteomes] -out outfile
```

## Get homologues for a given enzyme sequence.
```python
light getHomo -infile seq.fasta -ec [EC number] -out oufile
```

Gang Li<br/>
2018-10-30
