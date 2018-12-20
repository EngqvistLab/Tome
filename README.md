# Tome: Temperature optima for microorgianisms and enzymes
Estimation of enzyme temperature optima with optimal growth temperature of parent microorganism<br/>
This method has two submodules:
* (1) model for the prediction of optimal growth
temperature;
* (2) find the homologues for a given seqeuence with the same EC nubmer.

## Installation
```linux
pip install tome
```

## Depedences
* ncbi-blast-2.7.1+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* pandas
* Biopython
* numpy
* collections
* sklearn
* Pyhton 2.7 or Python 3


## Usage:
### 1. Prediction of optimal growth temperature
#### 1.1 Prediction of optimal growth temperature for one microorganism
```linux
tome predOGT -fasta test/proteomes/95_pyrococcus_horikoshii_archaea.fasta
```
Then you will get following results:<br/>
```
FileName	predOGT (C)
95_pyrococcus_horikoshii_archaea.fasta	94.0
```

#### 1.2 Predict optimal growth temperatures for a list of microorganisms
```linux
tome predOGT -indir test/proteomes/ -o test/proteomes/predicted_ogt.tsv
```
Then you will get an tab-sperated output file predicted_ogt.tsv with following
contents:<br/>
```
FileName	predOGT (C)
38_succinivibrio_dextrinosolvens_bacteria.fasta	38.27
95_pyrococcus_horikoshii_archaea.fasta	94.0
69_caldanaerobacter_subterraneus_bacteria.fasta	70.0
```
#### 1.3 train the model
In case there would be some warnings due to the versions of sklearn when loading
the model, one can use following command to train the model again:
```linux
./tome predOGT -train
```
Expected output after training is
```
A new model has beed successfully trained.
Model performance:
        RMSE: 2.159489340036136
          r2: 0.9552614628185572
  Pearsnon r:(0.9775886084277753, 0.0)
  Spearman r:SpearmanrResult(correlation=0.93331975456613, pvalue=0.0)

Saving the new model to replace the original one...
Done!
```

### 2. Get enzyme sequences
#### 2.1 Get enzymes for a given ec number.
For example, we want to get the enzymes with EC 3.2.1.1 with a temperature optima
higher 50 °C.
```linux
tome getEnzymes -ec 3.2.1.1 -temp_range 50,200 -outdir test/enzyme_without_seq/
```
Two output files will be generated: test/enzyme_without_seq/3.2.1.1_all.fasta and
test/enzyme_without_seq/3.2.1.1_all.xlsx
3.2.1.1_all.fasta contains all sequences for this EC number. This can be used for
mutisequence alignment with tools like Clustal Omega (https://www.ebi.ac.uk/Tools/msa/clustalo/)
enzyme_without_seq/3.2.1.1_all.xlsx contains following columns:
* uniprot id
* domain: Domain information of source organism (Archaea/Bacteria/Eukaryote)
* organism: name of source organism
* source: if the growth temperature is predicted or experimentally determined
* growth_temp: optimal growth temperature of source organism
* seqeunce: protein sequence

#### 2.2 Get homologous enzymes for an given ec number and sequence.
For example, we want to get all homologs of an enzyme with EC 3.2.1.1
from Photobacterium profundum (OGT = 13°C). We want those homologs with a temperature
optima higher 50 °C. The seqeuence for this enzyme is
```
>Q1Z0D7
MTSLFNTEYASTLSAPSVATNVILHAFDWPYSKVTENAKAIADNGYKAILVSPPLKSFHSKDGTQWWQRYQPQDYRVIDN
QLGNTNDFRTMVEILSLHDIDIYADIVFNHMANESHERSDLNYPNSNIISQYKDKREYFDSIKLFGDLSQPLFSKDDFLS
AFPIKDWKDPWQVQHGRISSGGSDPGLPTLKNNENVVKKQKLYLKALKKIGVKGFRIDAAKHMTLDHIQELCDEDITDGI
HIFGEIITDGGATKEEYELFLQPYLEKTTLGAYDFPLFHTVLDVFNKNASMASLINPYSLGSALENQRAITFAITHDIPN
NDVFLDQVMSEKNEQLAYCYILGRDGGVPLIYTDLDTSGIKNSRGKPRWCEAWNDPIMAKMIHFHNIMHCQPMVIIEQTL
DLLVFSRGHSGIVAINKGKTAVCYKLPAKYSEQDHTEIKEVINMEGVKLSPPSLSTEAGVILQLPAQSCAMLMV
```

```linux
tome getEnzymes -seq test/enzyme_with_seq/test.fasta -ec 3.2.1.1 -temp_range 50,200 -outdir test/enzyme_with_seq/
```
Five output files will be created:
* 3.2.1.1_all.fasta: the same file as described in Section 3
* 3.2.1.1_all.xlsx: the same file as described in Section 3
* blast_3.2.1.1.tsv: blast results in outfmt 6 format
* Q1Z0D7_homologs.fasta: a fasta file which contains seqeunces for all homologs of query enzyme
* Q1Z0D7_homologs.xlsx: an excel file with following columns:
  * uniprot id
  * Identity(%) from blast
  * Coverage(%) from blast
  * domain: Domain information of source organism (Archaea/Bacteria/Eukaryote)
  * organism: name of source organism
  * source: if the growth temperature is predicted or experimentally determined
  * growth_temp: optimal growth temperature of source organism
  * seqeunce: protein sequence

In this test case, 13 homologs with a temperature optima higher than 50 °C were found.

Gang Li<br/>
2018-11-23
