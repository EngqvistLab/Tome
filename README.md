# Tome: Temperature optima for microorgianisms and enzymes
Tome (Temperature optima for microorganisms and enzymes) is an open source suite
for two fundamental applications:
  * predict the optimal growth temperature from proteome sequences
  * get homologue enzymes for a given ec number with/without a sequence that have
  a temperature optima in a specified range.

## System
* macOS
* Linux

## Installation
##### (1). Download tome package
##### (2). Open your terminal
##### (3). Change directory to the tome package
```linux
cd [directory to tome, where setup.py is]
```
##### (4). Run following command
```linux
pip install -e .
```
##### (5) Now you can use 'tome' via command line.
There is a folder named 'test' in the package. One can use the instructions in
'Usage' section to test the package.

## Depedences
* ncbi-blast-2.7.1+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
(This is only mandatory for 'tome predOGT -seq')
* pandas
* Biopython
* numpy
* sklearn
* Pyhton 2.7 or Python 3


## Usage:
### 1. Prediction of optimal growth temperature
#### 1.1 Prediction of optimal growth temperature for one microorganism
```linux
tome predOGT --fasta test/proteomes/95_pyrococcus_horikoshii_archaea.fasta
```
Then you will get following results:<br/>
```
FileName	predOGT (C)
pyrococcus_horikoshii.fasta	94.0
```

#### 1.2 Predict optimal growth temperatures for a list of microorganisms. Fasta files must end with .fasta
```linux
tome predOGT --indir test/proteomes/ -o test/proteomes/predicted_ogt.tsv
```
Then you will get an tab-seperated output file predicted_ogt.tsv with following
contents:<br/>
```
FileName	predOGT (C)
succinivibrio_dextrinosolvens.fasta	38.27
pyrococcus_horikoshii.fasta	94.0
caldanaerobacter_subterraneus.fasta	70.0
```
#### 1.3 train the model
In case there would be some warnings due to the versions of sklearn when loading
the model, one can use following command to train the model again:
```linux
tome predOGT --train
```
Expected output after training is
```
A new model has beed successfully trained.
Model performance:
        RMSE: 2.159489340036136
          r2: 0.9552614628185572
  Pearson r:(0.9775886084277753, 0.0)
  Spearman r:SpearmanrResult(correlation=0.93331975456613, pvalue=0.0)

Saving the new model to replace the original one...
Done!
```

### 2. Get enzyme sequences
The first time this command is run two files totaling ~2.5 GB will be downloaded from the Zenodo data repository (https://zenodo.org/record/2539114#.XDyAUc_0l0s). These files contain the enzyme annotation data. 
#### 2.1 Get enzymes for a given ec number.
For example, we want to get the enzymes with EC 3.2.1.1 with a temperature optima
higher 50 °C.
```linux
tome getEnzymes --ec 3.2.1.1 --temp_range 50,200 --data_type Topt --outdir test/enzyme_without_seq/
```

Two output files will be generated: test/enzyme_without_seq/3.2.1.1_all.fasta and
test/enzyme_without_seq/3.2.1.1_all.tsv
3.2.1.1_all.fasta contains all sequences for this EC number. This can be used for
mutisequence alignment with tools like Clustal Omega (https://www.ebi.ac.uk/Tools/msa/clustalo/)
enzyme_without_seq/3.2.1.1_all.tsv contains following columns:
* uniprot id
* domain: Domain information of source organism (Archaea/Bacteria/Eukaryote)
* organism: name of source organism
* ogt: optimal growth temperature of source organism
* ogt_source: if the growth temperature is predicted or experimentally determined
* topt: temperature optima of the enzyme
* topt_source: if topt is predicted or experimentally determined
* seqeunce: protein sequence

One can also use following command to find enzymes from organisms with a OGT fall
into the temperature range specified with --temp_range. The output files have the same
format as above described.
```linux
tome getEnzymes --ec 3.2.1.1 --temp_range 50,200 --data_type OGT --outdir test/enzyme_without_seq/
```

#### 2.2 Get homologous enzymes for an given ec number and sequence.
For example, we want to get all homologs of an enzyme with EC 3.2.1.1
from Photobacterium profundum (OGT = 13°C). We want those homologs with a temperature
optima higher 50 °C. The sequence for this enzyme is
```
>Q1Z0D7
MTSLFNTEYASTLSAPSVATNVILHAFDWPYSKVTENAKAIADNGYKAILVSPPLKSFHSKDGTQWWQRYQPQDYRVIDN
QLGNTNDFRTMVEILSLHDIDIYADIVFNHMANESHERSDLNYPNSNIISQYKDKREYFDSIKLFGDLSQPLFSKDDFLS
AFPIKDWKDPWQVQHGRISSGGSDPGLPTLKNNENVVKKQKLYLKALKKIGVKGFRIDAAKHMTLDHIQELCDEDITDGI
HIFGEIITDGGATKEEYELFLQPYLEKTTLGAYDFPLFHTVLDVFNKNASMASLINPYSLGSALENQRAITFAITHDIPN
NDVFLDQVMSEKNEQLAYCYILGRDGGVPLIYTDLDTSGIKNSRGKPRWCEAWNDPIMAKMIHFHNIMHCQPMVIIEQTL
DLLVFSRGHSGIVAINKGKTAVCYKLPAKYSEQDHTEIKEVINMEGVKLSPPSLSTEAGVILQLPAQSCAMLMV
```
There should be only one sequence in the fasta file. If more than 1 sequence is provided,
only the first sequence would be used.
```linux
tome getEnzymes --seq test/enzyme_with_seq/test.fasta --ec 3.2.1.1 --temp_range 50,200 --data_type Topt --outdir test/enzyme_with_seq/
```
Two output files will be created:
* Q1Z0D7_homologs.fasta: a fasta file which contains sequences for all homologs of query enzyme
* Q1Z0D7_homologs.tsv: a tab-seperated file with following columns:
  * uniprot id
  * domain: Domain information of source organism (Archaea/Bacteria/Eukaryote)
  * organism: name of source organism
  * ogt: optimal growth temperature of source organism
  * ogt_source: if the growth temperature is predicted or experimentally determined
  * topt: temperature optima of the enzyme
  * topt_source: if topt is predicted or experimentally determined
  * Identity(%) from blast
  * Coverage(%) from blast
  * seqeunce: protein sequence

In this test case, 44 homologs with a temperature optima higher than 50 °C were found.

One can also use following command to find homologous enzymes from organisms with a OGT fall
into the temperature range specified with --temp_range. The output files have the same
format as above described.
```linux
tome getEnzymes --seq test/enzyme_with_seq/test.fasta --ec 3.2.1.1 --temp_range 50,200 --data_type OGT --outdir test/enzyme_with_seq/
```
In this case, 13 homologs from organisms with a OGT higher than 50 °C were found

## Help:
Use following commands you can get detailed information about the arguments of tome.
```linux
tome --help/-h
tome predOGT --help/-h
tome getEnzymes --help/-h
```
Or you can directly contact
Martin Engqvist: <martin.engqvist@chalmers.se> or Gang Li: <gangl@chalmers.se><br/>
<br/>
Gang Li<br/>
