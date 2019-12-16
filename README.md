# Tome: Temperature optima for microorganisms and enzymes
Tome (Temperature optima for microorganisms and enzymes) is an open source suite
for two fundamental applications:
  * predict the optimal growth temperature from proteome sequences
  * get homologue enzymes for a given ec number with/without a sequence that have
  a temperature optima in a specified range.

## Citation
If you have used `tome getEnzymes` in Tome v1.0 or `tome predOGT`, please cite
*Li, G., Rabe, K. S., Nielsen, J., & Engqvist, M. K. M. (2019). Machine learning applied to predicting microorganism growth temperatures and enzyme catalytic optima. ACS Synthetic Biology, 8, 1411–1420.*

If you have `tome getEnzymes` in Tome v2.0, please cite 
*Ref to be updated*

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
Get help message of `tome` with 
```
tome -h
```
```
usage: tome [-h] {predOGT,getEnzymes} ...

Tome (Temperature optima for microorganisms and enzymes) is an open source suite for two purposes: (1) predict the
optimal growth temperature from proteome sequences (predOGT) (2) get homologue enzymes for a given class id (EC number or
CAZy family) with/without a sequence (getEnzymes) A detailed list of options can be obtained by calling 'tome predOGT
--help'for predOGT or 'tome getEnzymes --help' for getEnzymes

positional arguments:
  {predOGT,getEnzymes}
    predOGT             Predict the optimal growth temperature from proteomes
    getEnzymes          Get (homologue) enzymes for a given EC number of CAZy class with/without a sequence

optional arguments:
  -h, --help            show this help message and exit
```

### 1. Prediction of optimal growth temperature
Get help message of `tome predOGT` with 
```
tome predOGT -h
```

```
usage: tome predOGT [-h] [--fasta] [--indir] [--train] [-o] [-p]

optional arguments:
  -h, --help       show this help message and exit
  --fasta          a fasta file containing all protein sequence of a proteome.
  --indir          a directory that contains a list of fasta files. Each fasta file is a proteome. Required for the
                   prediction of OGT for a list of microorganisms. Important: Fasta file names much end with .fasta
  --train          train the model again
  -o , --out       out file name
  -p , --threads   number of threads used for feature extraction, default is 1. if set to 0, it will use all cpus
                   available
```

#### Case 1.1 Prediction of optimal growth temperature for one microorganism
```linux
tome predOGT --fasta test/proteomes/pyrococcus_horikoshii.fasta
```
Then you will get following results:<br/>
```
FileName	predOGT (C)
pyrococcus_horikoshii.fasta	94.0
```

#### Case 1.2 Predict optimal growth temperatures for a list of microorganisms. Fasta files must end with .fasta
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
#### train the model
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
The first time `tome getEnzymes --database brenda` used, a SQLite3 database ~3.5 GB will be downloaded.  
The first time `tome getEnzymes --database cazy` used, a SQLite3 database ~601 MB will be downloaded.  
Those two databases are deposited on Zenodo (https://zenodo.org/record/3578468#.XffgbpP0nOQ). These files contain the enzyme annotation data.

If automatic download failed, then one can manually download following two files from Zenodo and put them into `tome/external_data/`: 
BRENDA: https://zenodo.org/record/3578468/files/brenda.sql
CAZy  : https://zenodo.org/record/3578468/files/cazy.sql
  
Get help message of `tome getEnzymes` with
```
tome getEnzymes -h
```

```
usage: tome getEnzymes [-h] [--database] [--class_id] [--temp_range] [--data_type] [--seq] [--outdir] [--evalue] [-p]

optional arguments:
  -h, --help          show this help message and exit
  --database          the dataset should be used. Should be either brenda or cazy
  --class_id , --ec   EC number or CAZy family. 1.1.1.1 for BRENDA, GH1 for CAZy, for instance.
  --temp_range        the temperature range that target enzymes should be in. For example: 50,100. 50 is lower bound and
                      100 is upper bound of the temperature.
  --data_type         [OGT or Topt], If OGT, Tome will find enzymes whose OGT of its source organims fall into the
                      temperature range specified by --temp_range. If Topt, Tome will find enzymes whose Topt fall into
                      the temperature range specified by --temp_range. Default is Topt
  --seq               input fasta file which contains the sequence of the query enzyme. Optional
  --outdir            directory for ouput files. Default is current working folder.
  --evalue            evalue used in ncbi blastp. Default is 1e-10
  -p , --threads      number of threads used for blast, default is 1. if set to 0, it will use all cpus available
```


#### Brenda Case 2.1 Get enzymes for a given ec number.
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

#### Brenda Case 2.2 Get homologous enzymes for an given ec number and sequence.
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

#### CAZy case 2.4 Get enzymes for a given CAZy family ID.
```
tome getEnzymes --database cazy --class_id GH1 --temp_range 50,60 --data_type Topt --outdir test/cazy_enzyme_without_seq/
```
Similar as for BRENDA, This will find all enzymes in CAZy database that (1) belongs to GH1 family, (2) have a predicted Topt between 50 and 60 °C.

#### CAZy case 2.4 Get homologous enzymes for a given CAZy family ID.
```
tome getEnzymes --database cazy --class_id GH1 --temp_range 50,60 --data_type Topt --outdir test/cazy_enzyme_with_seq/ --seq test/cazy_enzyme_with_seq/test.fasta
```
Similar as for BRENDA, This will find all enzymes in CAZy database that (1) belongs to GH1 family, (2) have a predicted Topt between 50 and 60 °C, (3) homologous to the given sequence in `test.fasta`
