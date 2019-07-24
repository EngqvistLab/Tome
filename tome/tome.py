#!/usr/bin/env python

#  This file is part of Tome
#
#  Tome is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Tome is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Tome.  If not, see <https://www.gnu.org/licenses/>


import sys
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from sklearn.externals import joblib
from collections import Counter
from multiprocessing import Pool, cpu_count
import numpy as np
import subprocess
import sqlite3
import argparse


# Esimtation of OGT for organism(s)
################################################################################
def print_out(line):
    sys.stdout.write(str(line)+'\n')

def parse_args():
    args = dict()
    for i in range(len(sys.argv)):
        item = sys.argv[i]
        if item.startswith('-'):
            item = item.replace('-','')
            try:args[item] = sys.argv[i+1]
            except:args[item] = ''

            if item == 'h': args['help'] = ''

    for i in range(len(sys.argv)):
        if 'tome' in sys.argv[i]:
            try:
                if sys.argv[i+1] in ['predOGT','getEnzymes']:
                    args['method'] = sys.argv[i+1]
            except: None
            break
    return args

def load_means_stds(predictor):
    means=dict()
    stds=dict()
    features=list()
    for line in open(predictor.replace('pkl','f'),'r'):
        if line.startswith('#'):continue
        cont=line.split()
        means[cont[0]]=float(cont[1])
        stds[cont[0]]=float(cont[2])
        features.append(cont[0])
    return means,stds,features

def train_model():
    from sklearn import svm
    from sklearn.metrics import r2_score
    from scipy.stats import spearmanr,pearsonr
    from sklearn.metrics import mean_squared_error as MSE

    path = os.path.dirname(os.path.realpath(__file__))
    predictor = os.path.join(path,'model/OGT_svr.pkl')
    def get_standardizer(X):
        mean,std=list(),list()
        for i in range(X.shape[1]):
            mean.append(np.mean(X[:,i]))
            std.append(float(np.var(X[:,i]))**0.5)
        return mean,std

    def standardize(X):
        Xs=np.zeros_like(X)
        n_sample,n_features=X.shape[0],X.shape[1]
        for i in range(n_features):
            Xs[:,i]=(X[:,i]-np.mean(X[:,i]))/float(np.var(X[:,i]))**0.5
        return Xs

    # load training dataset
    trainfile = os.path.join(path,'data/train.csv')
    df = pd.read_csv(trainfile,index_col=0)
    X = df.values[:,:-1]
    Y = df.values[:,-1].ravel()
    features = df.columns[:-1]

    Xs = standardize(X)
    model = svm.SVR(kernel='rbf',C = 64.0, epsilon = 1.0)
    model.fit(Xs,Y)

    # get model performance:
    p = model.predict(Xs)
    rmse = np.sqrt(MSE(Y,p))
    r2 = r2_score(Y,p)
    r_spearman = spearmanr(p,Y)
    r_pearson = pearsonr(p,Y)

    print_out('A new model has beed successfully trained.')
    print_out('Model performance:')
    print_out('        RMSE: '+ str(rmse))
    print_out('          r2: ' + str(r2))
    print_out('  Pearson r:' + str(r_pearson))
    print_out('  Spearman r:' + str(r_spearman))
    print_out('')

    # save model
    print_out('Saving the new model to replace the original one...')
    joblib.dump(model, predictor)

    fea = open(predictor.replace('pkl','f'),'w')
    means, stds = get_standardizer(X)
    fea.write('#Feature_name\tmean\tstd\n')
    for i in range(len(means)):
        fea.write('{0}\t{1}\t{2}\n'.format(features[i], means[i], stds[i]))
    fea.close()
    print_out('Done!')
    print_out('')

def load_model():
    path = os.path.dirname(os.path.realpath(__file__))
    predictor = os.path.join(path,'model/OGT_svr.pkl')
    try:
        model=joblib.load(predictor)
        means,stds,features = load_means_stds(predictor)
    except:
        print_out('Failed loading the model. Trying to train the model...')
        train_model()
        model=joblib.load(predictor)
        means,stds,features = load_means_stds(predictor)

    return model,means,stds,features

def do_count(seq):
    dimers = Counter()
    for i in range(len(seq)-1): dimers[seq[i:i+2]] += 1.0
    return dimers


def count_dimer(fasta_file,p):
    seqs = [str(rec.seq).upper() for rec in SeqIO.parse(fasta_file,'fasta')]

    if p == 0:num_cpus = cpu_count()
    else: num_cpus = p
    results = Pool(num_cpus).map(do_count, seqs)
    dimers = sum(results, Counter())
    return dict(dimers)

def get_dimer_frequency(fasta_file,p):
    dimers = count_dimer(fasta_file,p)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'
    dimers_fq = dict()

    # this is to remove dimers which contains letters other than these 20 amino_acids,
    # like *
    for a1 in amino_acids:
        for a2 in amino_acids:
            dimers_fq[a1+a2] = dimers.get(a1+a2,0.0)
    number_of_aa_in_fasta = sum(dimers_fq.values())
    for key,value in dimers_fq.items(): dimers_fq[key] = value/number_of_aa_in_fasta
    return dimers_fq

def predict(fasta_file,model,means,stds,features,p):
    dimers_fq = get_dimer_frequency(fasta_file,p)

    Xs = list()
    for fea in features:
        Xs.append((dimers_fq[fea]-means[fea])/stds[fea])

    Xs = np.array(Xs).reshape([1,len(Xs)])

    pred_ogt = model.predict(Xs)[0]
    return np.around(pred_ogt,decimals=2)


def predOGT(args):
    infile = args.fasta
    indir = args.indir

    outf = args.out

    model, means, stds, features = load_model()
    outf.write('FileName\tpredOGT (C)\n')

    if infile is not None:
        pred_ogt = predict(infile,model,means,stds,features,args.threads)
        outf.write('{0}\t{1}\n'.format(infile.split('/')[-1], pred_ogt))

    elif indir is not None:
        for name in os.listdir(indir):
            if name.startswith('.'): continue
            if not name.endswith('.fasta'): continue
            pred_ogt = predict(os.path.join(indir,name),model,means,stds,features,args.threads)
            outf.write('{0}\t{1}\n'.format(name, pred_ogt))
    else: sys.exit('Please provide at least a fasta file or a directory that contains \
    a list of fasta files')
    outf.close()

################################################################################

def download_external_data(link):
    realpath = os.path.dirname(os.path.realpath(__file__))
    external_data_path = os.path.join(realpath,'external_data/')
    print_out('Downloading data from {0}'.format(link))
    try:
        subprocess.call(['wget',link,'-P', external_data_path])
        #subprocess.call(['wget {0} -P {1}'.format(link,external_data_path)])
    except:
        file_name = link.split('/')[-1]
        file_name = os.path.join(external_data_path,file_name)
        subprocess.call(['curl',link,'-o',file_name])
        #subprocess.call(['curl {0} -o {}'.format(link,file_name)])

def build_df_from_fetch(cursor):
    # cursor: which has executed a SELECT command

    data = cursor.fetchall()
    cursor.execute('PRAGMA TABLE_INFO(annotation)')

    columns = [item[1] for item in cursor.fetchall()]
    df = pd.DataFrame(data=data,columns=columns)
    df = df.set_index('id',drop=False)
    return df


def select_based_on_class_id(cursor,class_id,class_column,target_column,temps):
    '''
    cursor:        cursor of the sql database which contains all annotation of enzymes
    class_id:      ec number for BRENDA and family id for CAZy
    class_column:  the column name of sequence id. 'ec' for BRENDA and 'family' for
                   CZAy
    target_column: the column name of temperature to be used. Either 'topt' or 'ogt'
    temps:         a list that contains lower and upper bound of temperatures

    return a dataframe with selected enzymes
    '''

    if class_column == 'ec':
        cmd = '''
                SELECT * FROM annotation
                WHERE {0} = '{1}' AND {2} > {3} AND {2} < {4}
        '''.format(class_column,class_id,target_column,temps[0],temps[1])
    if class_column == 'family':
        cmd = '''
                SELECT * FROM annotation
                WHERE {0} LIKE '%{1}%' AND {2} > {3} AND {2} < {4};
        '''.format(class_column,class_id,target_column,temps[0],temps[1])

    cursor.execute(cmd)

    df = build_df_from_fetch(cursor)
    print('{0} enzymes were selected\n'.format(df.shape[0]))
    return df

def build_fasta_from_dataframe(df,class_id,seqid_column,outfasta):
    '''
    df is a dataframe with all annotations of enzymes.
    class_id:      ec number for BRENDA and family id for CAZy
    seqid_column:  the column name of sequence id. 'uniprot_id' for BRENDA and 'genbank'
                   for CAZy
    outfasta: output file

    output file have a format of class_id.fasta
    '''
    print('Build a fasta file for selected enzymes:')

    fhand = open(outfasta,'w')

    for ind in df.index:
        record = SeqRecord(Seq(df.loc[ind,'sequence'],
                       IUPAC.protein),
                   id=df.loc[ind,seqid_column], name="",
                   description="")
        SeqIO.write([record],fhand,'fasta')
        #fhand.write('>{0}\n{1}\n'.format(df.loc[ind,seqid_column],df.loc[ind,'sequence']))
    fhand.close()
    print(outfasta,'\n')

def check_input_fasta(seqfile):
    '''
    Check if there is only one squence provided, if not produce a new fasta with
    the first sequence.

    '''
    seqs = [rec for rec in SeqIO.parse(seqfile,'fasta')]
    if len(seqs) == 1: new_seqfile = seqfile
    else:
        print('Warning: there are more than one query sequence provided. Only the first sequence will be used.')
        new_seqfile = seqfile+'.first.fasta'
        fhand = open(new_seqfile,'w')
        SeqIO.write([seqs[0]],fhand,'fasta')
        fhand.close()
    return new_seqfile


def run_blastp(class_id,seqfile,cpu_num,outdir,evalue):
    '''
    class_id:     ec number for BRENDA and family id for CAZy
    seqfile:      a fasta file that contains the query seqeunce. Only one sequence
                  is allowed. If multiple sequences provided, only the first one
                  would be used
    cpu_num:      number of threads
    outdir:       output directory
    evalue:       evalue cut-off used in blastp

    '''
    dbseq = os.path.join(outdir,'{0}_all.fasta'.format(class_id))
    db = os.path.join(outdir,'db')
    out = os.path.join(outdir,'blast_{}.tsv'.format(class_id))

    cmd = '''makeblastdb -dbtype prot -in {0} -out {1}
    blastp -query {2} -db {1} -outfmt 6 -num_threads {3} -out {4} -evalue {5} -max_hsps 1

    '''.format(dbseq,db,seqfile,cpu_num,out,evalue)
    print('Running blastp:')
    os.system(cmd)
    os.system('rm {0}*'.format(db))


def select_based_on_blast(df,class_id,seqid_column,outdir):
    '''
    df:            a data frame with selected enzymes based on class id and temperature range.
                   the index columns is id
    class_id:      ec number for BRENDA and family id for CAZy
    class_column:  the column name of sequence id. 'ec' for BRENDA and 'family' for
                   CZAy
    outdir:        output directory
    return a dataframe with three more columns: identity,coverage,evalue
    '''
    df = df.set_index(seqid_column,drop=False)

    blastRes = list()
    index = list()
    # blastRes = [ident,coverage,seq]
    # index = [seq_ids]
    blastfile = os.path.join(outdir,'blast_{}.tsv'.format(class_id))

    for line in open(blastfile):
        cont = line.split()
        target = cont[1]
        ident = float(cont[2])
        seq = df.loc[target,'sequence']
        cov = float(cont[3])/len(seq)*100
        eval = float(cont[-2])

        blastRes.append([ident,cov,eval])
        index.append(target)

    os.system('rm {0}'.format(blastfile))
    dfblast = pd.DataFrame(data=blastRes,index=index,columns=['identity(%)','coverage(%)','evalue'])

    dfmerged = pd.merge(df,dfblast,left_index=True,right_index=True,how='inner')
    print('{0} enzymes were selected'.format(dfmerged.shape[0]))
    dfmerged = dfmerged.set_index('id')
    return dfmerged

def check_database(args,params):
    realpath = os.path.dirname(os.path.realpath(__file__))
    dbfile = os.path.join(realpath,params['dbfiles'][args.database])
    if not os.path.isfile(dbfile):
        download_external_data(params['dblinks'][args.database])


def getEnzymes(args,**params):

    check_database(args,params)
    class_column = params['class_column'][args.database]
    seqid_column = params['seqid_column'][args.database]
    path = os.path.dirname(os.path.realpath(__file__))
    dbfile = os.path.join(path,params['dbfiles'][args.database])
    conn = sqlite3.connect(dbfile)
    cursor = conn.cursor()

    print('''
        Step 1: Select enzymes from {0} with given {1} as {2} and with a {3} between {4}
    '''.format(args.database,seqid_column,args.class_id,args.data_type,args.temp_range))

    temps = [float(item) for item in args.temp_range.split(',')]
    df = select_based_on_class_id(cursor,args.class_id,class_column,args.data_type.lower(),temps)
    build_fasta_from_dataframe(df,args.class_id,seqid_column,
    os.path.join(args.outdir,'{0}_{1}_all.fasta'.format(args.class_id,args.data_type.lower())))

    if args.seq is None:
        outtsv = os.path.join(args.outdir,'{0}_{1}_all.tsv'.format(args.class_id,args.data_type.lower()))
        print('Saved results as a tsv file: ',outtsv)
        df.to_csv(outtsv,sep='\t')
        return
    print('''
        Step 2: Search for homologs.
    ''')
    # check input seq
    seqfile = check_input_fasta(args.seq)
    run_blastp(args.class_id,seqfile,args.threads,args.outdir,args.evalue)

    print('''
        Step 3: Select enzymes based on homology.
    ''')
    df = select_based_on_blast(df,args.class_id,seqid_column,args.outdir)
    outtsv = os.path.join(args.outdir,'{0}_{1}_homologs.tsv'.format(args.class_id,args.data_type.lower()))
    print('Saved results as a tsv file: ',outtsv)
    df.to_csv(outtsv,sep='\t')

    build_fasta_from_dataframe(df,args.class_id,seqid_column,
    os.path.join(args.outdir,'{0}_{1}_homologs.fasta'.format(args.class_id,args.data_type.lower())))

def main():
    parser = argparse.ArgumentParser(prog='tome',description='''Tome (Temperature\
    optima for microorganisms and enzymes) is an open source suite for two purposes:\
    (1) predict the optimal growth temperature from proteome sequences (predOGT)
    (2) get homologue enzymes for a given class id (EC number or CAZy family) \
    with/without a sequence (getEnzymes) A detailed list of options can be obtained \
    by calling 'tome predOGT --help'for predOGT or 'tome getEnzymes --help' for getEnzymes''')

    subparsers = parser.add_subparsers(dest='command')


    ############################## OGT arguments ###############################
    parser_ogt = subparsers.add_parser('predOGT',help='Predict the optimal growth\
    temperature from proteomes')

    parser_ogt.add_argument('--fasta',help='a fasta file containing all protein \
    sequence of a proteome.',metavar='',default=None)

    parser_ogt.add_argument('--indir',help='a directory that contains a list of \
    fasta files. Each fasta file is a proteome. Required for the prediction of OGT\
    for a list of microorganisms. Important: Fasta file names much end with .fasta',
    metavar='',default=None)

    parser_ogt.add_argument('--train',help='train the model again',
    action='store_true')

    parser_ogt.add_argument('-o','--out',help='out file name',
    type=argparse.FileType('w', encoding='UTF-8'),default=sys.stdout,metavar='')

    parser_ogt.add_argument('-p','--threads',default=1,help='number of threads \
    used for feature extraction, default is 1. if set to 0, it will use all cpus available',
    metavar='')


    ############################## Topt arguments  #############################
    parser_topt = subparsers.add_parser('getEnzymes',help='Get (homologue) enzymes \
    for a given EC number of CAZy class with/without a sequence')

    parser_topt.add_argument('--database',metavar='',choices=['brenda','cazy'],
    help='the dataset should be used. Should be either brenda or cazy',
    default = 'brenda')

    parser_topt.add_argument('--class_id','--ec',help='EC number or CAZy family.\
    1.1.1.1 for BRENDA, GH1 for CAZy, for instance.',metavar='')

    parser_topt.add_argument('--temp_range',help='the temperature range that target\
    enzymes should be in. For example: 50,100. 50 is lower bound and 100 is upper\
    bound of the temperature.',metavar='')

    parser_topt.add_argument('--data_type',choices=['ogt','topt','OGT','Topt'],
    help = '[OGT or Topt], If OGT, Tome will find enzymes whose OGT of its source\
    organims fall into the temperature range specified by --temp_range. If Topt, \
    Tome will find enzymes whose Topt  fall into the temperature range specified\
    by --temp_range. Default is Topt',metavar='',default = 'Topt')

    parser_topt.add_argument('--seq',help='input fasta file which contains the \
    sequence of the query enzyme. Optional',metavar='',default=None)

    parser_topt.add_argument('--outdir',help='directory for ouput files. Default\
    is current working folder.', metavar='',default='./')

    parser_topt.add_argument('--evalue',help='evalue used in ncbi blastp. Default \
    is 1e-10',metavar='',type=float,default=1e-10)

    parser_topt.add_argument('-p','--threads',default=1,help='number of threads \
    used for blast, default is 1. if set to 0, it will use all cpus available',
    metavar='')

    args = parser.parse_args()



    params = {
    'dbfiles':{'brenda':'external_data/brenda.sql','cazy':'external_data/cazy.sql'},
    'dblinks':{'brenda':'http','cazy':'http'},
    'class_column':{'brenda':'ec','cazy':'family'},
    'seqid_column':{'brenda':'uniprot_id','cazy':'genbank'}
    }

    if args.command == 'predOGT':
        if args.train: train_model()
        else: predOGT(args)
    elif args.command == 'getEnzymes': getEnzymes(args,**params)
    else: parser.print_help()

if __name__ == "__main__":
    main()
