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


################################################################################

def download_external_data(link):
    realpath = os.path.dirname(os.path.realpath(__file__))
    external_data_path = os.path.join(realpath,'external_data/')
    print('Downloading data from {0}'.format(link))
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
    df = df.set_index('id')
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


def run_blastp(class_id,seqfile,cpu_num,outdir,evalue,data_type):
    '''
    class_id:     ec number for BRENDA and family id for CAZy
    seqfile:      a fasta file that contains the query seqeunce. Only one sequence
                  is allowed. If multiple sequences provided, only the first one
                  would be used
    cpu_num:      number of threads
    outdir:       output directory
    evalue:       evalue cut-off used in blastp
    data_type:    topt or ogt

    '''
    dbseq = os.path.join(outdir,'{0}_{1}_all.fasta'.format(class_id,data_type))
    db = os.path.join(outdir,'db')
    out = os.path.join(outdir,'blast_{0}_{1}.tsv'.format(class_id,data_type))

    cmd = '''makeblastdb -dbtype prot -in {0} -out {1}
    blastp -query {2} -db {1} -outfmt 6 -num_threads {3} -out {4} -evalue {5} -max_hsps 1

    '''.format(dbseq,db,seqfile,cpu_num,out,evalue)
    print('Running blastp:')
    os.system(cmd)
    os.system('rm {0}*'.format(db))


def select_based_on_blast(df,class_id,seqid_column,outdir,data_type):
    '''
    df:            a data frame with selected enzymes based on class id and temperature range.
                   the index columns is id
    class_id:      ec number for BRENDA and family id for CAZy
    class_column:  the column name of sequence id. 'ec' for BRENDA and 'family' for
                   CZAy
    outdir:        output directory
    return a dataframe with three more columns: identity,coverage,evalue
    '''
    df['id'] = df.index
    df = df.set_index(seqid_column,drop=False)

    blastRes = list()
    index = list()
    # blastRes = [ident,coverage,seq]
    # index = [seq_ids]
    blastfile = os.path.join(outdir,'blast_{0}_{1}.tsv'.format(class_id,data_type))

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
    realpath = os.path.dirname(os.path.realpath(__file__)).replace('core','')
    dbfile = os.path.join(realpath,params['dbfiles'][args.database])
    if not os.path.isfile(dbfile):
        download_external_data(params['dblinks'][args.database])


def main(args,**params):

    check_database(args,params)
    class_column = params['class_column'][args.database]
    seqid_column = params['seqid_column'][args.database]
    path = os.path.dirname(os.path.realpath(__file__)).replace('core','')
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
    run_blastp(args.class_id,seqfile,args.threads,args.outdir,args.evalue,args.data_type.lower())

    print('''
        Step 3: Select enzymes based on homology.
    ''')
    df = select_based_on_blast(df,args.class_id,seqid_column,args.outdir,args.data_type.lower())
    outtsv = os.path.join(args.outdir,'{0}_{1}_homologs.tsv'.format(args.class_id,args.data_type.lower()))
    print('Saved results as a tsv file: ',outtsv)
    df.to_csv(outtsv,sep='\t')

    build_fasta_from_dataframe(df,args.class_id,seqid_column,
    os.path.join(args.outdir,'{0}_{1}_homologs.fasta'.format(args.class_id,args.data_type.lower())))
