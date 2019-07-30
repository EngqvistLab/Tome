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
import argparse


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
        from tome.core import predOGT
        if args.train: predOGT.train_model()
        else: predOGT.main(args)
    elif args.command == 'getEnzymes':
        from tome.core import getEnzymes
        getEnzymes.main(args,**params)
    else: parser.print_help()

if __name__ == "__main__":
    main()
