'''
The script retrives protein seqeunces for annotated enzymes from Uniprot.
A fasta file which contains all proteins sequences would be generated.

Gang Li
2018-11-21
'''

import urllib2
import sys
from multiprocessing import Pool, cpu_count
import os
from Bio import SeqIO

def print_out(line):
    sys.stdout.write(str(line)+'\n')

def retrive_from_uniprot(uniprot_id):
    # get fasta content for given uniprot_id, return as string format
    url = 'https://www.uniprot.org/uniprot/{0}.fasta'.format(uniprot_id)
    response = urllib2.urlopen(url)
    rec = response.read()
    return rec

def load_uniprot_ids(infile,existed_ids):
    uniprot_ids = list()
    for line in open(infile):
        if line.startswith('ec'): continue
        uniprot_id = line.split()[1]
        if not existed_ids.get(uniprot_id,False):uniprot_ids.append(uniprot_id)
    return uniprot_ids

def get_existed_sequence_ids(outdir):
    existed_ids = dict()
    for name in os.listdir(outdir):
        if not name.endswith('fasta'): continue
        for rec in SeqIO.parse(os.path.join(outdir,name),'fasta'):
            uniprot_id = rec.id.split('|')[1]
            existed_ids[uniprot_id] = True
    return existed_ids

def get_out_file(outdir):
    k = 0
    for name in os.listdir(outdir):
        if not name.endswith('fasta'): continue
        k += 1
    outfile = os.path.join(outdir,'all_enzymes_{0}.fasta'.format(k))
    return outfile


def main():
    uniprot_id_file = './external_data/2_unid_growth_temp_mapping.tsv'
    outdir = './external_data/'

    existed_ids = get_existed_sequence_ids(outdir)
    uniprot_ids = load_uniprot_ids(uniprot_id_file,existed_ids)

    outfile = get_out_file(outdir)
    fhand = open(outfile,'w')


    numids = len(uniprot_ids)
    print_out('Number of total sequences to be downloaded:{0}'.format(numids))

    num_cpus = cpu_count()-1
    print_out('Number of cpus: {0}'.format(num_cpus))
    for i in range(numids/100+1):
        if i < numids/100:
            results = Pool(num_cpus).map(retrive_from_uniprot, uniprot_ids[100*i:100*(i+1)])
        else:
            results = Pool(num_cpus).map(retrive_from_uniprot, uniprot_ids[100*i:])

        for res in results: fhand.write(res+'\n')

        k = i*100
        print_out('{0} seqeunces have been downloaded'.format(k))
    print_out('complete!')
    fhand.close()

if __name__ == '__main__': main()
