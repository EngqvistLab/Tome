# Construct a fasta file to store all the annotated sequences.
# The fasta file has the following format:
# >uniprot_id
# sequence
#
# Gang Li
# 2018-11-22

from Bio import SeqIO
import os
import sys

def print_out(line):
    sys.stdout.write(str(line)+'\n')

def load_uniprot_ids(annofile):

    uniprot_ids = dict()

    for line in open(annofile):
        if line.startswith('ec'): continue
        uniprot_id = line.split()[1]
        uniprot_ids[uniprot_id] = True

    return uniprot_ids


def load_seqs_for_selected_ids(fa_file,uniprot_ids,seqs=dict()):
    for rec in SeqIO.parse(fa_file,'fasta'):
        uniprot_id = rec.id.split(';')[0]
        if uniprot_ids.get(uniprot_id,False): seqs[uniprot_id] = rec.seq
    return seqs

def main():
    annofile = '../external_data/2_unid_growth_temp_mapping.tsv'
    uniprot_ids = load_uniprot_ids(annofile)

    print_out('Number of sequences in annotation file: {0}'.format(len(uniprot_ids)))

    fa_dir = '../external_data/fasta/'

    seqs = dict()
    for name in os.listdir(fa_dir):
        if not name.endswith('fasta'): continue
        fa_file = os.path.join(fa_dir,name)
        seqs = load_seqs_for_selected_ids(fa_file,uniprot_ids)

    print_out('Number of sequences collected: {0}'.format(len(seqs)))

    outfile = '../external_data/all_enzyme_sequences.fasta'
    fhand = open(outfile,'w')
    for id,seq in seqs.items(): fhand.write('>{0}\n{1}\n'.format(id,seq))
    fhand.close()


if __name__ == '__main__':
    main()
