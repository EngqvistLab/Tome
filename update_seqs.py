'''
The script retrives protein seqeunces for annotated enzymes from Uniprot.
A fasta file which contains all proteins sequences would be generated.

Gang Li
2018-11-21
'''

import urllib2
import sys

def print_out(line):
    sys.stdout.write(str(line))

def retrive_from_uniprot(uniprot_id):
    # get fasta content for given uniprot_id, return as string format
    url = 'https://www.uniprot.org/uniprot/{0}.fasta'.format(uniprot_id)
    response = urllib2.urlopen(url)
    rec = response.read()
    return rec

def main():
    uniprot_id_file = './external_data/2_unid_growth_temp_mapping.tsv'
    outfile = './external_data/all_enzymes.fasta'

    fhand = open(outfile,'w')
    k = 0
    for line in open(uniprot_id_file):
        if line.startswith('ec'): continue
        uniprot_id = line.split()[1]
        rec = retrive_from_uniprot(uniprot_id)
        fhand.write(rec+'\n')
        k += 1
        if k%1000 == 0: print_out('{0} seqeunces have been downloaded'.format(k))
    print_out('complete!')
    fhand.close()

if __name__ == '__main__': main()
