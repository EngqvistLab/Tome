'''
The script retrives protein seqeunces for annotated enzymes from Uniprot.
A fasta file which contains all proteins sequences would be generated.

Gang Li
2018-11-21
'''

import urllib2
import sys
from multiprocessing import Pool, cpu_count

def print_out(line):
    sys.stdout.write(str(line)+'\n')

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
    uniprot_ids = list()
    numids = len(uniprot_ids)
    for line in open(uniprot_id_file):
        if line.startswith('ec'): continue
        uniprot_id = line.split()[1]
        uniprot_ids.append(uniprot_id)

    num_cpus = cpu_count()
    print_out('Number of cpus: {0}'.format(num_cpus))
    for i in range(numids/num_cpus+1):
        if i < numids/num_cpus:
            results = Pool(num_cpus).map(retrive_from_uniprot, uniprot_ids[20*i:20*(i+1)])
        else:
            results = Pool(num_cpus).map(retrive_from_uniprot, uniprot_ids[20*i:])

        for res in results: fhand.write(res+'\n')

        k = i*num_cpus
        if k%1000 == 0: print_out('{0} seqeunces have been downloaded'.format(k))
    print_out('complete!')
    fhand.close()

if __name__ == '__main__': main()
