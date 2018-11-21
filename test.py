import urllib,urllib2,sys
uniprot_id = 'P12345'
url = 'https://www.uniprot.org/uniprot/{0}.fasta'.format(uniprot_id)

response = urllib2.urlopen(url)

page = response.read(200000)

def print_out(line):
    sys.stdout.write(str(line))

print_out(page)
