import os

for i in range(10000):
    try: os.system('python update_seqs.py')
    except: None
