import sys
import os
import pandas as pd
from Bio import SeqIO
from sklearn.externals import joblib
from collections import Counter
from multiprocessing import Pool, cpu_count
import numpy as np

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

def load_model():
    path = os.path.dirname(os.path.realpath(__file__))
    predictor = os.path.join(path,'model/OGT_svr.pkl')
    try:
        model=joblib.load(predictor)
        means,stds,features = load_means_stds(predictor)
    except:
        sys.stdout.write('Failed loading the model. Trying to train the model...')
        from sklearn import svm
        from sklearn.metrics import r2_score
        from scipy.stats import spearmanr,pearsonr

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
        rmse = np.sqrt(MSE(p,Y))
        r2 = r2_score(Y,p)
        r_spearman = spearmanr(p,Y)
        r_pearson = pearsonr(p,Y)

        sys.stdout.write('A new model has beed successfully trained.')
        sys.stdout.write('Model performance:')
        sys.stdout.write('        RMSE: '+ str(rmse))
        sys.stdout.write('          r2: ' + str(r2))
        sys.stdout.write('  Pearsnon r:' + str(r_pearson))
        sys.stdout.write('  Spearman r:' + str(r_spearman))
        sys.stdout.write('')

        # save model
        sys.stdout.write('Saving new model to replace the original one...')
        joblib.dump(model, predictor)

        fea = open(predictor.replace('pkl','f'),'w')
        means, stds = get_standardizer(X)
        fea.write('#Feature_name\tmean\tstd\n')
        for i in range(len(mean)):
            fea.write('{0}\t{1}\t{2}\n'.format(features[i], mean[i], stds[i]))
        fea.close()
        sys.stdout.write('Done!')
        sys.stdout.write('')

    return model,means,stds,features

def do_count(seq):
    dimers = Counter()
    for i in range(len(seq)-1): dimers[seq[i:i+2]] += 1.0
    return dimers


def count_dimer(fasta_file):
    seqs = [str(rec.seq).upper() for rec in SeqIO.parse(fasta_file,'fasta')]

    num_cpus = cpu_count()
    results = Pool(num_cpus).map(do_count, seqs)
    dimers = sum(results, Counter())
    return dict(dimers)

def get_dimer_frequency(fasta_file):
    dimers = count_dimer(fasta_file)
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

def predict(fasta_file,model,means,stds,features):
    dimers_fq = get_dimer_frequency(fasta_file)

    Xs = list()
    for fea in features:
        Xs.append((dimers_fq[fea]-means[fea])/stds[fea])

    Xs = np.array(Xs).reshape([1,len(Xs)])

    pred_ogt = model.predict(Xs)[0]
    return np.around(pred_ogt,decimals=2)


def main():
    args = dict()
    for i in range(len(sys.argv)):
        item = sys.argv[i]
        if item.startswith('-'): args[item] = sys.argv[i+1]

    infile = args.get('-fasta', None)
    indir = args.get('-indir', None)

    if infile is None and indir is None:
        sys.exit('''
        Please check your inputs a gain.
        Usage:
            python predictOGT.py -fasta infile -o outfile
                                or
            python predictOGT.py -indir indir -o outfile

        -fasta: a fasta file that contains all proteins in the proteome
        -indir: a directory that contains all the fasta files

        Gang Li
        2018-06-26\n''')

    if args.get('-o', None) is None: outf = sys.stdout
    else: outf = open(args['-o'], 'w')

    model, means, stds, features = load_model()
    outf.write('FileName\tpredOGT (C)\n')

    if infile is not None:
        pred_ogt = predict(infile,model,means,stds,features)
        outf.write('{0}\t{1}\n'.format(infile.split('/')[-1], pred_ogt))

    else:
        for name in os.listdir(indir):
            if name.startswith('.'): continue
            pred_ogt = predict(os.path.join(indir,name),model,means,stds,features)
            outf.write('{0}\t{1}\n'.format(name, pred_ogt))


if __name__ == '__main__':
    main()
