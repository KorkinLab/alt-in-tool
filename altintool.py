import sys
sys.path.append(".")
from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import pickle
import sklearn as sk
from sklearn.externals import joblib
from semi_supervised import SemiSupervisedRF as ssrf
import pkg_resources
from pkg_resources import DistributionNotFound, VersionConflict

dependencies = [
	'pandas',
	'numpy',
	'scipy',
	'sklearn',
	'biopython'
]

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def parse_fasta(fasta_files, ignore_decimals=True):
	fasta_records = {}
	for f in fasta_files:
		record_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
		if ignore_decimals:
			for key in record_dict.keys():
				fasta_records[key.split('.')[0]] = record_dict[key]
		else:
			fasta_records = merge_two_dicts(fasta_records, record_dict)
	return fasta_records

def predict(main_iso, interactor, isoform, fasta, clf):
	main_iso_seq = None
	interactor_seq = None
	isoform_seq = None


	if main_iso in fasta.keys():
		main_iso_seq = fasta[main_iso].seq
	if interactor in fasta.keys():
		interactor_seq = fasta[interactor].seq
	if isoform_seq in fasta.keys():
		isoform_seq = fasta[isoform].seq

	if main_iso_seq is None:
		print "Missing sequence for {}".format(main_iso)	
		return np.nan
	if interactor_seq is None:
		print "Missing sequence for {}".format(interactor)	
		return np.nan
	if isoform is None:
		print "Missing sequence for {}".format(isoform)	
		return np.nan

	os.system("python scripts/extract_features.py {} {} {}".format(main_iso_seq, interactor_seq, isoform_seq))
	features = pickle.load(open('tmp/features.pickle', 'r'))
	features = np.array([features])
	interaction = clf.predict(features)
	os.system('rm tmp/*')
	return interaction[0]

if __name__ == "__main__":
	#pkg_resources.require(dependencies)
	interactors_file = sys.argv[1]
	output = sys.argv[-1]
	fasta_files = []
	if len(sys.argv) > 4:
		fasta_files = sys.argv[2:-1]
	else:
		fasta_files = [sys.argv[2]]
	fasta = parse_fasta(fasta_files)	
	interactors = pd.read_csv(interactors_file, sep='\t', header=0)
	clf = joblib.load('models/semi-supervised-rf.pkl') 
	out = open(output, 'w')
	for idx in interactors.index:
		main_iso = interactors.loc[idx, interactors.columns[0]]	
		interactor = interactors.loc[idx, interactors.columns[1]]	
		isoform = interactors.loc[idx, interactors.columns[2]]	
		interaction = predict(main_iso, interactor, isoform, fasta, clf)
		print interaction
		out.write("{}\t{}\t{}\t{}\n".format(main_iso, interactor, isoform, interaction))
		out.flush()
	#print fasta.keys()[0]#, fasta[fasta.keys()[0]]
	#print fasta[fasta.keys()[0]].seq
 
