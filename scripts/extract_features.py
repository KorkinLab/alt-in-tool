import os
import pickle
import pandas as pd
import numpy as np
import operator
import scipy as sp
import operator
import itertools
import subprocess
from sklearn import datasets
from sklearn.preprocessing import Imputer
from sklearn.metrics.pairwise import pairwise_distances
from Bio.PDB import *
from Bio import pairwise2
from sklearn.preprocessing import Imputer
from Bio.pairwise2 import format_alignment
from timeout import timeout
import sys
sys.path.append('.')
sys.path.append('./scripts')


### BioPython's pairwise sequence alignment. To change maximum allowed wait time change @timeout(seconds) ####
@timeout(3)
def align(main_iso_seq, isoform_seq, ratio):
	return pairwise2.align.globalmd(main_iso_seq, isoform_seq, 5, -4, -10, -10, -10 * ratio, -1.5)

class Config:
	emboss_path_ = ''
	interpro_path_ = ''

	def __init__(self, config_path):
		config_file = open(config_path)
		emboss_line = config_file.readline()
		
		self.emboss_path_ = emboss_line.split('=')[-1].split('\'')[1]
		if self.emboss_path_ == '':
		    self.emboss_path_ = os.environ['EMBOSS_PATH']
		if len(self.emboss_path_) > 0:
			if self.emboss_path_[-1] != '\\' or self.emboss_path_[-1] != '/':
				self.emboss_path_ += '/'
		   
		interpro_line = config_file.readline()
		self.interpro_path_ = interpro_line.split('=')[-1].split('\'')[1]
		if self.interpro_path_ == '':
		    self.interpro_path_ = os.environ['INTERPRO_PATH']
		if len(self.interpro_path_) > 0:
			if self.interpro_path_[-1] != '\\' or self.interpro_path_[-1] != '/':
				self.interpro_path_ += '/'
		print "EMBOSS PATH: {}".format(self.emboss_path_)
		print "INTERPRO PATH: {}".format(self.interpro_path_)

#calls pepstat, returns output filename (pepstat must be found by calling "pepstats")
def runpepstats(seq, config):
	fileloc = './tmp/seq'
	outfile = './tmp/{}.pepstats'.format('pepstats')
	f = open(fileloc, 'w+')
	f.write(seq)
	f.close()
	#print 'Running pepstats ...'
	os.system("{}pepstats -sequence ".format(config.emboss_path_)+fileloc+" -outfile "+outfile+" -auto") 
	#print 'Pepstats finished'
	return outfile

#Parse one amino acid pepstats output for matrix
def parsepepstats(stat):
	lines = stat
	result=list(range(0,106))
	result[0]	   = lines[0].split()[2]
	result[1]	   = lines[2].split()[3]
	result[2]       = lines[2].split()[6]
	result[3]       = lines[3].split()[4]
	result[4]       = lines[3].split()[7]
	result[5]       = lines[4].split()[3]
	result[6]       = lines[5].split()[5]
	result[7]       = lines[5].split()[7]
	result[8]       = lines[6].split()[5]
	result[9]       = lines[6].split()[7]
	for i in range(10,36):
		result[i]	= lines[i].split()[3]
	for i in range(10,36):
		result[i+26] 	= lines[i].split()[4]
	for i in range(10,36):
		result[i+52] 	= lines[i].split()[5]
	for i in range(38,47):
		result[i+50] = lines[i].split()[2]
	for i in range(38,47):
		result[i+59] = lines[i].split()[3]
	return (result)

#file run, consideration whether its just one protein file or multiple
def extract_biochemical_features(sequence, config):
	name = ['gene_name','mw','res','arw','charge','iep','A280MECR','A280MECC','A280ECR','A280ECC',\
	'A_N','B_N','C_N','D_N','E_N','F_N','G_N','H_N','I_N','J_N','K_N','L_N','M_N','N_N','O_N','P_N','Q_N','R_N','S_N','T_N','U_N','V_N','W_N','X_N','Y_N','Z_N',\
	'A_M','B_M','C_M','D_M','E_M','F_M','G_M','H_M','I_M','J_M','K_M','L_M','M_M','N_M','O_M','P_M','Q_M','R_M','S_M','T_M','U_M','V_M','W_M','X_M','Y_M','Z_M',\
	'A_D','B_D','C_D','D_D','E_D','F_D','G_D','H_D','I_D','J_D','K_D','L_D','M_D','N_D','O_D','P_D','Q_D','R_D','S_D','T_D','U_D','V_D','W_D','X_D','Y_D','Z_D',\
	'Tiny_N','Small_N','Aliphatic_N','Aromatic_N','Non-polar_N','Polar_N','Charged_N','Basic_N','Acidic_N',\
	'Tiny_M','Small_M','Aliphatic_M','Aromatic_M','Non-polar_M','Polar_M','Charged_M','Basic_M','Acidic_M']
	#group="group_177"
	i=runpepstats(sequence, config)
	f = open(i,"r")
	line = f.readlines()
	f.close()
	iteration=int(len(line)/48)
	for p in range(0,iteration):
		if p == 0:
			sample=line[0:48]
			tmp=parsepepstats(sample)
			output=pd.DataFrame([tmp])
		if p > 1:
			start=p*48
			end=start+48
			sample=line[start:end]
			tmp=parsepepstats(sample)
			output=output.append([tmp])
		#print(p)
	output.columns=name
	output = output.drop(['gene_name'], axis=1)
	output = output.apply(pd.to_numeric)
	#output.to_csv(group+".biochemical_feature.csv",sep=",",index=False)
	return np.squeeze(output.as_matrix())

def get_potential(data, domain1, domain2):
	domain1 = str(domain1)
	domain2 = str(domain2)
	if data.has_key(domain1):
		if data[domain1].has_key(domain2):
			return data[domain1][domain2]
	return 0

def load_stat_pot_data():
	data = pd.read_csv('./data/curated_statistical_potentials.csv', header=None)
	output = {}
	for i in data.index:
		pair, potential = data.loc[i, data.columns]
		domain1, domain2 = pair.split('-')
		if not output.has_key(domain1):
			output[domain1] = {}
		output[domain1][domain2] = potential
		if not output.has_key(domain2):
			output[domain2] = {}
			output[domain2][domain1] = potential
	return output

def get_max_potential(data, domain):
	domain = str(domain)
	max_pot = 0
	if data.has_key(domain):
		for key in data[domain].keys():
			pot = data[domain][key]
			if pot > min(max_pot):
				max_pot = pot
	return max_pot

def get_max_potentials(data, domains):
	max_pots = np.array([0, 0])
	for domain in domains:
		max_pot = get_max_potential(data, domain)
		if max_pot > min(max_pots):
			max_pots[np.argmin(max_pots)] = max_pot
	return max_pots

def top_3_interact_pot(data, domains1, domains2):
	inter_pot = np.array([0, 0, 0])
	for dom1 in domains1:
		for dom2 in domains2:
			pot = get_potential(data, dom1, dom2)
			if pot > min(inter_pot):
				inter_pot[np.argmin(inter_pot)] = pot
	return inter_pot



def make_fasta(seq, fname=None):
	if fname is None:
		fname = 'tmp/seq.fasta'
	else:
		fname = 'tmp/{}'.format(fname)
	f = open(fname, 'w')
	f.write('>SEQ\n')
	f.write(seq)
	f.close()
	return fname


def extract_interpro(seq, output_file, config):
	fasta = make_fasta(seq)
	#output = 'tmp/superfam_domains.interpro'
	#print('{}/interproscan.sh'.format(interpro_path) + ' -i ' + fasta + ' --output-file-base ' + output_file + ' -f ' + 'tsv ' + '--disable-precalc ' + '-appl ' + 'SUPERFAMILY-1.75')
	#cmd = subprocess.Popen(['{}/interproscan.sh'.format(interpro_path), '-i', fasta, '--output-file-base', output_file, '-f', 'tsv', '--disable-precalc', '-appl', 'SUPERFAMILY-1.75'])
	os.system('{}interproscan.sh'.format(config.interpro_path_) + ' -i ' + fasta + ' --output-file-base ' + output_file + ' -f ' + 'tsv ' + '--disable-precalc ' + '-appl ' + 'SUPERFAMILY-1.75')
	#interproscan.sh -i [input fasta file] --output-file-base [file base for output] -f tsv --disable-precalc -appl SUPERFAMILY-1.75	

def extract_domains(seq, config):
	#main_iso_interpro = 'tmp/main_iso.interpro'
	isoform_interpro = 'tmp/iso.interpro'
	#extract_interpro(main_iso_seq, main_iso_interpro)
	extract_interpro(isoform_seq, isoform_interpro, config)
	isoform_interpro += '.tsv'
	#subprocess.call(['parse_domain_information.py', isoform_interpro])#main_iso_interpro, isoform_interpro])
	os.system('python scripts/parse_domain_information.py ' + isoform_interpro)
	domains = pickle.load(open('tmp/parsed.pickle', 'r'))
	binding_sites = pickle.load(open('tmp/binding_sites.pickle', 'r'))
	os.system('rm {}'.format(isoform_interpro))
	return domains, binding_sites


def gap_lengths(aligned_seq):
        bits = (aligned_seq == '-')
        return np.array([sum(g) for b, g in itertools.groupby(bits) if b])


def extract_new_features(main_iso_seq, inter_seq, isoform_seq, config):
		
	#changes in domains
	domains_lost = 0
	domains_procured = 0
	domains_changed = 0

	main_iso_domains = []
	isoform_domains = []
	interactor_domains = []
	main_iso_binding_sites = []

	domains, binding_sites = extract_domains(main_iso_seq, config)
	if len(domains) > 0:
		main_iso_domains = domains[domains.keys()[0]]
		main_iso_binding_sites = binding_sites[binding_sites.keys()[0]]

	domains, binding_sites = extract_domains(isoform_seq, config)
	if len(domains) > 0:
		isoform_domains = domains[domains.keys()[0]]
		isoform_binding_sites = binding_sites[binding_sites.keys()[0]]

	interactor_domains, binding_sites = extract_domains(inter_seq, config)
	if len(domains) > 0:
		interactor_domains = domains[domains.keys()[0]]
		interactor_binding_sites = binding_sites[binding_sites.keys()[0]]


	for main_domain in main_iso_domains:
		if main_domain not in isoform_domains:
			domains_lost += 1

	for iso_domain in isoform_domains:
		if iso_domain not in main_domain:
			domains_procured += 1


	#Length
	length_delta_ratio = float(len(main_iso_seq) - len(isoform_seq))/len(main_iso_seq)
	length_delta = len(main_iso_seq) - len(isoform_seq)

	#Alignment
	n_termini = None
	c_termini = None
	alignment = None
	gap_sizes = None
	s1 = None
	s2 = None


	if len(main_iso_seq) >= len(isoform_seq):
		ratio = float(len(main_iso_seq)) / len(isoform_seq)
		alignments = align(main_iso_seq, isoform_seq, ratio)#pairwise2.align.globalmd(main_iso_seq, isoform_seq, 5, -4, -10, -10, -10 * ratio, -1.5)
		alignment = alignments[0]
		#print alignment

		s1 = np.array(list(alignment[0]))
		s2 = np.array(list(alignment[1]))

		gap_sizes = gap_lengths(s2)

		n_termini = 1 if sum(s2[:20] == '-') == 0 else (-1 if sum(s2[:20] == '-') > 12 else 0)
		c_termini = 1 if sum(s2[-20:] == '-') == 0 else (-1 if sum(s2[-20:] == '-') > 12 else 0)


	if len(main_iso_seq) < len(isoform_seq):
		ratio = float(len(isoform_seq)) / len(main_iso_seq)
		alignments = align(isoform_seq, main_iso_seq, ratio)#pairwise2.align.localmd(isoform_seq, main_iso_seq, 5, -4, -100, -10, -10 * ratio, -1.5)
		alignment = alignments[0]
		#print alignment

		s1 = np.array(alignment[0].split())
		s2 = np.array(alignment[1].split())

		gap_sizes = gap_lengths(s2)

		n_termini = 1 if sum(s2[:20] == '-') == 0 else (2 if sum(s2[:20] == '-') > 12 else 0)
		c_termini = 1 if sum(s2[-20:] == '-') == 0 else (2 if sum(s2[-20:] == '-') > 12 else 0)

	domain_or_linker = 0
	if len(main_iso_binding_sites) > 0:
		main_iso_domain_sites = main_iso_binding_sites
		for domain, binding_sites in main_iso_domain_sites.items():
			start = int(binding_sites[0])
			end = int(binding_sites[1])
			if '-' in alignment[start:end+1]:
				domain_or_linker = 1
				break

	max_gap_size = 0 if len(gap_sizes) == 0 else max(gap_sizes)
	mean_gap_size = 0 if len(gap_sizes) == 0 else np.mean(gap_sizes)
	large_change_num = sum(gap_sizes >= 10)
	small_change_num = sum(gap_sizes < 10)
	
	stat_pot_data = load_stat_pot_data()
	interact_pot_main_iso = top_3_interact_pot(stat_pot_data, main_iso_domains, interactor_domains)
	max_pot_main_iso = get_max_potentials(stat_pot_data, main_iso_domains)
	interact_pot_iso = top_3_interact_pot(stat_pot_data, isoform_domains, interactor_domains)
	max_pot_iso = get_max_potentials(stat_pot_data, isoform_domains)

	new_features = np.concatenate([np.squeeze(interact_pot_main_iso), np.squeeze(max_pot_main_iso), 
					np.squeeze(interact_pot_main_iso - interact_pot_iso), np.squeeze(max_pot_main_iso - max_pot_iso),
					[length_delta_ratio, length_delta, domains_lost, domains_procured, n_termini, c_termini, 
                         		max_gap_size, mean_gap_size, len(gap_sizes), large_change_num, small_change_num, domain_or_linker]])
	return new_features



def extract_features(main_sequence, interactor_sequence, isoform_sequence, config):
	#print 'Extract features'
	biochemical_main = extract_biochemical_features(main_sequence, config)
	biochemical_interactor = extract_biochemical_features(interactor_sequence, config)
	biochemical_isoform	= extract_biochemical_features(isoform_sequence, config)
	new_features = extract_new_features(main_iso_seq, interactor_seq, isoform_seq, config)
	#print biochemical_main


	#print np.shape(biochemical_main), np.shape(biochemical_interactor), np.shape(biochemical_isoform), np.shape(sequence_features)

	features = np.concatenate((biochemical_main, biochemical_interactor,
								biochemical_main - biochemical_isoform,
								new_features), axis=0)

	return features

if __name__ == "__main__":
	#main_iso = sys.argv[1]
	#interactor = sys.argv[2]
	#isoform = sys.argv[3]
	config = Config('config')
	main_iso_seq = sys.argv[1]
	interactor_seq = sys.argv[2]
	isoform_seq = sys.argv[3]
	#os.chdir('..')
	try:
		os.system('rm -r temp')
	except:
		pass
	try:
		os.system('rm tmp/*')
	except:
		pass
	features = extract_features(main_iso_seq, interactor_seq, isoform_seq, config)
	os.system('rm -r temp')
	os.system('rm tmp/*')
	pickle.dump(features, open('tmp/features.pickle', 'w'))
	#print features	
