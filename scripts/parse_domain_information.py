import pickle
import sys


def parse(files):
	parsed = {}
	binding = {}
	for fname in files:
		f = open(fname)
		for line in f:
			splitted = line.split('\t')
			split_size = len(splitted)
			id = splitted[0]
			domain = splitted[4]
			if id in parsed.keys():
				parsed[id].append(domain)
			else:
				parsed[id] = [domain]
			binding_sites = [splitted[6], splitted[7]]
			if id in binding.keys():
				binding[id][domain] = binding_sites
			else:
				binding[id] = {}
				binding[id][domain] = binding_sites
	pickle.dump(parsed, open('tmp/parsed.pickle', 'w'))
	pickle.dump(binding, open('tmp/binding_sites.pickle', 'w'))

if __name__ == "__main__":
	files = []
	if len(sys.argv) > 2:
		files =  sys.argv[1:]
	else:
		files = [sys.argv[1]]
	parse(files)


