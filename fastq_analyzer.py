import fileinput
import operator
import optparse
import sys
from collections import Counter
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import matplotlib.pyplot as plt
import numpy as np
import heapq

Encodings = {
	'Sanger': (33, 73),
    'Illumina-1.8': (33, 74),
    'Solexa': (59, 104),
    'Illumina-1.3': (64, 104),
	'Illumina-1.5': (66, 105)
}

def get_qual_range(quality_str):
	qualities = Counter(ord(char) for char in quality_str)
	min_q = min(qualities.keys())
	max_q = max(qualities.keys())
	return (min_q, max_q)


def get_encodings_in_range(rmin, rmax, ranges=Encodings):
	valid_encodings = []
	for encoding, (emin, emax) in ranges.items():
		if rmin >= emin and rmax <= emax:
			valid_encodings.append(encoding)
	return valid_encodings

def check_encoding(filepath):
	input = fileinput.input(filepath, openhook=fileinput.hook_compressed)
	gmin = 127
	gmax = 0
	encodings = []
	
	for i, line in enumerate(input):
		if (i+1) % 4 != 0:
			continue
		lmin, lmax = get_qual_range(line.rstrip())
		if lmin < gmin or lmax > gmax:
			gmin, gmax = min(lmin, gmin), max(lmax, gmax)
			encodings = get_encodings_in_range(gmin, gmax)
	input.close()
	print("{} : [{}, {}]".format(", ".join(encodings), gmin, gmax))

def analyze(filepath):
	record_iterator = SeqIO.parse(filepath, "fastq-sanger")
	i = 0
	data = []
	records = []
	for rec in record_iterator:
		i += 1
		records.append(rec)
		cg_count = rec.seq.count('C') + rec.seq.count('G')
		data.append((cg_count / len(rec.seq))*100)
	return records, data

def get_peaks(data):
	counter = Counter(data)
	points = [322, 417, 235]
	for key, value in counter.items():
		if value in points:
			print(value)
			print(key)

def get_reads(peaks, records, data, max_numb=5):
	reads = []
	for peak in peaks:
		index = 0
		numb = 0
		for d in data:
			if peak == d:
				reads.append(records[index])
				numb += 1
			if numb == max_numb:
				break
			index += 1
	return reads

def plot_data(data):
	counter = Counter(data)
	plt.plot(counter.keys(), counter.values(),'o')
	plt.xlabel("C/G dalis sekoje (%)", fontsize=15)
	plt.ylabel("Read'ų skačius", fontsize=15)
	plt.axis([0, 100, 0, 450])
	plt.xticks(np.arange(0, 100, 20))
	plt.yticks(np.arange(0, 450, 100))
	plt.tick_params(axis='both', which='major', labelsize=15)
	plt.grid(True, axis='y')
	plt.show()

def blast_search(reads):
	data = {}
	for read in reads:
		result_handle = NCBIWWW.qblast("blastn", "nt", read.seq, entrez_query='txid2[ORGN]', descriptions=1, alignments=1,hitlist_size=1)
		blast_record = NCBIXML.read(result_handle)
		for alignment in blast_record.alignments:
			org = alignment.title.split("|")[4]
			data[read.id] = org
			print("{} : {}", read.id, org)
	return data

if __name__ == "__main__":
	if(len(sys.argv) != 2):
		print("Please pass filepath as an argument")

	filepath = sys.argv[1]
	check_encoding(filepath)
	records, data = analyze(filepath)
	plot_data(data)

	# After analyzing the plot, peaks were chosen:
	peaks=[35.76158940397351, 53.64238410596026, 70.19867549668875]
	reads = get_reads(peaks, records, data)

	print("Blast search starting...")
	species = blast_search(reads)