import itertools
import regex as re
import pandas as pd
import numpy as np
from numpy import empty, array
np.set_printoptions(threshold=np.inf)
import csv
import pyfastx
from multiprocessing import Pool 
import ray
import sys
import argparse

parser = argparse.ArgumentParser(description='File input information')
parser.add_argument('--fasta', dest='fasta', type=str, help='Name of the input fasta file')
parser.add_argument('--repr', dest='representation', type=str, help='Choosen numeric representation')
parser.add_argument('--thread', dest='thread', type=int, help='Number of threads to use')
args = parser.parse_args()

print("File choosed:")
print(args.fasta)
print("Numeric Encoding:")
print(args.representation)
print("Number of threads to use")
print(args.thread)

def muri_fasta(sequence): #used by all numeric representations
	aux = pyfastx.Fasta(sequence)
	par = {}
	for ft in aux:
		par[ft.name] = str(ft)
	return(par)

def kmer_init(n): #used by natural vector
	dna = ["A", "C", "G", "T"]
	if(n == 1):
		return(dna)
	else:
		n = n-1
		return(["".join(x) for x in list(itertools.product(dna, kmer_init(n)))])

def kmer_order(n): #used by natural vector
	if(n==1):
		return(kmer_init(n))
	else:
		n=n-1
		return(kmer_order(n)+kmer_init(n+1))

@ray.remote
def natural(search_for, SEQ):  #used by natural and fast vector
	m=2
	pos = [(m.start()+1) for m in re.finditer(search_for, SEQ, overlapped=True)]
	ni = len(pos)
	if(ni == 0):
		return {"ni_"+search_for:0, "mu_"+search_for:0 , "D_2_"+search_for:0}
	else:
		mu = sum(pos)/ni
		return {"ni_"+search_for:ni, "mu_"+search_for:mu, "D_2_"+search_for: sum([((v-mu)**m) for v in pos]) / ((ni ** (m-1)) * (len(SEQ) ** (m-1))) }

def natural_count(SEQ):   #used by natural and fast vector
	f = {}
	aux = [natural.remote(i, SEQ) for i in ks]
	for elem in ray.get(aux):
		f = {**f, **elem}
	return f

def num_subsequences(seq, sub):    #used by magnus representation
	m, n = len(seq), len(sub)
	table = array([0] * n)
	for i in range(m):
		previous = 1
		for j in range(n):
			current = table[j]
			if seq[i] == sub[j]:
				table[j] += previous
			previous = current
	return table[n-1] if n else 1

def numberToBase(n, b):   #used by magnus representation
	if n == 0:
		return '0'
	digits = ''
	while n:
		digits+=str(int(n % b))
		n //= b
	return digits[::-1]

def mainmagnus(arg):    #used by magnus representation
	head, dna, index = arg
	dna = dna.upper()
	dnalength = len(dna)
	N=5
	mvl = int( (4**(N+1)-4)/3 ) 
	numwin=dnalength//N
	winarray = empty([numwin,mvl])
	for globalcounter in range(0,numwin):
		startindex=globalcounter*N
		endindex=startindex+N
		dnabase4 = ''
		magnusvec = array([0] * mvl)
		for i in range(startindex,endindex):
			if (dna[i]=='A'):
				dnabase4+='0'
			elif (dna[i]=='C'):
				dnabase4+='1'
			elif(dna[i]=='G'):
					dnabase4+='2'
			elif(dna[i]=='T' or dna[i]=='U'):
				dnabase4+='3'
			else:
				#print('DNA Sequence contains unallowed letters.')
				break
		counter = 0
		maxds = 4**N
		while counter<mvl:
			for ds in range(0,maxds):
				s = numberToBase(ds,4)
				while N-len(s)>=0:
					ps = int((4**len(s)-4)/3 +ds + 1 )
					alphas = num_subsequences(dnabase4,s)
					magnusvec[ps-1]=alphas
					s = '0'+s
					counter+=1
					
		winarray[globalcounter]=magnusvec
		
	summagnus=array([0]*mvl)
	for i in range(0,len(winarray)):
		summagnus=summagnus+winarray[i]
	
	with open(args.fasta+".magnus.csv", 'a', newline='') as file:
		writer = csv.writer(file)
		writer.writerow([head] + list(summagnus/numwin))

@ray.remote
def countk(search_for, SEQ):  #used by natural and fast vector
	pos = [(m.start()+1) for m in re.finditer(search_for, SEQ, overlapped=False)]
	return {search_for:(len(pos)/len(SEQ))}

def triplet_count(SEQ):   #used by natural and fast vector
	f = {}
	aux = [countk.remote(i, SEQ) for i in ks]
	for elem in ray.get(aux):
		f = {**f, **elem}
	return f

seq=""
db = muri_fasta(args.fasta)

if(args.representation == "Magnus"):
	with Pool(args.thread) as p:
		p.map(mainmagnus, iter([(h, s, i) for h, s, i in zip(db.keys(), db.values(), range(len(db)))]))
else:
	ray.init(num_cpus=args.thread)
	i = 1
	if(args.representation == "Fast"):
		ks = ["[AG]", "[CT]", "[AC]", "[GT]", "[GC]", "[AT]"]
		df = pd.DataFrame(columns = ["id", "ni_[AG]", "mu_[AG]", "D_2_[AG]", "ni_[CT]", "mu_[CT]", "D_2_[CT]", "ni_[AC]", "mu_[AC]", "D_2_[AC]", "ni_[GT]", "mu_[GT]", "D_2_[GT]", "ni_[GC]", "mu_[GC]", "D_2_[GC]", "ni_[AT]", "mu_[AT]", "D_2_[AT]"])
		df.to_csv("numeric("+args.fasta+").csv", index = True, header = True)
		for h, s in zip(db.keys(), db.values()):
			seq = s.upper()
			new = pd.DataFrame({'id':h,**natural_count(seq)}, index=[i])
			new.to_csv(args.fasta+".fast.csv", mode="a", index = True, header = False)
			i = i + 1
	else:
		if(args.representation == "4-mer"):
			ks = kmer_init(4)
			df = pd.DataFrame(columns = ["id"]+[k+x for x in ks for k in ["ni_", "mu_", "D_2_"]])
			df.to_csv(args.fasta+".4natvec.csv", index = True, header = True)
			for h, s in zip(db.keys(), db.values()):
				seq = s.upper()
				new = pd.DataFrame({'id':h,**natural_count(seq)}, index=[i])
				new.to_csv(args.fasta+".4natvec.csv", mode="a", index=True, header = False)
				i = i + 1
		else:
			if(args.representation == "6-mer"):
				ks = kmer_init(6)
				df = pd.DataFrame(columns = ["id"]+[k+x for x in ks for k in ["ni_", "mu_", "D_2_"]])
				df.to_csv(args.fasta+".6natvec.csv", index = True, header = True)
				for h, s in zip(db.keys(), db.values()):
					seq = s.upper()
					new = pd.DataFrame({'id':h,**natural_count(seq)}, index=[i])
					new.to_csv(args.fasta+".6natvec.csv", mode="a", index = True, header = False)
					i = i + 1
			else:
				if(args.representation == "c-4-mer"):
					ks = kmer_order(4)
					df = pd.DataFrame(columns = ["id"]+[k+x for x in ks for k in ["ni_", "mu_", "D_2_"]])
					df.to_csv(args.fasta+".c4natvec.csv", index = True, header = True)
					for h, s in zip(db.keys(), db.values()):
						seq = s.upper()
						new = pd.DataFrame({'id':h,**natural_count(seq)}, index=[i])
						new.to_csv(args.fasta+".c4natvec.csv", mode="a", index=True, header = False)
						i = i + 1
				else:
					if(args.representation == "c-6-mer"):
						ks = kmer_order(6)
						df = pd.DataFrame(columns = ["id"]+[k+x for x in ks for k in ["ni_", "mu_", "D_2_"]])
						df.to_csv(args.fasta+".c6natvec.csv", index = True, header = True)
						for h, s in zip(db.keys(), db.values()):
							seq = s.upper()
							new = pd.DataFrame({'id':h,**natural_count(seq)}, index=[i])
							new.to_csv(args.fasta+".c6natvec.csv", mode="a", index = True, header = False)
							i = i + 1
					else:
						ks = kmer_init(3)
						df = pd.DataFrame(columns = ["id"]+[k for k in ks])
						df.to_csv(args.fasta+".triplet.csv", index = True, header = True)
						for h, s in zip(db.keys(), db.values()):
							seq = s.upper()
							new = pd.DataFrame({'id':h,**triplet_count(seq)}, index=[i])
							new.to_csv(args.fasta+".triplet.csv", mode="a", index = True, header = False)
							i = i + 1
	ray.shutdown()

