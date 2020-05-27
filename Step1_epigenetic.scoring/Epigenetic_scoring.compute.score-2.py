#Select features first!
#output SNP scores
from __future__ import division
from multiprocessing import Pool
import re,sys,os,time
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--snp", help="SNP annotation results <xxx.SNP.anno>",required=True,type=str)
parser.add_argument("--feature", help="User-defined significant feature list <xxx.FC>",required=True,type=str)
parser.add_argument("--output", help="Output dictory name",required=False,default="output.SNP.scroing",type=str)
args = parser.parse_args()
start=time.time()

if args.output not in os.listdir(os.getcwd()):
		os.mkdir(args.output)
os.chdir(args.output)
f1=open("../"+args.feature)
w1={}
for line in f1:
	lines=line.strip().split("\t")
	if lines[0] not in w1:
		w1[lines[0]]={}
	w1[lines[0]][lines[1]+";"+lines[2]]=float(lines[3])
	print lines
f1.close()

f2=open("../"+args.snp)
f3=open(args.snp+".score","w")
type=w1.keys()
f3.write("SNP\t"+"\t".join(type)+"\n")
w2={}
for line in f2:
	lines=line.strip().split("\t")
	if lines[1] in w1:
		for p in lines[3].split(";"):
			if lines[2]+";"+p in w1[lines[1]]:
				if lines[0] not in w2:
					w2[lines[0]]={}
				if lines[1] not in w2[lines[0]]:
					w2[lines[0]][lines[1]]=0
				w2[lines[0]][lines[1]]=w2[lines[0]][lines[1]]+w1[lines[1]][lines[2]+";"+p]
f2.close()
for p in w2:
	mm=[]
	for p1 in type:
		mm.append(str(w2.get(p).get(p1,"0")))
	f3.write(p+"\t"+"\t".join(mm)+"\n")
f3.close()
print "Scoring finished"
end=time.time()
print '%s seconds' % (end - start)
