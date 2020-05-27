#prepare mysnps.vcf
from __future__ import division
import re,sys,os,time
import argparse
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
parser = argparse.ArgumentParser()
parser.add_argument("--vcf",type=str,required=True, help="Input SNP files for motif analysis<xxx.vcf format>")
parser.add_argument("--fa",type=str,required=True, help="All-chromosome.fa file for sequencing extracting")
parser.add_argument("--db",type=str,required=True, help="Path for motif database <xxx/fnGWAS/data/meme.motif>")
parser.add_argument("--output",type=str,required=False, default="motif.output",help="Output dictory for motif analysis")
parser.add_argument("--core", default=4,type=int,required=False, help="Prallel Threads for motif analysis. default 4, max 20")
args = parser.parse_args()
start=time.time()

if args.output not in os.listdir(os.getcwd()):
	os.mkdir(args.output)
os.chdir(args.output)

#extract region
f1=open("../"+args.vcf)
w={}
for line in f1:
	lines=line.strip().split('\t')
	p=re.match(".*CAF=([0-9\.]*[0-9])[^0-9]*([0-9][0-9\.]*[0-9]).*",line)
	if lines[0] not in w:
		w[lines[0]]={}
	if p:
		w[lines[0]][lines[2]]=[lines[1],lines[3],str(1-float(p.group(1))),lines[4],p.group(1)]
	if not p:
		w[lines[0]][lines[2]]=[lines[1],lines[3],"NA",lines[4],"NA"]

#extract sequence
record=SeqIO.to_dict(SeqIO.parse(args.fa,"fasta"))
f2=open(re.sub("\.vcf","",args.vcf)+".fa",'w')
for p in w:
	for q in w.get(p):
		n=w.get(p).get(q)
		n1=str(record["chr"+str(p)].seq[(int(n[0])-31):int(n[0])-1])
		n2=str(record["chr"+str(p)].seq[int(n[0]):(int(n[0])+30)])
		w1={};n3=n[1].split(',');n4=n[3].split(',')
		for n5 in n3:
			w1[n5]=n[2]
		for n5 in n4:
			w1[n5]=n[4]
		for d in w1:
			f2.write('>'+q+'-'+d+'-'+w1.get(d)+'\n')
			f2.write(n1+d+n2+'\n')
	print "chr"+p+" is done"
f1.close()
f2.close()
print "Sequence extracted: 1/3" 

#motif analysis(Prallel)
#extract sequence info
f=open(re.sub("\.vcf",".fa",args.vcf))
w6={};linecount=0
for line in f:
	if line[0]==">":
		motif=line.strip()
	else:
		w6[motif]=line.strip()
	linecount+=1
f.close()

#adjust suitable Prallel core counts
if int(linecount)<=500:
	core=1
elif int(args.core)>=(int(linecount)/2):
	core=int(linecount/500)+1
else:
	core=int(args.core)
if core>=20:
	core=20

#split sequence file
motif=w6.keys()
w5={}
size=int((linecount/2)/core)+1
for i in range(0,core):
	w5[str(i)]={}
	for j in range(i*size,i*size+size):
		if j<len(motif):
			w5[str(i)][motif[j]]=""

#Read motif name annotation
#w9:{name:TF:URL}
f=open(args.db+"/fimo.motifs.anno")
w9={}
for line in f:
	lines=line.strip().split("\t")
	kk=[]
	for p in re.split("::|;",lines[1]):
		kk.append(p)
	w9[lines[0]]={}
	for p in kk:
		w9[lines[0]][p]=lines[2]
f.close()

#output_split_sequence
for p1 in w5:
	f=open("split."+str(p1),"w")
	for p2 in w5.get(p1):
		f.write(p2+"\n"+w6.get(p2)+"\n")
	f.close()
print "Sequence split finished: 2/3"

def task(n):
	db="split"+str(n)
	os.environ["db"]=db
	os.environ["sequence"]="split."+str(n)
	os.environ["motif"]=args.db+"/fimo.meme"
	os.system("fimo --thresh 1e-4 --oc $db $motif $sequence")
	f2=open("./"+db+'/fimo.txt')
	f3=open('fimo-'+db+'.txt.anno','w')
	f3.write('meme-db\tTF\t'+f2.readline().strip()[1:]+'\tURL\n')
	for line in f2:
		lines=line.strip().split('\t')
		if lines[0] in w9:
			if float(lines[2])<31 and float(lines[3])>31:
				f3.write(re.sub("[$].*$","",lines[0])+'\t'+";".join(w9.get(lines[0]).keys())+"\t"+line.strip()+"\t"+w9.get(lines[0]).values()[0]+'\n')
	f2.close()
	f3.close()
	os.system("rm -f -r $db")

if  __name__ == '__main__':
	pool=Pool(core)
	pool.map(task,w5.keys())
	pool.close()
	pool.join()

os.system("head -1 fimo-split0.txt.anno>Result.motif.fimo.txt")
os.system("awk '$1!~/meme-db/' fimo-split*.txt.anno>>Result.motif.fimo.txt")
os.system("rm -f split.* fimo-split*.txt.anno")
print "Motif analysis finished: 3/3"

#anno TF & make w1:rs.GC-db-TF
#w:motif.db
#w1:rs.GC-db-TF
w={}
w1={}
f2=open("Result.motif.fimo.txt")
f2.readline()
for line in f2:
	lines=line.strip().split('\t')
	snp=lines[3].split("-")[0]
	allele=lines[3].split("-")[1]
	database=lines[0]
	if snp not in w1:
		w1[snp]={}
	if database not in w1[snp]:
		w1[snp][database]={}
	if lines[1] not in w1[snp][database]:
		w1[snp][database][lines[1]]={}
	w1[snp][database][lines[1]][allele]=""
	w[database]=""
f2.close()

#output all results
f4=open("meme-suit-motif-all.txt","w")
m=w.keys()
f4.write("SNP\t"+"\t".join(m)+"\n")
for p in w1:
	mm=[]
	for pp in m:
		if pp not in w1[p]:
			mm.append(".")
		else:
			mmm=[]
			for k in w1.get(p).get(pp).keys():
				for kk in w1.get(p).get(pp).get(k).keys():
					mmm.append(k+"("+kk+")")
			mm.append(", ".join(mmm))
	f4.write(p+"\t"+"\t".join(mm)+"\n")

#filtered no-unique TFs
p1=w1.keys()
for pp1 in p1:
	p2=w1.get(pp1).keys()
	for pp2 in p2:
		p3=w1.get(pp1).get(pp2).keys()
		for pp3 in p3:
			if len(w1.get(pp1).get(pp2).get(pp3).keys())>=2:
				del w1[pp1][pp2][pp3]
		if not w1[pp1][pp2]:
			del w1[pp1][pp2]
	if not w1[pp1]:
		del w1[pp1]

#write unique results
f5=open("Allele_specific_meme-suit-motif.txt","w")
f5.write("SNP\t"+"\t".join(m)+"\n")
for p in w1:
	mm=[]
	for pp in m:
		if pp not in w1[p]:
			mm.append(".")
		else:
			mmm=[]
			for k in w1.get(p).get(pp).keys():
				for kk in w1.get(p).get(pp).get(k).keys():
					mmm.append(k+"("+kk+")")
			mm.append(", ".join(sorted(mmm)))
	f5.write(p+"\t"+"\t".join(mm)+"\n")
w1.clear()

#write no-unique TFs with 2 or more hists
f5=open("Allele_specific_meme-suit-motif.txt")
f6=open("Allele_specific_meme-suit-motif-2hits.txt","w")
f6.write(f5.readline())
for line in f5:
	lines=line.strip().split("\t")
	mm=[]
	for p1 in lines[1:]:
		if p1!=".":
			for pp1 in p1.split(", "):
				if "::" not in pp1:
					mm.append(pp1)
				else:
					for ppp1 in pp1.split("::"):
						mm.append(ppp1)
	mm1="";mm2=""
	for p1 in lines[1:]:
		if p1==".":
			mm1=mm1+"\t."
		else:
			mm3={}
			for pp1 in p1.split(", "):
				if "::" in pp1:
					for ppp1 in pp1.split("::"):
						if mm.count(ppp1)>=2:
							mm3[pp1]=""
				elif mm.count(pp1)>=2:
					mm3[pp1]=""
			if mm3:
				mm2="yes"
				mm1=mm1+"\t"+", ".join(sorted(mm3.keys()))
			if not mm3:
				mm1=mm1+"\t."
	if mm2:
		f6.write(lines[0]+mm1+"\n")
f5.close()
f6.close()

#annotate results
f1=open("Allele_specific_meme-suit-motif.txt")
w1={}
data=f1.readline().strip().split("\t")
for p in data[1:]:
	w1[p]={}
for line in f1:
	lines=line.strip().split("\t")
	for i in range(1,len(lines)):
		if lines[i]!=".":
			for motif in lines[i].split(", "):
				TF=motif.split("(")[0]
				A=motif.split("(")[1].split(")")[0]
				if lines[0]+"-"+A not in w1[data[i]]:
					w1[data[i]][lines[0]+"-"+A]={}
				w1[data[i]][lines[0]+"-"+A][TF]=""
f1.close()

f=open("Allele_specific_meme-suit-motif.txt.anno","w")
f.write("SNP\tDatabase\tTF\tbinding.sequence(SNP.position)\tscore\tp-value\tURL\n")
w4={}
f2=open("Result.motif.fimo.txt")
f2.readline()
for line in f2:
	lines=line.strip().split("\t")
	SNP=re.sub("-NA$","",lines[3])
	if SNP in w1[lines[0]]:
		if lines[1] in w1[lines[0]][SNP]:
			if SNP+"\t"+lines[0]+"\t"+lines[1] not in w4:
				if lines[6]=="+":
					w4[SNP+"\t"+lines[0]+"\t"+lines[1]]=lines[10]+"(+,"+str(32-int(lines[4]))+")\t"+lines[7]+"\t"+lines[8]+"\t"+lines[-1]
				elif lines[6]=="-":
					w4[SNP+"\t"+lines[0]+"\t"+lines[1]]=lines[10]+"(-,"+str(int(lines[5])-30)+")\t"+lines[7]+"\t"+lines[8]+"\t"+lines[-1]
				pvalue=float(lines[8])
			else:
				if float(lines[8])<pvalue:
					if lines[6]=="+":
						w4[SNP+"\t"+lines[0]+"\t"+lines[1]]=lines[10]+"(+,"+str(32-int(lines[4]))+")\t"+lines[7]+"\t"+lines[8]+"\t"+lines[-1]
					elif lines[6]=="-":
						w4[SNP+"\t"+lines[0]+"\t"+lines[1]]=lines[10]+"(-,"+str(int(lines[5])-30)+")\t"+lines[7]+"\t"+lines[8]+"\t"+lines[-1]
				pvalue=float(lines[8])
f2.close()

for p in w4:
	f.write(p+"\t"+w4.get(p)+"\n")
f.close()
print "Motif analysis finished"
end=time.time()
print '%s seconds' % (end - start)
