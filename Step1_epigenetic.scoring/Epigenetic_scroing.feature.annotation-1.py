#Annotate SNP and output feature enrichment analysis result
from __future__ import division
import numpy as np
from scipy.stats import chi2_contingency
from multiprocessing import Pool
import re,sys,os,time
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--snp", help="GWAS SNP input file <xxx.bed format>",required=True,type=str)
parser.add_argument("--db", help="Path for database storing reference genotype data for LD analysis <./data/1000g.V3.EUR.genotype>",required=True,type=str)
parser.add_argument("--anno", help="file stored feature path <xxx.txt>",required=True,type=str)
parser.add_argument("--core", help="Prallel Threads <default:5,max=15>",required=False,default=5,type=int)
parser.add_argument("--output", help="Output dictory name",required=False,default="output.SNP.scroing",type=str)
args = parser.parse_args()
start=time.time()

if args.output not in os.listdir(os.getcwd()):
	os.mkdir(args.output)
os.chdir(args.output)

#extract GWAS merged loci
os.environ["snp"]=args.snp
os.environ["MHC"]=re.sub("/[^/]*$","",args.db)+"/MHC_NCBI_region.bed"
print re.sub("/[^\/]*$","",args.db)
os.system("awk 'BEGIN{OFS=\"\t\"}{print $1,($2-999999),($3+999999),$4}' ../$snp|awk 'BEGIN{OFS=\"\t\"}{if($2<=0)print $1,\"1\",$3,$4;if($2>0)print $0}'|sort -k1,1 -k2,2n>$snp\_1000k.bed")
os.system("bedtools merge -i $snp\_1000k.bed -c 4 -o collapse -delim \"_\"|bedtools intersect -v -a - -b $MHC>$snp\_1000k.merged.noHLA.bed")
os.system("bedtools merge -i $snp\_1000k.bed -c 4 -o collapse -delim \"_\"|bedtools intersect -a - -b $MHC -wa|bedtools subtract -a - -b $MHC>>$snp\_1000k.merged.noHLA.bed")
os.system("rm -f $snp\_1000k.bed")

f=open(args.snp+"_1000k.merged.noHLA.bed")
m=[]
for line in f:
	m.append(line)
f.close()
def run(m):
	lines=m.strip().split("\t")
	os.environ["chr"]=re.sub("chr","",lines[0])
	os.environ["start"]=lines[1]
	os.environ["end"]=lines[2]
	f=open(lines[0]+"_"+lines[1]+"_"+lines[2]+".snps","w")
	for p in lines[-1].split("_"):
		f.write(p+"\n")
	f.close()
	os.environ["snps"]=lines[0]+"_"+lines[1]+"_"+lines[2]+".snps"
	os.environ["genotype"]=args.db+"/chr"+str(re.sub("chr","",lines[0]))
	os.system("plink --bfile $genotype --chr $chr --from-bp $start --to-bp $end --make-bed --out $snps")
	os.system("plink --bfile $snps --r2 --ld-snp-list $snps --ld-window-kb 1000 --ld-window-r2 0.8 --out $snps")
	os.system("awk 'BEGIN{FS=\"\t\";OFS=\"\t\"}{print \"chr\"$1,$4-1,$4+1,$2}' $snps.bim>$snps\.bed.background")
	os.system("awk 'BEGIN{OFS=\"\t\"}{if(NR>1)print \"chr\"$4,$5-1,$5+1,$6}' $snps.ld|sort|uniq>$snps\.bed.positive")
	os.system("rm -f $snps.nosex $snps.bim $snps.bed $snps.fam $snps.log $snps.ld $snps")

if  __name__ == '__main__':
	pool=Pool(args.core)
	pool.map(run,m)
	pool.close()
	pool.join()
os.environ["snp"]=args.snp
os.system("cat chr*.background>$snp\.background.SNPs.bed")
os.system("cat chr*.positive>$snp\.positive.SNPs.bed")
os.system("rm -f chr*.background chr*.positive")

##annotate positive SNPs
print "Annotate positive SNPs"
SNP=args.snp+".positive.SNPs.bed"
#read bed file
f1=open(SNP)
linecount=0
w1={}
for line in f1:
	w1[line.strip()]=""
	linecount+=1
f1.close()

#adjust suitable Prallel core counts
if int(linecount)<=100:
	core=1
elif int(args.core)>=(int(linecount)/2):
	core=int(linecount/100)+1
else:
	core=int(args.core)
if core>=15:
	core=15

#split iuput bed file
snp=w1.keys()
w2={}
size=int(linecount/core)+1
for i in range(0,core):
	w2[str(i)]={}
	for j in range(i*size,i*size+size):
		if j<len(snp):
			w2[str(i)][snp[j]]=""

#output split.bed file
for p1 in w2:
	f=open(SNP+"."+args.anno+"."+str(p1),"w")
	for p2 in w2.get(p1):
		f.write(p2+"\n")
	f.close()
print "SNP file splited: 1/3"

#read feature list
f4=open("../"+args.anno)
w3={}
for line in f4:
	w3[line.strip()]=""
f4.close()

#def task
def task(i):
	os.environ["snp"]=SNP+"."+args.anno+"."+str(i)
	for m in w3:
		type=m.split("/")[-1].split("_")[0]
		cell=re.sub("^[^_][^_]*_|.bed","",m.split("/")[-1])
		features=m.split("/")[-1]
		os.environ["feature"]=m
		output=os.popen("bedtools intersect -a $feature -b $snp -wo")
		w4={};w5={}
		for line in output:
			lines=line.strip().split("\t")
			if lines[3] not in w4:
				w4[lines[3]]=0
			w4[lines[3]]=w4[lines[3]]+1
			if lines[7] not in w5:
				w5[lines[7]]={}
			w5[lines[7]][lines[3]]=""
		f5=open(SNP+"."+args.anno+"."+str(i)+"."+features+".count.temp","w")
		for p in w4:
			f5.write(type+"\t"+cell+"\t"+p+"\t"+str(w4.get(p))+"\n")
		f5.close()
		f6=open(SNP+"."+args.anno+"."+str(i)+"."+features+".SNP.anno.temp","w")
		for p in w5:
			f6.write(p+"\t"+type+"\t"+cell+"\t"+";".join(w5.get(p).keys())+"\n")
		f6.close()
	os.system("rm -f $snp")

if  __name__ == '__main__':
	pool=Pool(core)
	pool.map(task,w2.keys())
	pool.close()
	pool.join()
print "Annotated 2/3"

#merge result
os.environ["count"]=SNP+"."+args.anno
os.system("cat $count.*.count.temp>$count.count")
os.system("cat $count.*.SNP.anno.temp>$count\.SNP.anno")
os.system("rm -f $count.*.count.temp $count.*.SNP.anno.temp")
f6=open(SNP+"."+args.anno+".count")
w4={}
for line in f6:
	lines=line.strip().split("\t")
	if lines[0]+"\t"+lines[1]+"\t"+lines[2] not in w4:
		w4[lines[0]+"\t"+lines[1]+"\t"+lines[2]]=0
	w4[lines[0]+"\t"+lines[1]+"\t"+lines[2]]=w4[lines[0]+"\t"+lines[1]+"\t"+lines[2]]+int(lines[-1])
f6.close()
f7=open(SNP+"."+args.anno+".count","w")
for p in w4:
	f7.write(p+"\t"+str(w4.get(p))+"\n")
f7.close()
print "Results were merged 3/3!!!"

##annotate background SNPs
print "Annotate background SNPs"
SNP=args.snp+".background.SNPs.bed"
#read bed file
f1=open(SNP)
linecount=0
w1={}
for line in f1:
	w1[line.strip()]=""
	linecount+=1
f1.close()

#adjust suitable Prallel core counts
if int(linecount)<=100:
	core=1
elif int(args.core)>=(int(linecount)/2):
	core=int(linecount/100)+1
else:
	core=int(args.core)
if core>=15:
	core=15

#split iuput bed file
snp=w1.keys()
w2={}
size=int(linecount/core)+1
for i in range(0,core):
	w2[str(i)]={}
	for j in range(i*size,i*size+size):
		if j<len(snp):
			w2[str(i)][snp[j]]=""

#output split.bed file
for p1 in w2:
	f=open(SNP+"."+args.anno+"."+str(p1),"w")
	for p2 in w2.get(p1):
		f.write(p2+"\n")
	f.close()
print "SNP file splited: 1/3"

#read feature list
f4=open("../"+args.anno)
w3={}
for line in f4:
	w3[line.strip()]=""
f4.close()

#def task
def task(i):
	os.environ["snp"]=SNP+"."+args.anno+"."+str(i)
	for m in w3:
		type=m.split("/")[-1].split("_")[0]
		cell=re.sub("^[^_][^_]*_|.bed","",m.split("/")[-1])
		features=m.split("/")[-1]
		os.environ["feature"]=m
		output=os.popen("bedtools intersect -a $feature -b $snp -c")
		w4={}
		for line in output:
			lines=line.strip().split("\t")
			if lines[-2] not in w4:
				w4[lines[3]]=0
			w4[lines[-2]]=w4[lines[-2]]+int(lines[-1])
		f5=open(SNP+"."+args.anno+"."+str(i)+"."+features+".count.temp","w")
		for p in w4:
			f5.write(type+"\t"+cell+"\t"+p+"\t"+str(w4.get(p))+"\n")
		f5.close()
	os.system("rm -f $snp")

if  __name__ == '__main__':
	pool=Pool(core)
	pool.map(task,w2.keys())
	pool.close()
	pool.join()
print "Annotated 2/3"

#merge result
os.environ["count"]=SNP+"."+args.anno
os.system("cat $count.*.count.temp>$count.count")
os.system("rm -f $count.*.count.temp")
f6=open(SNP+"."+args.anno+".count")
w6={}
for line in f6:
	lines=line.strip().split("\t")
	if lines[0]+"\t"+lines[1]+"\t"+lines[2] not in w6:
		w6[lines[0]+"\t"+lines[1]+"\t"+lines[2]]=0
	w6[lines[0]+"\t"+lines[1]+"\t"+lines[2]]=w6[lines[0]+"\t"+lines[1]+"\t"+lines[2]]+int(lines[-1])
f6.close()
f7=open(SNP+"."+args.anno+".count","w")
for p in w6:
	f7.write(p+"\t"+str(w6.get(p))+"\n")
f7.close()
print "Results were merged 3/3!!!"

#Feature.enrichment.analysis
os.environ["N1"]=args.snp+".positive.SNPs.bed"
os.environ["N2"]=args.snp+".background.SNPs.bed"
m1=os.popen("wc -l $N1")
for line in m1:
	N11=int(re.split("[ \t][ \t]*",line.strip())[0])
m1.close()
m1=os.popen("wc -l $N2")
for line in m1:
	N22=int(re.split("[ \t][ \t]*",line.strip())[0])
m1.close()
w7={}
for p in w4:
	w7[p]=""
for p in w6:
	w7[p]=""
f=open(args.snp+".fearture.anno.summry","w")
f.write("Type\tCell\tfeature\tpositive.SNPs\tpositive.total\tbackground.SNPs\tbackground.total\tFC\tPvalue\n")

#output feature.FC & Pvalue
for p in w7:
	n1=w4.get(p,0)
	n2=w6.get(p,0)
	data=np.array([[n1,N11-n1],[n2,N22-n2]])
	Pvalue=chi2_contingency(data)[1]
	FC=str((n1*N22)/(N11*n2))
	f.write(p+"\t"+str(n1)+"//"+str(N11)+"\t"+str(n2)+"//"+str(N22)+"\t"+FC+"\t"+str(Pvalue)+"\n")
f.close()

print "Feature enrichment analysis finished"
end=time.time()
print '%s seconds' % (end - start)
