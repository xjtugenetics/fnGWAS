#Assign target gene for user-defined SNPs input
#prepare user-defined data in ./fnGWAS/data (see README)
#dowlaod annotate_variation.pl (http://annovar.openbioinformatics.org/en/latest/user-guide/download/) and copy it here
from __future__ import division
import re,sys,os,time
import argparse
from multiprocessing import Pool,Manager
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
parser = argparse.ArgumentParser()
parser.add_argument("--SNP",type=str,required=True,help="SNP input file<.avinput format>")
parser.add_argument("--db",type=str,required=True,help="Path for all annotation file<./fnGWAS/data>")
parser.add_argument("--output",type=str,required=False,default="Output_target_gene",help="Output dictory")
args = parser.parse_args()
start=time.time()

if "annotate_variation.pl" not in os.listdir(os.getcwd()):
	print "annotate_variation.pl not exit!!!"
	quit()

if args.output not in os.listdir(os.getcwd()):
		os.mkdir(args.output)
os.chdir(args.output)

#Adjust potential wrong ref in input avinput files
record=SeqIO.to_dict(SeqIO.parse(args.db+"/All-chromosome.fa","fasta"))
f1=open("../"+args.SNP)
f2=open(args.SNP+".adjust","w")
for line in f1:
	lines=line.strip().split("\t")
	ref=record["chr"+str(lines[0])].seq[int(lines[1])-1]
	if lines[3]==ref:
		f2.write(line)
	else:
		f2.write("\t".join(lines[0:3])+"\t"+lines[4]+"\t"+lines[3]+"\t"+lines[5]+"\n")
f1.close()
f2.close()
os.environ["snp"]=args.SNP
os.system("mv $snp.adjust $snp")

#Extracted coding/splcing SNPs
#w1{noncodingSNP: region}
os.environ["db"]=args.db
os.system("../annotate_variation.pl -geneanno -buildver hg19 $snp $db -dbtype wgEncodeGencodeBasicV19")
f1=open(args.SNP+".variant_function")
f2=open(args.SNP+".coding.promoter.gene","w")
f3=open(args.SNP+".noncoding.bed","w")
w1={}
for line in f1:
	lines=line.strip().split("\t")
	mm=lines[0].split(";")
	if "exonic" not in mm and "splicing" not in mm:
		f3.write("chr"+lines[2]+"\t"+str(int(lines[3])-1)+"\t"+str(int(lines[3])+1)+"\t"+lines[-1]+"\n")
		w1[lines[-1]]="chr"+"\t".join(lines[2:5])+"\t"+lines[-1]
	if ";" not in lines[0]:
		if lines[0]=="splicing":
			f2.write("chr"+"\t".join(lines[2:5])+"\t"+lines[-1]+"\t"+lines[0]+"\t"+lines[1].split("(")[0]+"\n")
	else:
		mm=lines[0].split(";")
		mm1=lines[1].split(";")
		for i in range(0,2):
			if mm[i]=="splicing":
				f2.write("chr"+"\t".join(lines[2:5])+"\t"+lines[-1]+"\t"+mm[i]+"\t"+mm1[i].split("(")[0]+"\n")
f1.close()
f3.close()
f1=open(args.SNP+".exonic_variant_function")
for line in f1:
	lines=line.strip().split("\t")
	f2.write("chr"+"\t".join(lines[3:6])+"\t"+lines[-1]+"\t"+lines[1]+"\t"+re.sub(":.*$","",lines[2])+"\n")
f1.close()

#w2:local pairs
w2={}
mm=os.popen("bedtools intersect -a $snp.noncoding.bed -b $db/gencode.v19.transcript_TSS1KB_promoter.bed -wo")
for line in mm:
	lines=line.strip().split("\t")
	f2.write("\t".join(lines[:4])+"\tpromoter\t"+lines[7]+"\n")
	if lines[3] not in w2:
		w2[lines[3]]={}
	w2[lines[3]][lines[7]]=""
mm.close()
f2.close()
end=time.time()
print "Local gene extracted" + "....."+'%s seconds' % (end - start)

#Anno cis.QTL
w3={};QTL=[]
data1=os.listdir(args.db+"/QTL/")
for mm1 in data1:
	f4=open(args.db+"/QTL/"+mm1)
	mm2=re.sub(".QTL","",mm1)
	mm3=mm2+"("+f4.readline().strip().split("\t")[-1]+")"
	QTL.append(mm3)
	for line in f4:
		lines=line.strip().split("\t")
		if lines[0] in w1:
			if lines[0] not in w3:
				w3[lines[0]]={}
			if lines[1] not in w3[lines[0]]:
				w3[lines[0]][lines[1]]={}
			w3[lines[0]][lines[1]][mm3]=lines[-1]
	f4.close()
end=time.time()
print "cis-QTL annotation finished"+"....."+'%s seconds' % (end - start)

#extract local QTL.pairs
f2=open(args.SNP+".noncoing.local.target.genes","w")
f2.write("SNP\tGene\t"+"\t".join(data1)+"\n")
for p1 in w3:
	if p1 in w2:
		for p2 in w3[p1]:
			if p2 in w2[p1]:
				mm=[]
				for p3 in QTL:
					mm.append(w3.get(p1).get(p2).get(p3,"-"))
				f2.write(p1+"\t"+p2+"\t"+"\t".join(mm)+"\n")
f2.close()
end=time.time()
print "Local gene extracted"+"....."+'%s seconds' % (end - start)

##Anno distal target gene
#filter SNPs
f=open("QTL.pairs","w")
for p1 in w3:
	for p2 in w3[p1]:
		f.write(p1+";"+p2+"\n")
f.close()
#extract SNP.gene within 1M
os.system("awk 'BEGIN{OFS=\"\t\"}{if($2<=999999)print $1,\"1\",$3+999999,$4;if($2>999999)print $1,$2-999999,$3+999999,$4}' $snp.noncoding.bed|bedtools intersect -a - -b $db/gencode.v19.transcript_TSS1KB_promoter.bed -wo|awk 'BEGIN{OFS=\"\t\"}{print $1,$3-999999-1,$3-999999+1,$5,$6,$7,$4\";\"$8}'|grep -w -Ff QTL.pairs ->$snp.1M.gene.bedpe")

#HiC annotation
w4={};HiC=[]
data2=os.listdir(args.db+"/HiC/")
for mm1 in data2:
	os.environ["mm1"]=mm1
	lines2=os.popen("bedtools pairtopair -a $snp.1M.gene.bedpe -b $db/HiC/$mm1 -type both")
	HiC.append(re.sub(".bedpe","",mm1))
	for line in lines2:
		lines=line.strip().split("\t")
		mm2=lines[6].split(";")
		if mm2[1] not in w2.get(mm2[0],{}):
			if mm2[0] not in w4:
				w4[mm2[0]]={}
			if mm2[1] not in w4[mm2[0]]:
				w4[mm2[0]][mm2[1]]={}
			pair=lines[-1]+":"+lines[7]+":"+lines[8]+"-"+lines[9]+","+lines[10]+":"+lines[11]+"-"+lines[12]
			data=re.sub(".bedpe","",mm1)
			if data not in w4[mm2[0]][mm2[1]]:
				w4[mm2[0]][mm2[1]][data]={}
			w4[mm2[0]][mm2[1]][data][pair]=""
	lines2.close()
end=time.time()
print "HiC annotation finished" + "....."+'%s seconds' % (end - start)

#Anno Distal.pairs with both QTL&HiC
f3=open(args.SNP+".noncoing.distal.target.genes","w")
f3.write("SNP\tGene\t"+"\t".join(HiC)+"\t"+"\t".join(QTL)+"\n")
for p1 in w4:
	if p1 in w3:
		for p2 in w4[p1]:
			if p2 in w3[p1]:
				mm=[]
				for p3 in HiC:
					if p3 in w4[p1][p2]:
						mm.append(";".join(w4[p1][p2][p3].keys()))
					elif p3 not in w4[p1][p2]:
						mm.append("-")
				for p4 in QTL:
					mm.append(w3.get(p1).get(p2).get(p4,"-"))
				f3.write(p1+"\t"+p2+"\t"+"\t".join(mm)+"\n")
f3.close()
end=time.time()
print "Distal gene extracted"+"....."+'%s seconds' % (end - start)
