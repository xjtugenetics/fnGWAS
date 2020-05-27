#Rscript SNP2Gene_coloc.analysis.R SNP_gene.pairs
library("coloc")
argv <- commandArgs(TRUE)
base=getwd()

data=read.table(argv[1],header=T)
data$"overlapping.SNPs_count"=0
data$"PP.H0.abf"="NA"
data$"PP.H1.abf"="NA"
data$"PP.H2.abf"="NA"
data$"PP.H3.abf"="NA"
data$"PP.H4.abf"="NA"
for (i in seq(1,nrow(data))){
	gwas<-read.table(as.character(data[i,5]),header=T)
	QTL<-read.table(as.character(data[i,6]),header=T)
	gwas=gwas[gwas$"p">0 & gwas$"p"<1 & gwas$"MAF">0 & gwas$"MAF"<1,]
	QTL=QTL[QTL$"p">0 & QTL$"p"<1 & QTL$"MAF">0 & QTL$"MAF"<1,]
	gwasQTL=merge(gwas,QTL,by.x="SNP",by.y="SNP")
	if(nrow(gwasQTL)>10){
		H<-coloc.abf(dataset1=list(snp=gwas$"SNP",N=gwas$"N", pvalues=gwas$"p", type="quant", MAF=gwas$"MAF"),dataset2=list(snp=QTL$"SNP",N=QTL$"N", pvalues=QTL$"p", type="quant", MAF=QTL$"MAF"))
		m=as.character(H[[1]])
		data[i,7]=m[1]
		data[i,8]=m[2]
		data[i,9]=m[3]
		data[i,10]=m[4]
		data[i,11]=m[5]
		data[i,12]=m[6]
	}
	rm(gwas);rm(QTL);rm(gwasQTL);gc()
}
print("caculation finihsed!!")

data=data[,-c(5,6)]
write.table(data,file=paste(argv[1],"coloc.txt",sep="."),sep="\t",quote=F,row.names=F,col.names=T)
