#Rscript picolo_pipeline.R Top.GWAS.QTL.Pvalue.pairs
library("piccolo")
argv <- commandArgs(TRUE)
data=read.table(argv[1],header=T,colClasses="character")
data$"H3"="NA"
data$"H4"="NA"
for(i in c(1:nrow(data))){
	res1 <- try({m1=pics.download(rsid=data[i,2], pvalue=-log10(as.numeric(data[i,3])))},silent = T)
	if (inherits(res1, "try-error")){
		print(paste(data[i,2],"missing",sep=" "))
	}
	res2 <- try({m2=pics.download(rsid=data[i,5], pvalue=-log10(as.numeric(data[i,6])))},silent = T)
	if (inherits(res2, "try-error")){
		print(paste(data[i,5],"missing",sep=" "))
	}
	if (ncol(res1)==6 & ncol(res2)==6){
		res3=merge(res1,res2,by="Linked_SNP")
		if(nrow(res3)>1){
			result=as.character(pics.coloc.lite(res1,res2))
			data[i,5]=result[1]
			data[i,6]=result[2]
			print(paste(data[i,1],data[i,3],"..pairs..over",sep="..."))
		}
		if(nrow(res3)<=1){
			print(paste(data[i,1],data[i,3],"..pairs..missing",sep="..."))
		}
	}
}
write.table(data,file=paste(argv[1],".piccolo.result",sep=""),quote=F,row.names=F,col.names=T)
