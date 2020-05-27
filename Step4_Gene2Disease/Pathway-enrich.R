#Rscript Pathway-gene-list.R gene-list
#http://www.bioconductor.org/packages/devel/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#introduction
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(ReactomePA)
library('Rgraphviz')
argv <- commandArgs(TRUE)
data <- read.table(argv[1],header=F)
data2 = bitr(data[,1], fromType="SYMBOL", toType="ENTREZID", annoDb="org.Hs.eg.db")

#GO
goMF <- enrichGO(gene=data2$"ENTREZID",organism="human",ont= "MF",pAdjustMethod = "fdr",pvalueCutoff=1,qvalueCutoff=1, readable=TRUE)
goCC <- enrichGO(gene=data2$"ENTREZID",organism="human",ont= "CC",pAdjustMethod = "fdr",pvalueCutoff=1,qvalueCutoff=1, readable=TRUE)
goBP <- enrichGO(gene=data2$"ENTREZID",organism="human",ont= "BP",pAdjustMethod = "fdr",pvalueCutoff=1,qvalueCutoff=1, readable=TRUE)
write.table(summary(goMF),file="GO-enrich-MF.txt",quote=F,sep="\t",row.names=F)
write.table(summary(goCC),file="GO-enrich-CC.txt",quote=F,sep="\t",row.names=F)
write.table(summary(goBP),file="GO-enrich-BP.txt",quote=F,sep="\t",row.names=F)
pdf("Enrich-analysis.pdf",width=8,height=4,family="GB1")
dotplot(goMF,showCategory=10,title="GO_MF")
dotplot(goCC,showCategory=10,title="GO_CC")
dotplot(goBP,showCategory=10,title="GO_BP")

#KEGG
kk <- enrichKEGG(data2$"ENTREZID",organism='human',pvalueCutoff=1,qvalueCutoff = 1,pAdjustMethod = "fdr",readable=T,use_internal_data = TRUE)
write.table(summary(kk),file="KEGG-enrich.txt",quote=F,sep="\t",row.names=F)
dotplot(kk,showCategory=10,title="KEGG")

#Disese ontology
DO=enrichDO(data2$"ENTREZID", ont = "DO", pvalueCutoff = 1, pAdjustMethod = "fdr",minGSSize = 1, qvalueCutoff = 1,readable = TRUE)
write.table(summary(DO),file="DO-enrich.txt",quote=F,sep="\t",row.names=F)
dotplot(DO,showCategory=10,title="Disese.ontology")

#Rectome Pathway
Pathway=enrichPathway(data2$"ENTREZID", pvalueCutoff = 1, pAdjustMethod = "fdr",minGSSize = 1, qvalueCutoff = 1,readable = TRUE)
write.table(summary(Pathway),file="Pathway-enrich.txt",quote=F,sep="\t",row.names=F)
dotplot(Pathway,showCategory=10,title="Rectome.Pathway")
dev.off()
