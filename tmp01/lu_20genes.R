library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)

setwd("~/gse162498/lu")
Naive_CD4_T<-readRDS("Naive_CD4_T.rds")
CD8_T<-readRDS("CD8_T.rds")
a<-read.csv("b")
a1<-t(a)
markers.to.plot<-c("ACTG1","ACTN1","ADGRE5","ADGRG1","ADORA2A","AHNAK","ANXA1","ANXA6","APMAP","APOL2","ARL4A","ARMH1","ARRDC2","ATF3","ATF4","BIRC3","BTG2","C12orf75","CASP8","CCL3","CCL3L3","CCL5","CCNL1","CCR5","CCR7","CCR8","CD127","CD103","CD25","CD27","CD274","CD28","CD3","CD38","CD4","CD44","CD45RA","CD45RO","CD55","CD57","CD6","CD63","CD69","CD7","CD8","CD83","CD8A","CD8B","CLK1","CMC1","CREM","CSRNP1","CST7","CTSW","CXCL10","CXCL11","CXCL9","CYBA","CYCS","CYTIP","DDIT3","DDX3X","DGKA","DNAJA1","DNAJB1","DUSP1","DUSP2","EDIL3","EEF1B2","EFHD2","EGR1","EIF1","EIF4A2","EIF4A3","EIF4G2","EIF5","EOMES","EZR","FAM102A","FAM177A1","FCGR3A","FCMR","FCRL6","FGR","FLNA","FOS","FOSB","FOXP3","FTH1","FUS","GABRA2","GP5","GPR183","GZMA","GZMB","GZMH","GZMK","H3-3B","HLA-A","HLA-B","HNRNPH1","HOPX","HSP90AA1","HSPA1A","HSPA1B","HSPA5","HSPA6","HSPA8","HSPD1","HSPH1","ICOS","IDO1","IER5","IFI44L","IFI6","IFNG","IL12RB1","IL18RAP","IL32","IL7R","IRF7","ISG15","ISG20","ITGAL","ITGB2","ITIH5","IVNS1ABP","JAML","JUNB","JUND","KDM6B","KLF2","KLF6","KLHL6","KLRD1","KLRF1","LAG3","LAIR2","LAYN","LDHA","LDHB","LEF1","LMNA","LSR","LTB","LYZ","MAL","MCL1","MCM5","MT2A","MX1","MX2","MYADM","MYC","MYO1F","MYOM2","NELL2","NEU1","NFKBIA","NKG7","NR4A1","NXF1","OAF","OXNAD1","PABPC1","PD1","PDCD1","PDCD4","PDL1","PFKFB3","PFN1","PIK3IP1","PLA2G2D","PLEK","PLK3","PMAIP1","PNRC1","PPP1R15A","PPP1R15B","PRF1","PRSS23","PTGER4","PTP4A1","RAP1GAP2","RARRES3","RGCC","RGS1","RGS2","RHOB","RPLP2","RPS26","SARAF","SAT1","SDCBP","SELL","SESN1","SF1","SFPQ","SGK1","SIK1","SLC2A3","SLC4A10","SMAP2","SOCS3","SOX10","SPOCK2","SPON2","SRGN","SRSF3","SRSF5","SRSF6","SRSF7","SSBP4","STAT1","SYNE1","SYNGR3","TAGAP","TAMALIN","TBET","TCF7","TIGIT","TIM3","TLE5","TLR3","TMEM123","TNFAIP3","TNFRSF18","TNFRSF9","TOB1","TOP2A","TRA2B","TRABD2A","TRBV10-2","TRBV28","TRDC","TRDV3","TRGC2","TRGV9","TRNAU1AP","TSC22D3","TSPYL2","TTC38","TUBA4A","TXK","TYROBP","UBB","UBE2S","USP36","VIM","XAF1","XCL1","XCL2","YBX3","YPEL5","ZBTB37","ZFAND5","ZFP36","ZFP36L2","ZNF331")
c<-DotPlot(CD8_T,features = markers.to.plot)
d<-DotPlot(Naive_CD4_T,features = markers.to.plot)
c1<-c$data
d1<-d$data
write.table(c1,"c1")
write.table(d1,"d1")


a<-read.csv("b1")
a1<-t(a)
write.table(a1,"a1")
markers.to.plot<-c("CD69","CD8A" ,"CD8B" ,"FOS" ,"GZMA" ,"GZMH" ,"MYO1F" ,"SYNE1","TSC22D3" ,"TSPYL2" ,"C12orf75" ,"CCR7", "CD27" ,"CD6" ,"CST7" ,"CTSW", "DUSP2", "FCMR" ,"GZMK" ,"HOPX" ,"HSPD1", "IL7R" ,"KLF2", "LEF1", "LTB" ,"MAL", "SAT1" ,"TAGAP" ,"CCL5" ,"DUSP1" ,"EZR" ,"KLF6", "NFKBIA" ,"NKG7" ,"PRF1")


CD8_T<-readRDS("CD8_T01.rds")
a<-DotPlot(CD8_T,features = markers.to.plot)
a1<-a$data
write.table(a1,"a1_cd8T")

ex_CD8_T<-readRDS("Exhausted_CD8.rds")
a<-DotPlot(ex_CD8_T,features = markers.to.plot)
a1<-a$data
write.table(a1,"a1_ex_cd8T")

ex_CD48_T<-readRDS("Exhausted_CD48.rds")
a<-DotPlot(ex_CD48_T,features = markers.to.plot)
a1<-a$data
write.table(a1,"a1_ex_cd48T")

Naive_CD8<-readRDS("Naive_CD8.rds")
a<-DotPlot(Naive_CD8,features = markers.to.plot)
a1<-a$data
write.table(a1,"a1_Naive_CD8")