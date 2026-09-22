library(tximport)
library(DESeq2)
library(pheatmap)
library(dendextend)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(org.Dr.eg.db)


#### DESeq2 + PCA ####

fileenamesIso <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/Juex/me_analysis/isoforms/MIR", pattern = "*isoforms.results", full.names = TRUE)

tx2gene <- read.csv(fileenamesIso[1], sep = "\t", stringsAsFactors = F)
condition<-c(rep("AmiRWTsham",4),rep("BmiRWTcryo",4),rep("CmiRMutsham",4),rep("DmiRMutcryo",4))
reps<-rep(c("rep1","rep2","rep3","rep4"),times=c(4,4,4,4))


data<-data.frame(condition=condition, reps=reps, stringsAsFactors = T)
samplenames <- c(paste(rep("AmiRWTsham", 4), 1:4, sep = ""), 
                 paste(rep("BmiRWTcryo", 4), 1:4, sep = ""), 
                 paste(rep("CmiRMutsham", 4), 1:4, sep = ""), 
                 paste(rep("DmiRMutcryo", 4), 1:4, sep = ""))
row.names(data) <- samplenames


txi <- tximport(fileenamesIso, type = "rsem", tx2gene = tx2gene[,1:2])
#design is unpaired, all distinct animals and heart pools
dds <- DESeqDataSetFromTximport(txi, data, design=~condition)
dds <- dds[rowSums(counts(dds))>10, ]
dds <- DESeq(dds)

#save for github:
#saveRDS(dds, "dds_miRzdKO_manuscript.rds")

rld <- rlog(dds, blind=FALSE)

data <- plotPCA(rld,intgroup=c("condition"),ntop=500,returnData=TRUE)
names = rownames(colData(rld))
percentVar <- round(100 * attr(data, "percentVar"))

data$condition <- sapply(sapply(as.character(data$condition), strsplit, "miR"), "[[", 2)

data$conditionII <- data$condition
data$conditionII <- gsub("Mut", "miR-dKo ", data$conditionII)
data$conditionII <- gsub("WT", "WT ", data$conditionII)

data$Genotype <- "WT"
data$Genotype[grepl("Mut", data$condition)] <- "miR-zdKO"

data$Injury <- "Sham"
data$Injury[grepl("cryo", data$conditionII)] <- "Cryo-injury"

ggplot(data, aes(PC1, PC2, shape=Genotype, color=Injury)) +
  geom_point(size=5, alpha = 0.8) +
  scale_color_manual(values = c(`Cryo-injury` = "blueviolet", Sham = "chartreuse4")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_x_continuous(limits = c(-20, 20)) +
  theme_bw() +
  theme(text = element_text(size=24))

#### FPKM ####
filenames <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/Juex/me_analysis/gene_results/miR/", pattern = "*genes.results", full.names = TRUE)
samplenames <- c(paste(rep("AmiRWTsham", 4), 1:4, sep = ""), 
                 paste(rep("BmiRWTcryo", 4), 1:4, sep = ""), 
                 paste(rep("CmiRMutsham", 4), 1:4, sep = ""), 
                 paste(rep("DmiRMutcryo", 4), 1:4, sep = ""))

genes<-read.table(filenames[1],header=TRUE,sep="\t",stringsAsFactors = FALSE)[,1]
fpkm <-do.call(cbind,lapply(filenames,function(fn)read.table(fn,header=TRUE,sep="\t",stringsAsFactors = FALSE)[,7]))
fpkm <- data.frame(genes,fpkm,stringsAsFactors = FALSE)
colnames(fpkm) <- c("ENSEMBL",samplenames)

#mean and standard error
fpkm_list_AmiRWTsham <- as.list(as.data.frame(t(fpkm[,2:5])))
fpkm_list_BmiRWTcryo <- as.list(as.data.frame(t(fpkm[,6:9])))
fpkm_list_CmiRWTcryo <- as.list(as.data.frame(t(fpkm[,10:13])))
fpkm_list_DmiRWTcryo <- as.list(as.data.frame(t(fpkm[,14:17])))

fpkm_mean_AmiRWTsham <- sapply(fpkm_list_AmiRWTsham, mean)
fpkm_mean_BmiRWTcryo <- sapply(fpkm_list_BmiRWTcryo, mean)
fpkm_mean_CmiRWTcryo <- sapply(fpkm_list_CmiRWTcryo, mean)
fpkm_mean_DmiRWTcryo <- sapply(fpkm_list_DmiRWTcryo, mean)

std <- function(x) sd(x)/sqrt(length(x))
fpkm_se_AmiRWTsham <- sapply(fpkm_list_AmiRWTsham, std)
fpkm_se_BmiRWTcryo <- sapply(fpkm_list_BmiRWTcryo, std)
fpkm_se_CmiRWTcryo <- sapply(fpkm_list_CmiRWTcryo, std)
fpkm_se_DmiRWTcryo <- sapply(fpkm_list_DmiRWTcryo, std)

fpkm_mean_treatment <- data.frame("WTsham_meanFPKM" = fpkm_mean_AmiRWTsham, "WTsham_seFPKM" = fpkm_se_AmiRWTsham,
                                  "WTcryo_meanFPKM" = fpkm_mean_BmiRWTcryo, "WTcryo_seFPKM" = fpkm_se_BmiRWTcryo,
                                  "MutSham_meanFPKM" = fpkm_mean_CmiRWTcryo, "MutSham_seFPKM" = fpkm_se_CmiRWTcryo,
                                  "Mutcryo_meanFPKM" = fpkm_mean_DmiRWTcryo, "Mutcryo_seFPKM" = fpkm_se_DmiRWTcryo)

trial <- as.list(as.data.frame(t(fpkm_mean_treatment)))
fpkm_max_treatment <- as.numeric(sapply(trial, max))
fpkmDetail <- cbind(fpkm, fpkm_mean_treatment, fpkm_max_treatment)




#### Q1, extent of changes between WT and mutant in sham: ####

#PCA suggests minor changes only, number of DEGs:

#obtain default LFCs and p values:
sham_WT_mut_res <- as.data.frame(results(dds, contrast=c("condition","CmiRMutsham","AmiRWTsham")))

#filter to DEGs:
exprs_sham <- filter(fpkmDetail, (MutSham_meanFPKM>1 | WTsham_meanFPKM>1))

#20 passing our threshold
sham_WT_mut_DE <- filter(sham_WT_mut_res, rownames(sham_WT_mut_res) %in% exprs_sham$ENSEMBL, padj <0.05, abs(log2FoldChange) > log2(1.5))
#51 if lower a bit
sham_WT_mut_DE_weak <- filter(sham_WT_mut_res, rownames(sham_WT_mut_res) %in% exprs_sham$ENSEMBL, padj <0.1, abs(log2FoldChange) > log2(1.1))

#top upreg is CARMN

#could re-check with just the relevant samples (high variability in others may reduce power) but probs unnecessary

#20 passing our threshold
sham_WT_mut_exprs <- filter(sham_WT_mut_res, rownames(sham_WT_mut_res) %in% exprs_sham$ENSEMBL)

sham_WT_mut_exprs$DE <- "No"
sham_WT_mut_exprs$DE[rownames(sham_WT_mut_exprs) %in% rownames(sham_WT_mut_DE)] <- "DE"

ggplot(sham_WT_mut_exprs) + aes(x = log2FoldChange, y = -log10(padj), color = DE) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") + 
  #coord_cartesian(xlim = c(-3.9,3.9), ylim = c(0,20)) +
  scale_color_manual(values = c("No" = "grey60", "DE" = "steelblue")) +
  theme_minimal() +
  xlab("LFC miR-dKo vs. WT\n(cryo-injury)") +
  ylab("-log10(p)")


#note high fold changes:
#CARMN is strongest up-regulated gene and by far the most significant change
#cyfip2 and ablim3 also very highly upregulated - they are neighbours and approx 5-6MB upstream of CARMNs
#afap1|1a next down the list, next to ablim3 
#worth doing a full look at chr14 changing genes when full annotation is done later


#
#### Q2a, extent of changes between WT and mutant in cryo: ####

#PCA suggests minor changes only, number of DEGs:

#obtain default LFCs and p values:
cryo_WT_mut_res <- as.data.frame(results(dds, contrast=c("condition","DmiRMutcryo","BmiRWTcryo")))

#filter to DEGs:
exprs_cryo <- filter(fpkmDetail, (Mutcryo_meanFPKM>1 | WTcryo_meanFPKM>1))

#239 passing our threshold
cryo_WT_mut_exprs <- filter(cryo_WT_mut_res, rownames(cryo_WT_mut_res) %in% exprs_cryo$ENSEMBL)
cryo_WT_mut_DE <- filter(cryo_WT_mut_res, rownames(cryo_WT_mut_res) %in% exprs_cryo$ENSEMBL, padj <0.05, abs(log2FoldChange) > log2(1.5))

cryo_WT_mut_exprs$DE <- "No"
cryo_WT_mut_exprs$DE[rownames(cryo_WT_mut_exprs) %in% rownames(cryo_WT_mut_DE)] <- "DE"

ggplot(cryo_WT_mut_exprs) + aes(x = log2FoldChange, y = -log10(padj), color = DE) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") + 
  coord_cartesian(xlim = c(-4.5,5), #ylim = c(0,15)
  ) +
  scale_color_manual(values = c("No" = "grey60", "DE" = "steelblue")) +
  theme_minimal() +
  xlab("LFC miR-dKo vs. WT\n(cryo-injury)") +
  ylab("-log10(p)")

#note high fold changes:
#CARMN is again strongest up-regulated gene and by far the most significant change
#cyfip2 and ablim3 and afap1l1a also same pattern though less significance
#cyfip2 particularly strong change still, others bit weaker


#
#### Q2b, unique changes/genotype-specific sham-cryo responses between WT and mutant (crude) ####

#crude approach, unique DEGs to either (clustering approach to be trialled after to see if extra useful detail)

#obtain default LFCs and p values:
WT_sham_cryo_res <- as.data.frame(results(dds, contrast=c("condition","BmiRWTcryo","AmiRWTsham")))
Mut_sham_cryo_res <- as.data.frame(results(dds, contrast=c("condition","DmiRMutcryo","CmiRMutsham")))

#filter to DEGs:
exprs_WT <- filter(fpkmDetail, (WTsham_meanFPKM>1 | WTcryo_meanFPKM>1))
exprs_mut <- filter(fpkmDetail, (MutSham_meanFPKM>1 | Mutcryo_meanFPKM>1))

#1351 cryo-responsive genes in WT
WT_sham_cryo_DE <- filter(WT_sham_cryo_res, rownames(WT_sham_cryo_res) %in% exprs_WT$ENSEMBL, 
                          padj <0.05, abs(log2FoldChange) > log2(1.5))
length(unique(rownames(WT_sham_cryo_DE)))

#1769 cryo-responsive genes in Mut
Mut_sham_cryo_DE <- filter(Mut_sham_cryo_res, rownames(Mut_sham_cryo_res) %in% exprs_mut$ENSEMBL, 
                           padj <0.05, abs(log2FoldChange) > log2(1.5))
length(unique(rownames(Mut_sham_cryo_DE)))

#genes up in both:
WT_sham_cryo_DEup <- filter(WT_sham_cryo_DE, log2FoldChange > log2(1.5))#935
Mut_sham_cryo_DEup <- filter(Mut_sham_cryo_DE, log2FoldChange > log2(1.5))#1110

sum(rownames(WT_sham_cryo_DEup) %in% rownames(Mut_sham_cryo_DEup))#524 up in both

(1110-524) + 524 + (935-524)
length(unique(c(rownames(WT_sham_cryo_DEup), rownames(Mut_sham_cryo_DEup))))

BothUp <- rownames(WT_sham_cryo_DEup)[rownames(WT_sham_cryo_DEup) %in% rownames(Mut_sham_cryo_DEup)]

#genes down in both:
WT_sham_cryo_DEdown <- filter(WT_sham_cryo_DE, log2FoldChange < -log2(1.5))#416
Mut_sham_cryo_DEdown <- filter(Mut_sham_cryo_DE, log2FoldChange < -log2(1.5))#659

sum(rownames(WT_sham_cryo_DEdown) %in% rownames(Mut_sham_cryo_DEdown))#129 down in both

(416-129) + 129 + (659-129)#946
length(unique(c(rownames(WT_sham_cryo_DEdown), rownames(Mut_sham_cryo_DEdown))))

BothDown <- rownames(WT_sham_cryo_DEdown)[rownames(WT_sham_cryo_DEdown) %in% rownames(Mut_sham_cryo_DEdown)]

#hence, overlap of 653 genes between both
524+129

DownUp <- rownames(WT_sham_cryo_DEdown)[rownames(WT_sham_cryo_DEdown) %in% rownames(Mut_sham_cryo_DEup)]
UpDown <- rownames(WT_sham_cryo_DEup)[rownames(WT_sham_cryo_DEup) %in% rownames(Mut_sham_cryo_DEdown)]

1351-653-6 #698 unique gene regulations WT (6 are opposing)
1769-653-6 #1116 unique changing genes mut (6 are opposing)

653+6+692+1110
length(unique(c(rownames(WT_sham_cryo_DE), rownames(Mut_sham_cryo_DE))))

#figue it out:

cryo_responsive_DE <- filter(fpkmDetail, ENSEMBL %in% c(rownames(WT_sham_cryo_DE), rownames(Mut_sham_cryo_DE)))

cryo_responsive_DE$UpWT[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEup)] <- "Yes"
cryo_responsive_DE$DownWT[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEdown)] <- "Yes"
cryo_responsive_DE$UpMut[cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEup)] <- "Yes"
cryo_responsive_DE$DownMut[cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEdown)] <- "Yes"


cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEup)] <- "UpWTOnly"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEdown)] <- "DownWTOnly"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEup)] <- "UpMutOnly"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEdown)] <- "DownMutOnly"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEup) &
                            cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEup)] <- "UpBoth"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEdown) &
                            cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEdown)] <- "DownBoth"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEup) &
                            cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEdown)] <- "Opposing"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEdown) &
                            cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEup)] <- "Opposing"
table(cryo_responsive_DE$Shared)
sum(table(cryo_responsive_DE$Shared))
sum(table(cryo_responsive_DE$Shared)[c(2,5)])
sum(table(cryo_responsive_DE$Shared)[c(2,5)])

#uniquely changing genes therefore heavily outweigh shared changes

#enriched in WT only (possible miR targets):
WT_sham_cryo_DEupOnly <- filter(WT_sham_cryo_DEup, !rownames(WT_sham_cryo_DEup) %in% rownames(Mut_sham_cryo_DEup))



#### Q2bi, genotype-specific response via an interaction term ####

ddsInt <- DESeqDataSetFromTximport(txi, data, design=~Genotype + Injury + Genotype:Injury)
ddsInt <- ddsInt[rowSums(counts(ddsInt))>10, ]
ddsInt <- DESeq(ddsInt)

resultsNames(ddsInt)
resultsNames(dds)

#for default genotype of Mut (alphabetically before WT):
Mut_sham_cryo_resInt <- as.data.frame(results(ddsInt, contrast=c("Injury","Cryo-injury","Sham")))
#the p vals look essentially identical?
dim(Mut_sham_cryo_resInt)
#dim(Mut_sham_cryo_res)

#now for WT
WT_sham_cryo_resInt <- as.data.frame(results(ddsInt, list( c("Injury_Sham_vs_Cryo.injury","GenotypeWT.InjurySham") ) ))
#again, identical, with the sign change as expected purely due to alphabetical listing
WT_sham_cryo_resInt$log2FoldChange <- -WT_sham_cryo_resInt$log2FoldChange

# the interaction term, answering: is the cryo-injury effect *different* across genotypes?
GenotypeInjuryInteraction <- as.data.frame(results(ddsInt, name = "GenotypeWT.InjurySham"))

#number of genes with a good p:
dim(filter(GenotypeInjuryInteraction, padj > 0.05))
dim(filter(GenotypeInjuryInteraction, padj < 0.05)) 
#153, but how many are also sig/decent fold change in individual genotypes sham - cryo?

#example of genes diff across genotypes:
plotCounts(ddsInt, "ENSDARG00000092976", intgroup=c("Genotype", "Injury"))#loss in wt, gain in mut
plotCounts(ddsInt, "ENSDARG00000058940", intgroup=c("Genotype", "Injury"))#loss in wt, stable in mut
plotCounts(ddsInt, "ENSDARG00000061120", intgroup=c("Genotype", "Injury"))#loss in wt, stable in mut

plotCounts(ddsInt, "ENSDARG00000109596", intgroup=c("Genotype", "Injury"))#stable in wt, gained in mut

#this is indeed useful to have a p val to show these genes, tho it seems v stringent (only 98...?)

#filter to DEGs:
exprs_WT <- filter(fpkmDetail, (WTsham_meanFPKM>1 | WTcryo_meanFPKM>1))
exprs_mut <- filter(fpkmDetail, (MutSham_meanFPKM>1 | Mutcryo_meanFPKM>1))

#1351 cryo-responsive genes in WT
WT_sham_cryo_DE <- filter(WT_sham_cryo_resInt, rownames(WT_sham_cryo_resInt) %in% exprs_WT$ENSEMBL, 
                          padj <0.05, abs(log2FoldChange) > log2(1.5))
length(unique(rownames(WT_sham_cryo_DE)))

#1769 cryo-responsive genes in Mut
Mut_sham_cryo_DE <- filter(Mut_sham_cryo_resInt, rownames(Mut_sham_cryo_resInt) %in% exprs_mut$ENSEMBL, 
                           padj <0.05, abs(log2FoldChange) > log2(1.5))
length(unique(rownames(Mut_sham_cryo_DE)))

#genes up in both:
WT_sham_cryo_DEup <- filter(WT_sham_cryo_DE, log2FoldChange > log2(1.5))#935
Mut_sham_cryo_DEup <- filter(Mut_sham_cryo_DE, log2FoldChange > log2(1.5))#1110

sum(rownames(WT_sham_cryo_DEup) %in% rownames(Mut_sham_cryo_DEup))#524 up in both

(1110-524) + 524 + (935-524)
length(unique(c(rownames(WT_sham_cryo_DEup), rownames(Mut_sham_cryo_DEup))))

BothUp <- rownames(WT_sham_cryo_DEup)[rownames(WT_sham_cryo_DEup) %in% rownames(Mut_sham_cryo_DEup)]

#genes down in both:
WT_sham_cryo_DEdown <- filter(WT_sham_cryo_DE, log2FoldChange < -log2(1.5))#416
Mut_sham_cryo_DEdown <- filter(Mut_sham_cryo_DE, log2FoldChange < -log2(1.5))#659

sum(rownames(WT_sham_cryo_DEdown) %in% rownames(Mut_sham_cryo_DEdown))#129 down in both

(416-129) + 129 + (659-129)#946
length(unique(c(rownames(WT_sham_cryo_DEdown), rownames(Mut_sham_cryo_DEdown))))

BothDown <- rownames(WT_sham_cryo_DEdown)[rownames(WT_sham_cryo_DEdown) %in% rownames(Mut_sham_cryo_DEdown)]

#hence, overlap of 653 genes between both
524+129

DownUp <- rownames(WT_sham_cryo_DEdown)[rownames(WT_sham_cryo_DEdown) %in% rownames(Mut_sham_cryo_DEup)]
UpDown <- rownames(WT_sham_cryo_DEup)[rownames(WT_sham_cryo_DEup) %in% rownames(Mut_sham_cryo_DEdown)]

1351-653-6 #698 unique gene regulations WT (6 are opposing)
1769-653-6 #1116 unique changing genes mut (6 are opposing)

653+6+692+1110
length(unique(c(rownames(WT_sham_cryo_DE), rownames(Mut_sham_cryo_DE))))

#label all
cryo_responsive_DE <- filter(fpkmDetail, ENSEMBL %in% c(rownames(WT_sham_cryo_DE), rownames(Mut_sham_cryo_DE)))

cryo_responsive_DE$UpWT[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEup)] <- "Yes"
cryo_responsive_DE$DownWT[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEdown)] <- "Yes"
cryo_responsive_DE$UpMut[cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEup)] <- "Yes"
cryo_responsive_DE$DownMut[cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEdown)] <- "Yes"


cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEup)] <- "UpWTOnly"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEdown)] <- "DownWTOnly"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEup)] <- "UpMutOnly"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEdown)] <- "DownMutOnly"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEup) &
                            cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEup)] <- "UpBoth"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEdown) &
                            cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEdown)] <- "DownBoth"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEup) &
                            cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEdown)] <- "Opposing"
cryo_responsive_DE$Shared[cryo_responsive_DE$ENSEMBL %in% rownames(WT_sham_cryo_DEdown) &
                            cryo_responsive_DE$ENSEMBL %in% rownames(Mut_sham_cryo_DEup)] <- "Opposing"
table(cryo_responsive_DE$Shared)
sum(table(cryo_responsive_DE$Shared))
sum(table(cryo_responsive_DE$Shared)[c(2,5)])
sum(table(cryo_responsive_DE$Shared)[c(2,5)])


#524 both up genes, could be some here:
sum(BothUp %in% rownames(filter(GenotypeInjuryInteraction, padj > 0.05)))
522/524
#0 are significant per genotype
sum(BothUp %in% rownames(filter(GenotypeInjuryInteraction, padj < 0.05)))

#129 both down genes, none in there
sum(BothDown %in% rownames(filter(GenotypeInjuryInteraction, padj > 0.05)))
sum(BothDown %in% rownames(filter(GenotypeInjuryInteraction, padj < 0.05)))

cryo_responsive_DE$IntTerm <- "Non-significant"
#cryo_responsive_DE$IntTerm[cryo_responsive_DE$ENSEMBL %in% 
#                             rownames(filter(GenotypeInjuryInteraction, padj < 0.1))] <- "Weak Significant"
cryo_responsive_DE$IntTerm[cryo_responsive_DE$ENSEMBL %in% 
                             rownames(filter(GenotypeInjuryInteraction, padj < 0.05))] <- "Significant"

table(cryo_responsive_DE$Shared)
table(filter(cryo_responsive_DE, IntTerm == "Significant")$Shared) #primarily the exclusive genes in the venn
#table(filter(cryo_responsive_DE, IntTerm == "Weak Significant")$Shared) #primarily the exclusive genes in the venn


#stacked bar plot rather than venn, up and down genes in each genotype
#stacks of light to deep red

#cryo-injury responsive gene sets:
#a) induced in both genotypes
#b) induced in 1 genotype 
#c) induced in 1 genotype (genotype:injury p < 0.05)

#d-f equivalent for repressed

WT_sham_cryo_DE$InducedRepressed <- "Induced"
WT_sham_cryo_DE$InducedRepressed[WT_sham_cryo_DE$log2FoldChange <0] <- "Repressed"

WT_sham_cryo_DE$GenotypeSpecificityLevel[WT_sham_cryo_DE$InducedRepressed == "Induced" &
                                           rownames(WT_sham_cryo_DE) %in% rownames(Mut_sham_cryo_DEup)] <- "Induced in both genotypes"
WT_sham_cryo_DE$GenotypeSpecificityLevel[WT_sham_cryo_DE$InducedRepressed == "Induced" &
                                           !rownames(WT_sham_cryo_DE) %in% rownames(Mut_sham_cryo_DEup)] <- "Induced in 1 genotype only"
WT_sham_cryo_DE$GenotypeSpecificityLevel[WT_sham_cryo_DE$InducedRepressed == "Induced" &
                                           !rownames(WT_sham_cryo_DE) %in% rownames(Mut_sham_cryo_DEup) &
                                           rownames(WT_sham_cryo_DE) %in% rownames(filter(GenotypeInjuryInteraction, padj < 0.05))] <- "Induced in 1 genotype only (genotype:injury p < 0.05)"

WT_sham_cryo_DE$GenotypeSpecificityLevel[WT_sham_cryo_DE$InducedRepressed == "Repressed" &
                                           rownames(WT_sham_cryo_DE) %in% rownames(Mut_sham_cryo_DEdown)] <- "Repressed in both genotypes"
WT_sham_cryo_DE$GenotypeSpecificityLevel[WT_sham_cryo_DE$InducedRepressed == "Repressed" &
                                           !rownames(WT_sham_cryo_DE) %in% rownames(Mut_sham_cryo_DEdown)] <- "Repressed in 1 genotype only"
WT_sham_cryo_DE$GenotypeSpecificityLevel[WT_sham_cryo_DE$InducedRepressed == "Repressed" &
                                           !rownames(WT_sham_cryo_DE) %in% rownames(Mut_sham_cryo_DEdown) &
                                           rownames(WT_sham_cryo_DE) %in% rownames(filter(GenotypeInjuryInteraction, padj < 0.05))] <- "Repressed in 1 genotype only (genotype:injury p < 0.05)"
table(WT_sham_cryo_DE$GenotypeSpecificityLevel)

Mut_sham_cryo_DE$InducedRepressed <- "Induced"
Mut_sham_cryo_DE$InducedRepressed[Mut_sham_cryo_DE$log2FoldChange <0] <- "Repressed"

Mut_sham_cryo_DE$GenotypeSpecificityLevel[Mut_sham_cryo_DE$InducedRepressed == "Induced" &
                                            rownames(Mut_sham_cryo_DE) %in% rownames(WT_sham_cryo_DEup)] <- "Induced in both genotypes"
Mut_sham_cryo_DE$GenotypeSpecificityLevel[Mut_sham_cryo_DE$InducedRepressed == "Induced" &
                                            !rownames(Mut_sham_cryo_DE) %in% rownames(WT_sham_cryo_DEup)] <- "Induced in 1 genotype only"
Mut_sham_cryo_DE$GenotypeSpecificityLevel[Mut_sham_cryo_DE$InducedRepressed == "Induced" &
                                            !rownames(Mut_sham_cryo_DE) %in% rownames(WT_sham_cryo_DEup) &
                                            rownames(Mut_sham_cryo_DE) %in% rownames(filter(GenotypeInjuryInteraction, padj < 0.05))] <- "Induced in 1 genotype only (genotype:injury p < 0.05)"

Mut_sham_cryo_DE$GenotypeSpecificityLevel[Mut_sham_cryo_DE$InducedRepressed == "Repressed" &
                                            rownames(Mut_sham_cryo_DE) %in% rownames(WT_sham_cryo_DEdown)] <- "Repressed in both genotypes"
Mut_sham_cryo_DE$GenotypeSpecificityLevel[Mut_sham_cryo_DE$InducedRepressed == "Repressed" &
                                            !rownames(Mut_sham_cryo_DE) %in% rownames(WT_sham_cryo_DEdown)] <- "Repressed in 1 genotype only"
Mut_sham_cryo_DE$GenotypeSpecificityLevel[Mut_sham_cryo_DE$InducedRepressed == "Repressed" &
                                            !rownames(Mut_sham_cryo_DE) %in% rownames(WT_sham_cryo_DEdown) &
                                            rownames(Mut_sham_cryo_DE) %in% rownames(filter(GenotypeInjuryInteraction, padj < 0.05))] <- "Repressed in 1 genotype only (genotype:injury p < 0.05)"
table(Mut_sham_cryo_DE$GenotypeSpecificityLevel)

WT_sham_cryo_DE$Genotype <- "WT"
Mut_sham_cryo_DE$Genotype <- "miR-zdKO"

plot_injuryDEGs_injurySpecific <- rbind(WT_sham_cryo_DE[,c(7:9)], 
                                        Mut_sham_cryo_DE[,c(7:9)])

plot_injuryDEGs_injurySpecific$GenotypeSpecificityLevel <- as.character(plot_injuryDEGs_injurySpecific$GenotypeSpecificityLevel)
plot_injuryDEGs_injurySpecific$GenotypeSpecificityLevel <- as.factor(plot_injuryDEGs_injurySpecific$GenotypeSpecificityLevel)
plot_injuryDEGs_injurySpecific$GenotypeSpecificityLevel <- factor(plot_injuryDEGs_injurySpecific$GenotypeSpecificityLevel,
                                                                  levels(plot_injuryDEGs_injurySpecific$GenotypeSpecificityLevel)[c(2,1,3,5,4,6
                                                                                                                                    #5,4,6,2,1,3
                                                                  )])

plot_injuryDEGs_injurySpecific$Genotype <- as.factor(plot_injuryDEGs_injurySpecific$Genotype)
plot_injuryDEGs_injurySpecific$Genotype <- factor(plot_injuryDEGs_injurySpecific$Genotype,
                                                  levels(plot_injuryDEGs_injurySpecific$Genotype)[c(2,1)])

ggplot(plot_injuryDEGs_injurySpecific) + aes(x = InducedRepressed, fill = GenotypeSpecificityLevel) +
  geom_bar(position = "stack") +
  facet_wrap(~Genotype) +
  scale_fill_manual(values = c(`Induced in both genotypes` = "plum1",
                               `Induced in 1 genotype only` = "mediumpurple1",
                               `Induced in 1 genotype only (genotype:injury p < 0.05)` = "blueviolet",
                               `Repressed in both genotypes` = "chartreuse",
                               `Repressed in 1 genotype only` = "chartreuse3",
                               `Repressed in 1 genotype only (genotype:injury p < 0.05)` = "chartreuse4")) +
  ylab("No. Genes") +
  xlab("") +
  theme_bw() +
  theme(text = element_text(size =20),
        axis.text.x = element_text(size =14)
  )


#
#### Q3, what are the uniquely changing genes representative of - cluster all cryo-response and put into context ####

#cluster all cryo-injury response genes, isolate genotype-specific responses within the broader clustering

cryo_responsive_DE <- filter(fpkmDetail, ENSEMBL %in% c(rownames(WT_sham_cryo_DE), rownames(Mut_sham_cryo_DE)))

#optional, just do the genotype:injury genes
cryo_responsive_DE <- filter(fpkmDetail, ENSEMBL %in% c(rownames(WT_sham_cryo_DE), rownames(Mut_sham_cryo_DE)),
                             ENSEMBL %in% rownames(filter(GenotypeInjuryInteraction, padj < 0.05)))

trial <- unique(c(rownames(WT_sham_cryo_DE), rownames(Mut_sham_cryo_DE)))

trial[!trial %in% cryo_responsive_DE$ENSEMBL]

mat <- cryo_responsive_DE[,2:17] #2461 genes
rownames(mat) <- cryo_responsive_DE[,1]
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
mat <- t(apply(mat, 1, cal_z_score))

myColor <- colorRampPalette(c("steelblue", "white", "red"))(50)
myBreaks <- c(seq(min(mat), 0, 
                  length.out=ceiling(50/2)), 
              seq(max(mat)/50, 
                  max(mat), 
                  length.out=floor(50/2)))


data_cols<-data

data_cols$genotype <- "WT"
data_cols$genotype[grepl("Mut", data_cols$condition)] <- "miR_ko"
data_cols$cryoInjury <- "Sham"
data_cols$cryoInjury[grepl("cryo", data_cols$condition)] <- "Cryo"

rownames(data_cols) == colnames(mat)

#data_colsHeat$Hours <- factor(data_colsHeat$Hours, levels(data_colsHeat$Hours)[c(1,3,4,2)])

my_colour_anno = list(
  Injury = c(`Cryo-injury` = "blueviolet", Sham = "chartreuse4"),
  Genotype = c(WT = "coral1", `miR-zdKO` = "steelblue")
)

library(pheatmap)
p <- pheatmap(mat,annotation_colors = my_colour_anno,gaps_col = c(8),
              annotation_col = data.frame("Injury" = data_cols[,9], row.names = rownames(data_cols)),
              show_colnames = F, 
              show_rownames = F, 
              cluster_cols = F,
              cluster_rows = T,
              cutree_rows = 4,#cant split C2 even with 12 clusters
              treeheight_col = 0, 
              treeheight_row = 45,
              color = myColor, 
              breaks = myBreaks,
              border_color = NA)


#pattern per cluster:
mat_dist <- hclust(dist(mat))
library(dendextend)
mat_dist_cut7 <- cut_lower_fun(as.dendrogram(mat_dist), h = 6.8)
#6.8 gives 7 clusters with all genes

mat_dist_cut7<- cut_lower_fun(as.dendrogram(mat_dist), h = 6)
#6 gives 4 with interaction genes

Cluster1 <- mat_dist_cut7[[1]]#
Cluster2 <- mat_dist_cut7[[2]]#
Cluster3 <- mat_dist_cut7[[3]]#
Cluster4 <- mat_dist_cut7[[4]]#
Cluster5 <- mat_dist_cut7[[5]]#
Cluster6 <- mat_dist_cut7[[6]]#
Cluster7 <- mat_dist_cut7[[7]]#

#double check genes match if needed:
miR143_145_cryo_DE <- filter(miR143_145_cryo_RNAseq_summary, !DEG_Cluster == "None")

table(miR143_145_cryo_DE$DEG_Cluster) #match

#general patterns:
Clusters <- list(Cluster1, Cluster2, Cluster3, Cluster4)#, Cluster5, Cluster6, Cluster7)
ClusterPlots <- list()

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#titles - assigned later:
clusterTitles <- c("C1(ECM/Wnt)", "C2(Mitosis/ECM)", "C3(CM/Muscle)",
                   "C4(Respiration)", "C5(Respiration/CM)", "C6", "C7(CM/Muscle)")

clusterTitles <- c("C1", "C2", "C3",
                   "C4")


for (i in 1:4) {
  mat_Cluster1 <- data.frame(t(mat[Clusters[[i]],]))
  #mean per condition:
  trial <- sapply(mat_Cluster1[grep("WTsham", rownames(mat_Cluster1)),], mean)
  triali <- sapply(mat_Cluster1[grep("WTcryo", rownames(mat_Cluster1)),], mean)
  trialii <- sapply(mat_Cluster1[grep("Mutsham", rownames(mat_Cluster1)),], mean)
  trialiii <- sapply(mat_Cluster1[grep("Mutcryo", rownames(mat_Cluster1)),], mean)
  
  trend_Cluster1 <- c(rep("Sham", length(Clusters[[i]])), rep("Cryo-injury", length(Clusters[[i]])))
  trend_Cluster1 <- data.frame("cryoInjury" = factor(trend_Cluster1, levels = levels(as.factor(trend_Cluster1))[c(2,1)]), 
                               "genotype" = "WT",
                               rbind(data.frame("meanZ" = trial), 
                                     data.frame("meanZ" = triali)))
  
  trend_Cluster1 <- rbind(trend_Cluster1, 
                          data.frame("cryoInjury" = factor(trend_Cluster1$cryoInjury, levels = levels(as.factor(trend_Cluster1$cryoInjury))[c(2,1)]), 
                                     "genotype" = "miR-zdKO",
                                     rbind(data.frame("meanZ" = trialii), 
                                           data.frame("meanZ" = trialiii))))
  
  trend_Cluster1$genotype <- factor(trend_Cluster1$genotype)
  trend_Cluster1$genotype <- factor(trend_Cluster1$genotype, levels = levels(trend_Cluster1$genotype)[c(2,1)])
  
  ClusterPlots[[i]] <- ggplot(trend_Cluster1, aes(y=meanZ, x=genotype, fill = cryoInjury)) +
    geom_split_violin(alpha = 0.7) +
    theme_minimal() +
    ggtitle(label =  clusterTitles[i]) +
    #ggtitle(label =  paste("Cluster", i, " - ", length(Clusters[[i]]), " DEGs", sep = "")) +
    scale_fill_manual(values = my_colour_anno[[1]]) +
    scale_y_continuous(breaks = c(-1, 0 , 1)) +
    theme(axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.y = element_text(size =15),
          axis.title.x = element_blank(), 
          title = element_text(size=14),
          legend.position = "none") +
    labs(y = "") 
}

library(grid)
library(gridExtra)
grid.arrange(#top = textGrob("\nEC/nonEC Cluster Profiles\n", gp = gpar(fontface = 'bold', fontsize = 18)),
  #bottom = textGrob("Day\n", gp = gpar(fontface = 'bold', fontsize = 15)),
  #left = textGrob("\nAverage Z score/gene", gp = gpar(fontface = 'bold', fontsize = 15), rot = 90),
  ClusterPlots[[1]], ClusterPlots[[2]], ClusterPlots[[3]], ClusterPlots[[4]],# ClusterPlots[[5]], ClusterPlots[[6]], ClusterPlots[[7]],
  ncol = 4)
#pdf("ETV2_cluster7.pdf", width = 5, height = 17)
#dev.off()

ClusterPlots[[1]]
ClusterPlots[[2]]
ClusterPlots[[3]]
ClusterPlots[[4]]
ClusterPlots[[5]]
ClusterPlots[[7]]

egoall <- list()
egoCC <- list()

background <- filter(fpkmDetail, fpkm_max_treatment>1)$ENSEMBL

for (i in 1:4){
  egoall[[i]] <- enrichGO(gene         = Clusters[[i]],
                          universe      = background,
                          keyType       = "ENSEMBL",
                          OrgDb         = org.Dr.eg.db,
                          ont           = "all",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.1,
                          qvalueCutoff  = 0.1,
                          readable      = TRUE)
}


ego_C1 <- data.frame(egoall[[1]])
ego_C2 <- data.frame(egoall[[2]])
ego_C3 <- data.frame(egoall[[3]])
ego_C4 <- data.frame(egoall[[4]])
ego_C5 <- data.frame(egoall[[5]])
ego_C6 <- data.frame(egoall[[6]])
ego_C7 <- data.frame(egoall[[7]])

#C1 and C2 are broadly, on in mutants - 0 terms
ego_onMutC1C2 <- enrichGO(gene         = unlist(Clusters[1:2]),
                          universe      = background,
                          keyType       = "ENSEMBL",
                          OrgDb         = org.Dr.eg.db,
                          ont           = "all",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.1,
                          qvalueCutoff  = 0.1,
                          readable      = TRUE)

#C3 and C4 are on in WT:
ego_onMutC3C4 <- enrichGO(gene         = unlist(Clusters[3:4]),
                          universe      = background,
                          keyType       = "ENSEMBL",
                          OrgDb         = org.Dr.eg.db,
                          ont           = "all",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.1,
                          qvalueCutoff  = 0.1,
                          readable      = TRUE)


egoallSimple <- list()

for (i in c(1:5,7)){
  egoallSimple[[i]] <- simplify(egoall[[i]], cutoff = 0.7)
}

lapply(egoallSimple, dim)

ego_C1 <- data.frame(egoallSimple[[1]])
ego_C2 <- data.frame(egoallSimple[[2]])
ego_C3 <- data.frame(egoallSimple[[3]])
ego_C4 <- data.frame(egoallSimple[[4]])
ego_C5 <- data.frame(egoallSimple[[5]])
ego_C6 <- data.frame(egoallSimple[[6]])
ego_C7 <- data.frame(egoallSimple[[7]])

#original runs saved here, if re-running then the Go signatures change
#write.csv(ego_C1, "ego_C1i.csv", row.names = F)
#write.csv(ego_C2, "ego_C2i.csv", row.names = F)
#write.csv(ego_C3, "ego_C3i.csv", row.names = F)
#write.csv(ego_C4, "ego_C4i.csv", row.names = F)
#write.csv(ego_C5, "ego_C5i.csv", row.names = F)
#write.csv(ego_C6, "ego_C6i.csv", row.names = F)
#write.csv(ego_C7, "ego_C7i.csv", row.names = F)

egoallPlots <- list()

for (i in 1:7){
  egoallPlots[[i]] <-  dotplot(egoall[[i]], showCategory = 10) +
    theme(axis.text.y = element_text(size =12))
}

dotplot(egoall[[4]], showCategory = 10, label_format = function(x) stringr::str_wrap(x, width=30)) +
  theme(axis.text.y = element_text(size =10)) 

grid.arrange(#top = textGrob("\nEC-favoured Clusters GO enrichment\n", gp = gpar(fontface = 'bold', fontsize = 18)),
  #bottom = textGrob("\n", gp = gpar(fontface = 'bold', fontsize = 15)),
  #left = textGrob("\n", gp = gpar(fontface = 'bold', fontsize = 15), rot = 90),
  egoallPlots[[1]], egoallPlots[[2]], egoallPlots[[3]], egoallPlots[[4]], ncol = 4)

grid.arrange(#top = textGrob("\nEC-favoured Clusters GO enrichment\n", gp = gpar(fontface = 'bold', fontsize = 18)),
  #bottom = textGrob("\n", gp = gpar(fontface = 'bold', fontsize = 15)),
  #left = textGrob("\n", gp = gpar(fontface = 'bold', fontsize = 15), rot = 90),
  egoallPlots[[4]], egoallPlots[[5]], egoallPlots[[7]], ncol = 3)

#custom dot plot (used simplified for this):

#dotplot is % on x axis, color for P:
#ranks by Bp then CC then MF so pretty shit system
#focus on BP? then can include "regeneration" from C1, but omit mitochondria terms C4 (v strong)
#focus on p? then no regeneration from C1, no prolif from C2
#go to 10 GO terms? too wordy for a letter
#cherrypick?
dotplot(egoallSimple[[1]], showCategory = 8)
dotplot(egoallSimple[[2]], showCategory = 8)

#too much info for a letter tho, need to get one (color by gene ratio), p val on x axis

ownPlots <- list()

egoallSimple_df <- list(ego_C1, ego_C2, ego_C3,ego_C4,ego_C5,ego_C6,ego_C7)

for (i in c(1:5,7)){
  egoallSimple_df[[i]]$selectHits <- as.numeric(sapply(strsplit(egoallSimple_df[[i]]$GeneRatio, "\\/"), "[[", 1))
  egoallSimple_df[[i]]$select <- as.numeric(sapply(strsplit(egoallSimple_df[[i]]$GeneRatio, "\\/"), "[[", 2))
  egoallSimple_df[[i]]$geneRatio <- egoallSimple_df[[i]]$selectHits/egoallSimple_df[[i]]$select*100
  
  egoallSimple_df[[i]] <- egoallSimple_df[[i]][order(egoallSimple_df[[i]]$Description),]
  
  egoallSimple_df[[i]]$DescriptionII <- stringr::str_wrap(egoallSimple_df[[i]]$Description, width = 20)
  
  egoallSimple_df[[i]]$Description <- factor(egoallSimple_df[[i]]$Description, labels = egoallSimple_df[[i]]$DescriptionII)
  egoallSimple_df[[i]]$Description <- factor(egoallSimple_df[[i]]$Description,
                                             levels = levels(egoallSimple_df[[i]]$Description)[order(egoallSimple_df[[i]]$p.adjust, decreasing = T)])
  
  egoallSimple_df[[i]] <- egoallSimple_df[[i]][order(egoallSimple_df[[i]]$p.adjust, -egoallSimple_df[[i]]$geneRatio),]
  
  egoallSimple_df[[i]] <- egoallSimple_df[[i]][,c(3,10,16)]
  
  colnames(egoallSimple_df[[i]])[3] <- "% of cluster\n with GO term"
}


#can fix the geneRatio color though for all, or just plot the p

for (i in c(1:5,7)){
  ownPlots[[i]] <- ggplot(egoallSimple_df[[i]][1:5,]) + aes(x = -log10(p.adjust), y = Description) +
    geom_bar(stat = "identity", fill = "steelblue", color = "grey60")+
    theme_bw() +
    scale_x_continuous(breaks = seq(0,30,2)) +
    theme(text = element_text(size=24)) +
    xlab("") +
    ylab("")
}

egoallSimple_df[[1]]
levels(egoallSimple_df[[2]]$Description)

ownPlots[[1]]#needs regen in there
ownPlots[[2]]#needs prolif in there, diff scale
ownPlots[[3]]#diff scale
ownPlots[[4]]#diff scale
ownPlots[[5]]#diff scale
ownPlots[[7]]

#for those which are mis-represented by restricting to only 5, do some cherry pick to show interpretation of all GO signature:

#C1, add regen, lose G protein, oxygene trans
selectedTerms <- as.character(egoallSimple_df[[1]][1:5,]$Description)
selectedTerms <- c(selectedTerms[-c(3,5)], "fatty acid ligase\nactivity", "regeneration")

ownPlots[[1]] <- ggplot(
  filter(egoallSimple_df[[1]], Description %in% selectedTerms)) + aes(x = -log10(p.adjust), y = Description) +
  geom_bar(stat = "identity", fill = "steelblue", color = "grey60")+
  scale_x_continuous(breaks = seq(0,30,2)) +
  theme_bw() +
  theme(text = element_text(size=24)) +
  ylab("") +
  xlab("")

#C2, keep ECM organ, add 2x terms, one diff and one cell cycle term
selectedTerms <- as.character(egoallSimple_df[[2]][1:5,]$Description)
selectedTerms <- c(selectedTerms[-c(1:2)], "mitotic cell cycle", "glycosaminoglycan\nbinding")

ownPlots[[2]] <- ggplot(
  filter(egoallSimple_df[[2]], Description %in% selectedTerms)) + aes(x = -log10(p.adjust), y = Description) +
  geom_bar(stat = "identity", fill = "steelblue", color = "grey60")+
  scale_x_continuous(breaks = seq(0,30,10)) +
  theme_bw() +
  theme(text = element_text(size=24)) +
  ylab("") 
xlab("")

ownPlots[[3]] <- ggplot(
  egoallSimple_df[[3]][1:5,]) + aes(x = -log10(p.adjust), y = Description) +
  geom_bar(stat = "identity", fill = "steelblue", color = "grey60")+
  scale_x_continuous(breaks = seq(0,30,1)) +
  theme_bw() +
  theme(text = element_text(size=24)) +
  ylab("") +
  xlab("")

#C4, swap out transporter complex
selectedTerms <- as.character(egoallSimple_df[[4]][1:5,]$Description)
selectedTerms <- c(selectedTerms[-c(5)], "cytochrome complex")

ownPlots[[4]] <- ggplot(
  filter(egoallSimple_df[[4]], Description %in% selectedTerms)) + aes(x = -log10(p.adjust), y = Description) +
  geom_bar(stat = "identity", fill = "steelblue", color = "grey60")+
  theme_bw() +
  scale_x_continuous(breaks = seq(0,30,4)) +
  theme(text = element_text(size=24)) +
  ylab("") +
  xlab("")

#C5, swap out transporter complex, channel activity
selectedTerms <- as.character(egoallSimple_df[[5]][1:5,]$Description)
selectedTerms <- c(selectedTerms[-c(1,4)], "blood vessel\nendothelial cell\ndifferentiation")

ownPlots[[5]] <- ggplot(
  filter(egoallSimple_df[[5]], Description %in% selectedTerms)) + aes(x = -log10(p.adjust), y = Description) +
  geom_bar(stat = "identity", fill = "steelblue", color = "grey60")+
  scale_x_continuous(breaks = seq(0,30,1)) +
  theme_bw() +
  theme(text = element_text(size=24)) +
  ylab("") +
  xlab("")

ownPlots[[1]]#needs regen in there
ownPlots[[2]]#needs prolif in there
ownPlots[[3]]
ownPlots[[4]]
ownPlots[[5]]
ownPlots[[7]]


#facet wrap approach - for reference 
egoallSimple_df_select <- list()

for (i in c(1:7)){
  egoallSimple_df_select[[i]] <- egoallSimple_df[[i]][1:10,]
}

egoallSimple_df_all <- bind_rows(egoallSimple_df_select, .id = "cluster")
egoallSimple_df_all <- filter(egoallSimple_df_all, !is.na(egoallSimple_df_all$p.adjust))

ggplot(egoallSimple_df_all) + aes(x = -log10(p.adjust), y = Description, fill = `% of cluster\n with GO term`) +
  geom_bar(stat = "identity")+
  facet_wrap(~cluster, ncol =7, scales = "free_x") +
  theme_bw() +
  ylab("")

#nope...


#### add in gene annotation, build summary table ####

#import gencode annotation -zebrafish
GRCz11_gtf_table <- read.csv("\\\\cmvm.datastore.ed.ac.uk//cmvm/scs/groups/Baker-lab/Juex/12-17_ljx_annotation_zebrafish.csv",
                             header = T)
miR143_145_cryo_RNAseq_summary <- merge(unique(GRCz11_gtf_table[,c(2:8)]),fpkmDetail, by.x = "gene_id", by.y = "ENSEMBL", all.x = T)

WT_sham_cryo_res$ENSEMBL <- rownames(WT_sham_cryo_res)
Mut_sham_cryo_res$ENSEMBL <- rownames(Mut_sham_cryo_res)

miR143_145_cryo_RNAseq_summary <- merge(miR143_145_cryo_RNAseq_summary, WT_sham_cryo_res[,c(2,6,7)], by.x = "gene_id", by.y = "ENSEMBL", all.x = T)
miR143_145_cryo_RNAseq_summary <- merge(miR143_145_cryo_RNAseq_summary, Mut_sham_cryo_res[,c(2,6,7)], by.x = "gene_id", by.y = "ENSEMBL", all.x = T)
colnames(miR143_145_cryo_RNAseq_summary)
colnames(miR143_145_cryo_RNAseq_summary)[c(33:36)] <- c("LFC_WT_injury", "Padj_WT_injury", "LFC_Mut_injury", "Padj_Mut_injury")

#add cluster info
miR143_145_cryo_RNAseq_summary$DEG_Cluster <- "None"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$gene_id %in% Cluster1] <- "BiggerInMutant-cxcl8/wnt4/slit3"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$gene_id %in% Cluster2] <- "SimilarlyUp-mki67/col1a1/pcna/myh7"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$gene_id %in% Cluster3] <- "OpposingUpMutant-tpm4a/lrrc10"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$gene_id %in% Cluster4] <- "WeakerDownMut-vegfaa"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$gene_id %in% Cluster5] <- "WeakMutantNoChange-tnnt2a/synpo2lb/ryr2b"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$gene_id %in% Cluster6] <- "OpposingDownMutant"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$gene_id %in% Cluster7] <- "OnlyDownMutant-actn2b/myoz2"


miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$DEG_Cluster == "BiggerInMutant-cxcl8/wnt4/slit3"] <- "C1-BiggerInMutant-cxcl8/wnt4/slit3"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$DEG_Cluster == "SimilarlyUp-mki67/col1a1/pcna/myh7"] <- "C2-SimilarlyUp-mki67/col1a1/pcna/myh7"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$DEG_Cluster == "OpposingUpMutant-tpm4a/lrrc10"] <- "C3-OpposingUpMutant-tpm4a/lrrc10"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$DEG_Cluster == "WeakerDownMut-vegfaa"] <- "C4-WeakerDownMut-vegfaa"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$DEG_Cluster == "WeakMutantNoChange-tnnt2a/synpo2lb/ryr2b"] <- "C5-WeakMutantNoChange-tnnt2a/synpo2lb/ryr2b"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$DEG_Cluster == "OpposingDownMutant"] <- "C6-OpposingDownMutant"
miR143_145_cryo_RNAseq_summary$DEG_Cluster[miR143_145_cryo_RNAseq_summary$DEG_Cluster == "OnlyDownMutant-actn2b/myoz2"] <- "C7-OnlyDownMutant-actn2b/myoz2"

#write.csv(miR143_145_cryo_RNAseq_summary, "miR143_145_cryo_RNAseq_summary.csv", row.names = F)

fpkm_allG <- filter(miR143_145_cryo_RNAseq_summary, fpkm_max_treatment >1)
fpkm_allGDE <- filter(fpkm_allG, (LogFC >log2(1.5) | LogFC < -log2(1.5)) & padj <0.05)
table(fpkm_allGDE$gene_id)


#plotting some initial targets - volcano for PCGs:


#### complete summary table for GEO - where do genotype:injury DEGs appear in clustering ####

miR143_145_cryo_RNAseq_summary <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/users/mbennet5/Mukesh CARMN/mirMutantAnalysis/miR143_145_cryo_RNAseq_summary.csv")

colnames(miR143_145_cryo_RNAseq_summary)

miR143_145_cryo_DE <- filter(miR143_145_cryo_RNAseq_summary, !DEG_Cluster == "None")

miR143_145_cryo_DE_GI <- filter(miR143_145_cryo_RNAseq_summary, !DEG_Cluster == "None", padj_GenotypeInjury < 0.05)[,-c(2:5,8:23,25,27,29,31)]
miR143_145_cryo_DE_GI <- miR143_145_cryo_DE_GI[,c(1:2,13,3:12,14:15)]

#for mukesh to pick out genes of interest for text/discussion
#write.csv(miR143_145_cryo_DE_GI, "miR143_145_cryo_DE_GI.csv", row.names = F)

table(miR143_145_cryo_DE$DEG_Cluster)
table(filter(miR143_145_cryo_DE, gene_id %in% rownames(filter(GenotypeInjuryInteraction, padj < 0.05)))$DEG_Cluster)
table(filter(miR143_145_cryo_DE, gene_id %in% rownames(filter(GenotypeInjuryInteraction, padj < 0.05)))$DEG_Cluster)/table(miR143_145_cryo_DE$DEG_Cluster) * 100
table(filter(miR143_145_cryo_DE, gene_id %in% rownames(filter(GenotypeInjuryInteraction, padj < 0.05)))$DEG_Cluster)/98*100

#simplest: % found in each cluster:
plot_int_cluster <- as.data.frame(table(filter(miR143_145_cryo_DE, gene_id %in% rownames(filter(GenotypeInjuryInteraction, padj < 0.05)))$DEG_Cluster)/table(miR143_145_cryo_DE$DEG_Cluster) * 100)
plot_int_cluster$cluster <- sapply(strsplit(as.character(plot_int_cluster$Var1), "-"), "[[", 1)

ggplot(plot_int_cluster) + aes(x = cluster, y = Freq) + 
  geom_bar(stat = "identity", color = "grey60", fill = "grey90") +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  xlab("") +
  ylab("% Genotype-specific\nCryo-responsive genes")

#add in column of genotype-specific p:
GenotypeInjuryInteraction$EnsID <- rownames(GenotypeInjuryInteraction)
trial <- merge(miR143_145_cryo_RNAseq_summary, GenotypeInjuryInteraction, by.x = "gene_id", by.y = "EnsID", all.x = T)

#write.csv(miR143_145_cryo_RNAseq_summary, "miR143_145_cryo_RNAseq_summary_2026.csv", row.names = F)


#
#### Q4, plausible mechanisms, TargetScan targets ####

#available at TargetScan online (accessed Q3 2024)

#miR-143 - (3p strand)
miR143_targets <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/users/mbennet5/Mukesh CARMN/mirMutantAnalysis/targetscan_6_2_zf_mir-143.csv")
#miR-145 - (5p strand)
miR145_targets <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/users/mbennet5/Mukesh CARMN/mirMutantAnalysis/targetscan_6_2_zf_mir-145.csv")

#match via ENSDARG:
miR143_targets$gene_id <- gsub("\\.[0-9]*", "", miR143_targets$Representative.3..UTR)
miR145_targets$gene_id <- gsub("\\.[0-9]*", "", miR145_targets$Representative.3..UTR)

#targets which are different between genotypes in cryo are the key set to examine
cryo_WT_mut_DE_anno$EnsID_merge <- gsub("\\.[0-9]*", "", cryo_WT_mut_DE_anno$EnsID)

DEGs_forPrediction <- cryo_WT_mut_DE_anno

#cbx4 is a match on gene name only, ENSDARG00000070807.1 in TargetScan rep UTR column and ENSDARG00000099441 in our RNAseq annotation
#ENSDARG00000070807 no longer found in ZFIN, it may be historical ID for cbx4

cryo_WT_mut_DE_targets <- filter(DEGs_forPrediction, EnsID_merge %in% c(miR143_targets$gene_id, miR145_targets$gene_id) |
                                   gene_name %in% c(miR143_targets$Target.gene, miR145_targets$Target.gene), log2FoldChange >0)

#save a table:
trial <- merge(cryo_WT_mut_DE_targets, miR143_targets[,c(1,9)], by.x = "gene_name", by.y = "Target.gene", all.x = T)
trial <- merge(trial, miR143_targets[,c(11,9)], by.x = "EnsID", by.y = "gene_id", all.x = T)
trial$Total.context..score.x[is.na(trial$Total.context..score.x)] <- trial$Total.context..score.y[is.na(trial$Total.context..score.x)]
cryo_WT_mut_DE_targets <- trial[,-16]
colnames(cryo_WT_mut_DE_targets)[15] <- "miR-143_TargetScan_TotalContextScore"

trial <- merge(cryo_WT_mut_DE_targets, miR145_targets[,c(1,9)], by.x = "gene_name", by.y = "Target.gene", all.x = T)
trial <- merge(trial, miR145_targets[,c(11,9)], by.x = "EnsID", by.y = "gene_id", all.x = T)
trial$Total.context..score.x[is.na(trial$Total.context..score.x)] <- trial$Total.context..score.y[is.na(trial$Total.context..score.x)]
cryo_WT_mut_DE_targets <- trial[,-17]
colnames(cryo_WT_mut_DE_targets)[16] <- "miR-145_TargetScan_TotalContextScore"

#save
#write.csv(cryo_WT_mut_DE_targets, "MiRNA_143_145_Targets_CryoInjury_2026.csv", row.names = F)

#key info:
cryo_WT_mut_DE_targets <- cryo_WT_mut_DE_targets[,c(1,2,4,8,15,16)]
View(cryo_WT_mut_DE_targets)


#add in cluster/genotype:injury info where exists:
cryo_WT_mut_DE_targets <- merge(cryo_WT_mut_DE_targets, miR143_145_cryo_DE[,c(1,37,39)], by.x = "EnsID", by.y = "gene_id", all.x = T)

cryo_WT_mut_DE_nonTargets <- filter(DEGs_forPrediction, EnsID_merge %in% c(miR143_targets$gene_id, miR145_targets$gene_id) |
                                      gene_name %in% c(miR143_targets$Target.gene, miR145_targets$Target.gene), log2FoldChange <0)
#is there a confidence thresh at which targets begin to be enriched amongst DEGs? should be genotype-associated DEGs during cryo ideally

#build fisher's test to confirm targets in expected gene sets and add confidence to predictions
dim(filter(DEGs_forPrediction, log2FoldChange >0))
dim(filter(DEGs_forPrediction, log2FoldChange <0))
dim(DEGs_forPrediction)

#usual more stringent DEG threshold
36/116 #31.0% mutant-induced are in targetscan (n.b. 15 previously reported)
26/123 #21.1% mutant-repressed are in targetscan

a <- 36
b <- 116
c <- 36+26
d <- 239

#some decent bias to mut-induced via Fisher's - especially for a low gene pool
fisher.test(data.frame("Induced" = c(a,b-a),
                       "Repressed" = c(c-a,d-c-(b-a)), row.names = c("TargetScan", "Not")), alternative = "greater")

#
#### heatmap of predicted targets ####

MiRNA_143_145_Targets_CryoInjury <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/users/mbennet5/Mukesh CARMN/mirMutantAnalysis/MiRNA_143_145_Targets_CryoInjury_2026.csv")

#miR-143:
miR143_targets <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/users/mbennet5/Mukesh CARMN/mirMutantAnalysis/targetscan_6_2_zf_mir-143.csv")
miR145_targets <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/users/mbennet5/Mukesh CARMN/mirMutantAnalysis/targetscan_6_2_zf_mir-145.csv")

#match via ENSDARG:
miR143_targets$gene_id <- gsub("\\.[0-9]*", "", miR143_targets$Representative.3..UTR)
miR145_targets$gene_id <- gsub("\\.[0-9]*", "", miR145_targets$Representative.3..UTR)

#annotate a heatmap with miR-143, miR-145 target confidence, FPKM, annotation in heatmap, GO/KEGG...?

#any highlighted details from Kesh too - version to be finalised

MiRNA_143_145_Targets_CryoInjury_fpkm <- filter(fpkmDetail, ENSEMBL %in% MiRNA_143_145_Targets_CryoInjury$EnsID)

mat <- MiRNA_143_145_Targets_CryoInjury_fpkm[,2:17] #36 genes
rownames(mat) <- miR143_145_cryo_RNAseq_summary$gene_name[match(MiRNA_143_145_Targets_CryoInjury_fpkm[,1], miR143_145_cryo_RNAseq_summary$gene_id)]
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
mat <- t(apply(mat, 1, cal_z_score))

myColor <- colorRampPalette(c("steelblue", "white", "red"))(50)
myBreaks <- c(seq(min(mat), 0, 
                  length.out=ceiling(50/2)), 
              seq(max(mat)/50, 
                  max(mat), 
                  length.out=floor(50/2)))

data_cols<-data

data_cols$genotype <- "WT"
data_cols$genotype[grepl("Mut", data_cols$condition)] <- "miR_ko"
data_cols$cryoInjury <- "Sham"
data_cols$cryoInjury[grepl("cryo", data_cols$condition)] <- "Cryo-injury"

rownames(data_cols) == colnames(mat)

#data_colsHeat$Hours <- factor(data_colsHeat$Hours, levels(data_colsHeat$Hours)[c(1,3,4,2)])

my_colour_anno = list(
  Injury = c(`Cryo-injury` = "blueviolet", Sham = "chartreuse4"),
  Genotype = c(WT = "coral1", `miR-zdKO` = "steelblue"),
  miR.143_Q = c("TargetScan Hit" = "grey90", "Strong TargetScan Hit(<-0.1)" = "grey30"),
  miR.145_Q = c("TargetScan Hit" = "grey90", "Strong TargetScan Hit(<-0.1)" = "grey30")
)


data_rows <- MiRNA_143_145_Targets_CryoInjury[,c(2,14:16)]
data_rows <- merge(data_rows, miR143_145_cryo_RNAseq_summary[,c(1,37)], by.x = "EnsID_merge", by.y = "gene_id") 
rownames(data_rows) <- data_rows$gene_name

#strength colours (grey-scale), anything >-0.1 in sep category
data_rows$miR.143_Q[grepl(">", data_rows$miR.143_TargetScan_TotalContextScore)] <- "TargetScan Hit"
data_rows$miR.143_Q[as.numeric(data_rows$miR.143_TargetScan_TotalContextScore) > -0.1] <- "TargetScan Hit"
data_rows$miR.143_Q[as.numeric(data_rows$miR.143_TargetScan_TotalContextScore) <= -0.1] <- "Strong TargetScan Hit(<-0.1)"
#strength colours (grey-scale), anything >-0.1 in sep category
data_rows$miR.145_Q[grepl(">", data_rows$miR.145_TargetScan_TotalContextScore)] <- "TargetScan Hit"
data_rows$miR.145_Q[as.numeric(data_rows$miR.145_TargetScan_TotalContextScore) > -0.1] <- "TargetScan Hit"
data_rows$miR.145_Q[as.numeric(data_rows$miR.145_TargetScan_TotalContextScore) <= -0.1] <- "Strong TargetScan Hit(<-0.1)"


rownames(mat) <- tolower(rownames(mat))

library(pheatmap)
p <- pheatmap(mat[rownames(mat) %in% miR143_145_cryo_DE_GI$gene_name,],
              annotation_colors = my_colour_anno,gaps_col = c(8),
              annotation_col = data.frame("Injury" = data_cols[,9], row.names = rownames(data_cols)),
              annotation_row = data_rows[,6:7],
              show_colnames = F, 
              show_rownames = T, 
              cluster_cols = F,
              cluster_rows = T,annotation_legend = T,fontsize_row = 20,
              #cutree_rows = 2,
              treeheight_col = 0, 
              treeheight_row = 20,
              color = myColor, 
              breaks = myBreaks,
              border_color = NA)

#### Q5, plausible mechanisms, local effects ####

#show carmn and neighbours in the WT v mut volcano plots
sham_WT_mut_DE$EnsID <- rownames(sham_WT_mut_DE)

#miR143_145_cryo_RNAseq_summary <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/users/mbennet5/Mukesh CARMN/mirMutantAnalysis/miR143_145_cryo_RNAseq_summary.csv")
sham_WT_mut_DE_anno <- merge(sham_WT_mut_DE, miR143_145_cryo_RNAseq_summary[,1:7], by.x = "EnsID", by.y = "gene_id")

#genes near carmn:
filter(sham_WT_mut_DE_anno, gene_name == "BX088707.3")
carmn_locus <- filter(sham_WT_mut_DE_anno, seqnames == 14, start > 38715673 - 5000000, start < 38744637 + 5000000)$gene_name

sham_WT_mut_DE_anno$DE <- "No"
sham_WT_mut_DE_anno$DE[rownames(sham_WT_mut_DE_anno) %in% rownames(sham_WT_mut_DE_anno)] <- "DE"

sham_WT_mut_DE_anno$Notable <- NA
sham_WT_mut_DE_anno$Notable[sham_WT_mut_DE_anno$gene_name %in% carmn_locus] <- "carmn_locus"

sham_WT_mut_DE_anno$gene_name[sham_WT_mut_DE_anno$gene_name == "BX088707.3"] <- "carmn"

ggplot(sham_WT_mut_DE_anno) + aes(x = log2FoldChange, y = -log10(padj), color = DE, label = gene_name) +
  geom_point(alpha = 0.5, size = 3) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") + 
  ggrepel::geom_label_repel(data = filter(sham_WT_mut_DE_anno, !is.na(Notable)), 
                            nudge_x = -2, , nudge_y = 10, force = 600, size = 6) + 
  #coord_cartesian(xlim = c(-3.9,3.9), ylim = c(0,20)) +
  scale_color_manual(values = c("No" = "grey60", "DE" = "steelblue")) +
  theme_minimal() +
  theme(text = element_text(size = 20)) +
  xlab("LFC miR-dKo vs. WT\n(Sham)") +
  ylab("-log10(p)")


#show carmn and neighbours in the WT v mut volcano plots
cryo_WT_mut_DE$EnsID <- rownames(cryo_WT_mut_DE)

#miR143_145_cryo_RNAseq_summary <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/users/mbennet5/Mukesh CARMN/mirMutantAnalysis/miR143_145_cryo_RNAseq_summary.csv")
cryo_WT_mut_DE_anno <- merge(cryo_WT_mut_DE, miR143_145_cryo_RNAseq_summary[,1:7], by.x = "EnsID", by.y = "gene_id")

#genes near carmn:
filter(cryo_WT_mut_DE_anno, gene_name == "BX088707.3")
carmn_locus <- filter(cryo_WT_mut_DE_anno, seqnames == 14, start > 38715673 - 5000000, start < 38744637 + 5000000)$gene_name

cryo_WT_mut_DE_anno$DE <- "No"
cryo_WT_mut_DE_anno$DE[rownames(cryo_WT_mut_DE_anno) %in% rownames(cryo_WT_mut_DE_anno)] <- "DE"

cryo_WT_mut_DE_anno$Notable <- NA
cryo_WT_mut_DE_anno$Notable[cryo_WT_mut_DE_anno$gene_name %in% carmn_locus] <- "carmn_locus"

cryo_WT_mut_DE_anno$gene_name[cryo_WT_mut_DE_anno$gene_name == "BX088707.3"] <- "carmn"

ggplot(cryo_WT_mut_DE_anno) + aes(x = log2FoldChange, y = -log10(padj), color = DE, label = gene_name) +
  geom_point(alpha = 0.5, size = 2.5) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") + 
  ggrepel::geom_label_repel(data = filter(cryo_WT_mut_DE_anno, !is.na(Notable)), 
                            nudge_x = -2, , nudge_y = 20, force = 100, size = 6) + 
  #coord_cartesian(xlim = c(-3.9,3.9), ylim = c(0,20)) +
  scale_color_manual(values = c("No" = "grey60", "DE" = "steelblue")) +
  theme_minimal() +
  theme(text = element_text(size = 20)) +
  xlab("LFC miR-dKo vs. WT\n(Cryo)") +
  ylab("-log10(p)")

