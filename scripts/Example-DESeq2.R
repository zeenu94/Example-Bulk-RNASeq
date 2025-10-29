library(dplyr)
library(DESeq2)

key <- read.csv("rosmap.txt",sep="\t")
colnames(key)[6]<- "SampleID"


assay <- read.csv("assay_rnaSeq_metadata.csv", sep = ",")
biospecimen <- read.csv("biospecimen_metadata.csv",sep=",")
clinical <- read.csv("clinical.csv")

meta_rna <- biospecimen %>%
  inner_join(clinical, by = "individualID") %>%
  inner_join(assay, by = "specimenID")



# Clean and filter metadata\
meta_rna <- meta_rna_spec
meta_rna <- meta_rna %>% filter(AD == 1 | AD==0)
meta_rna$pmi[meta_rna$pmi == "missing or unknown"] <- NA
meta_rna$pmi <- as.numeric(meta_rna$pmi)
meta_rna$age_death[meta_rna$age_death == "90+"] <- 90
meta_rna$age_death <- as.numeric(meta_rna$age_death)

meta_rna$tissue[meta_rna$tissue=="dorsolateral prefrontal cortex"]="DLPFC"
meta_rna$tissue[meta_rna$tissue=="Head of caudate nucleus"]="CN"
meta_rna$tissue[meta_rna$tissue=="posterior cingulate cortex"]="PCC"

meta_rna$libraryBatch[is.na(meta_rna$libraryBatch)] <- names(which.max(table(meta_rna$libraryBatch)))


library(data.table)

# DLPFC

dlp <- fread("dlpfc_readsCount.gct.gz")
mat_dlp <- as.matrix(dlp)
rownames(mat_dlp) <- dlp$gene_ID
dim(mat_dlp)
storage.mode(mat_dlp) <- "integer"

common_samples <- intersect(colnames(mat_dlp), meta_rna$specimenID)
bulk_expr <- mat_dlp[, common_samples]
meta_rna <- meta_rna[meta_rna$specimenID %in% common_samples, ]
meta_rna <- meta_rna[match(colnames(bulk_expr), meta_rna$specimenID), ]
rownames(meta_rna) <- meta_rna$specimenID

meta_rna <- meta_rna[, c("specimenID", "AD", "msex", "pmi", "age_death","RIN",  "tissue","libraryBatch")]
meta_rna <- meta_rna[complete.cases(meta_rna), ]
bulk_expr <- bulk_expr[,colnames(bulk_expr) %in% meta_rna$specimenID]


meta_rna$msex        <- factor(meta_rna$msex)

meta_rna$libraryBatch<- factor(meta_rna$libraryBatch)
meta_rna$AD <- factor(meta_rna$AD,levels = c(0,1))

dds <- DESeqDataSetFromMatrix(
  countData = bulk_expr_subset,
  colData   = meta_rna,
  design    = ~ AD +  msex + age_death + RIN  + libraryBatch
)

# Filter very-low counts (DESeq2 suggestion)
keep_genes <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep_genes, ]
dds <- DESeq(dds)

## 7) TL effect: Long vs Short
res <- results(dds, contrast = c("AD","1","0"))
res$gene <- rownames(res)
res <- res[order(res$padj), ]
res
write.table(res,"results.csv", col.names = T,row.names = T, sep = "\t",quote = F)

