### TCGA Transcriptomics-Proteomics Integration

#############################################
# Step 1 - Set working directory
#############################################
setwd("C:/Users/Tony Lin/Desktop/Wiita_lab/Projects/Proteomics_project/Phosphoproteomics_MMCL")


#############################################
# Step 2a - Read in data files (phospho)
#############################################
#####
#1. Read phosphosites data
b.phospho = read.table(file = 'TCGA-CPTAC_BRCA/data-phospho/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv',
                  sep = '\t', header = TRUE, colClasses = "character", row.names = 1)

o.phospho = read.table(file = 'TCGA-CPTAC_OV/data-phospho/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv',
                       sep = '\t', header = TRUE, colClasses = "character", row.names = 1)

#2. Trim header
names(b.phospho) = sub("\\.Log\\.Ratio$", "", names(b.phospho))
names(o.phospho) = sub("\\.01A\\.Log\\.Ratio$", "", names(o.phospho))
#####


#############################################
# Step 2b - Read in data files (RNAseq)
#############################################
#####
#1. Collect file names
file_names_RNA = function(file_name, data_type) {
  # file_names = a character string detailing name of file
  # data_type = a character string detailing type of RNA file names to extract (i.e. "counts", "FPKM")
  mapID = read.table(file = file_name, sep = "\t", header = TRUE, colClasses = "character")
  fileID = data.frame(file_name = mapID$file_name, label = mapID$cases.0.submitter_id)
  fileID = fileID[grepl(data_type, fileID$file_name), ]   #can be updated to look at FPKM
  fileID$file_name = sub("\\.gz$", "", fileID$file_name)
  return(fileID)   #Output a data frame containing a label and a file_name column
}
b.fileID = file_names_RNA("TCGA-CPTAC_BRCA/GDC_extract/TCGA-BRCA_submitterID_fileID_mapping.tsv",
                          data_type = "counts")
o.fileID = file_names_RNA("TCGA-CPTAC_OV/GDC_extract/TCGA-OV_submitterID_fileID_mapping.tsv",
                          data_type = "counts")

#2. Read RNAseq files from fileID
read_RNA = function(path, fileID) {
  # path = a character string describing the directory of uncompressed RNAseq files
  # fileID = output from file_names_RNA function
  rna = apply(fileID, 1,
               function(x) {
                 df = read.table(paste0(path, x[1]),
                                 sep = "\t", colClasses = c("character", "numeric"))
                 names(df) = c("Gene", x[2])   # Labeling column names
                 return(df)
               })
  rna = Reduce(function(df1, df2) merge(df1, df2, by = "Gene", all = TRUE), rna)
  return(rna)   #returns a data frame with each column as a sample and each row as a gene
}
b.rna = read_RNA("TCGA-CPTAC_BRCA/data-RNAseq-uncompressed/", b.fileID)
o.rna = read_RNA("TCGA-CPTAC_OV/data-RNAseq-uncompressed/", o.fileID)
#####
#save.image("raw_data.RData") #saves workspace


#############################################
# Step 3 - Standardize Gene Identifier (RNAseq)
#############################################
#load("raw_data.RData")
#####
# Remove rows with invalid gene names
b.rna = b.rna[grepl("^ENSG", b.rna$Gene), ]
o.rna = o.rna[grepl("^ENSG", o.rna$Gene), ]

# Convert ENSEMBL gene ID to HGNC gene symbol
require(curl)
require(rjson)
require(gProfileR)
ENSG = gsub("\\..*", "", b.rna$Gene)   #Breast and ovarian data carry same list of identifiers
ENSG = split(ENSG, ceiling(seq_along(ENSG)/5000))   #Query 50 Gene IDs at a time
url.rna = lapply(ENSG, function(x) gconvert(x, target = "HGNC", filter_na = FALSE,
                                            mthreshold = 1)[c("alias", "name")])
##bioDBnet mapping method (faster but lot of missing gene symbols)
#url.rna = lapply(ENSG,
#                 function(x) {
#                   URL = paste0("https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi?method=db2db&input=ensemblgeneid&inputValues=",
#                          paste(x, collapse = ","),
#                          "&outputs=genesymbol&taxonId=9606&format=col")
#                   curl_download(URL, "ENSG2GeneSymbol.json")
#                   return(as.data.frame(fromJSON(file = "ENSG2GeneSymbol.json"),
#                                        stringsAsFactors = FALSE))
#                 })   #Loop through ENSG to extract 500 gene symbols at a time from bioDBnet
url.rna = Reduce(rbind, url.rna)   #Collapse list into a single data frame

# Filter rows for valid gene symbols
b.rna$Symbol = url.rna[[2]]
o.rna$Symbol = url.rna[[2]]
b.rna = b.rna[b.rna$Symbol != "N/A", ]
o.rna = o.rna[o.rna$Symbol != "N/A", ]
#####


#############################################
# Step 4 - Standardize Gene Identifier (phospho)
#############################################
#####
# Convert RefSeq protein ID to gene symbol for breast data
require(curl)
require(rjson)
b.RefSeq = gsub("\\..*", "", row.names(b.phospho))
b.RefSeq = split(b.RefSeq, ceiling(seq_along(b.RefSeq)/500))   #Query 500 Gene IDs at a time
b.url.phospho = lapply(b.RefSeq,
                 function(x) {
                   URL = paste0("https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi?method=db2db&input=refseqproteinaccession&inputValues=",
                                paste(x, collapse = ","),
                                "&outputs=genesymbol&taxonId=9606&format=col")
                   curl_download(URL, "RefSeq2GeneSymbol.json")
                   return(as.data.frame(fromJSON(file = "RefSeq2GeneSymbol.json"),
                                        stringsAsFactors = FALSE))
                 })   #Loop through b.RefSeq to extract 500 gene symbols at a time
b.url.phospho = Reduce(rbind, b.url.phospho)   #Collapse list into a single data frame
b.phospho$Symbol = b.url.phospho[[2]]
b.phospho$Gene = gsub(":.*", "", row.names(b.phospho))
b.phospho$Site = toupper(paste0(b.phospho$Symbol, "_", gsub(".*:", "", row.names(b.phospho))))

# Convert RefSeq protein ID to gene symbol for ovarian data
o.RefSeq = gsub("\\..*", "", row.names(o.phospho))
o.RefSeq = split(o.RefSeq, ceiling(seq_along(o.RefSeq)/500))   #Query 500 Gene IDs at a time
o.url.phospho = lapply(o.RefSeq,
                       function(x) {
                         URL = paste0("https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi?method=db2db&input=refseqproteinaccession&inputValues=",
                                      paste(x, collapse = ","),
                                      "&outputs=genesymbol&taxonId=9606&format=col")
                         curl_download(URL, "RefSeq2GeneSymbol.json")
                         return(as.data.frame(fromJSON(file = "RefSeq2GeneSymbol.json"),
                                              stringsAsFactors = FALSE))
                       })   #Loop through b.RefSeq to extract 500 gene symbols at a time
o.url.phospho = Reduce(rbind, o.url.phospho)   #Collapse list into a single data frame
o.phospho$Symbol = o.url.phospho[[2]]
o.phospho$Gene = gsub(":.*", "", row.names(o.phospho))
o.phospho$Site = toupper(paste0(o.phospho$Symbol, "_", gsub(".*:", "", row.names(o.phospho))))
#####


#############################################
# Step 5 - Standardize Sample Name
#############################################
#####
# Standardize BRCA phospho labels
names(b.phospho)[1:3] = c("normal1", "normal2", "normal3")
names(b.phospho) = sub("\\.01A", "", names(b.phospho))
names(b.phospho) = gsub("\\.", "-", names(b.phospho))

# Standardize BRCA RNA labels
names(b.rna) = sub("^TCGA-", "", names(b.rna))

# Standardize OV phospho labels
names(o.phospho) = sub("^X", "", names(o.phospho))
names(o.phospho) = sub("\\.", "-", names(o.phospho))

# Standardize OV RNA labels
names(o.rna) = sub("^TCGA-", "", names(o.rna))
#####
#save.image("clean_data.RData") #saves workspace


#############################################
# Step 6 - Convert RNA read counts to CPM
#############################################
#load("clean_data.RData") #saves workspace
#####
require(edgeR)
require(RColorBrewer)

## Convert data frame to DGEList data class
b.counts = DGEList(counts = b.rna[, 2:106], genes = b.rna[c("Gene", "Symbol")])
o.counts = DGEList(counts = o.rna[, 2:51], genes = o.rna[c("Gene", "Symbol")])

## Filtering weakly expressed genes or non-informative features
b.keep <- rowSums(cpm(b.counts) > 1) >= (0.2 * 105)   # 20% detection for retention
b.counts <- b.counts[b.keep, , keep.lib.sizes = FALSE]
o.keep <- rowSums(cpm(o.counts) > 1) >= (0.2 * 50)   # 20% detection for retention
o.counts <- o.counts[o.keep, , keep.lib.sizes = FALSE]


## Calculate normalization factors to account for compositional differences among libraries
b.counts <- calcNormFactors(b.counts)
o.counts <- calcNormFactors(o.counts)

## Compute normalized log2-transformed CPM
b.logcpm = cpm(b.counts, normalized.lib.sizes = TRUE, log = TRUE)
rownames(b.logcpm) = b.counts$genes$Symbol
o.logcpm = cpm(o.counts, normalized.lib.sizes = TRUE, log = TRUE)
rownames(o.logcpm) = o.counts$genes$Symbol

## Visualize relationsihps between cell lines using multidimentional scaling (MDS)
b.clinical = read.csv(file = "TCGA-CPTAC_BRCA/BRCA_TCGA_tumor_clinical_data.csv")
b.clinical$Complete.TCGA.ID = sub("^TCGA-", "", b.clinical$Complete.TCGA.ID)
b.clinical = b.clinical[match(colnames(b.counts$counts), 
                              b.clinical$Complete.TCGA.ID), ]   #Row order corresponds to b.counts

# Visualize RNAseq breast data by ER status
##CHECK OUT DOCUMENTATION FOR DETAILS ON BCV vs LOG-FC METHOD
##https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/plotMDS.DGEList
brewer = brewer.pal(4, "Dark2")
cairo_ps("ER_status_bcv.eps")
plotMDS(b.counts, col = brewer[b.clinical$ER.Status], method = "bcv", 
        pch = 20, main = "ER Status")
legend("bottomleft", legend = c("Negative", "Indeterminate", "Positive"), pch = 20, 
       col = c("#D95F02", "#1B9E77", "#7570B3"))
dev.off()

# Visualize RNAseq breast data by HER2 status
cairo_ps("HER2_status_bvc.eps")
plotMDS(b.counts, col = brewer[b.clinical$HER2.Final.Status], method = "bcv",
        pch = 20, main = "HER2 Status")
legend("bottomleft", legend = c("Negative", "Equivocal", "Positive"), pch = 20, 
       col = c("#D95F02", "#7570B3", "#1B9E77"))
dev.off()

# Visualize RNAseq breast data by PR status
cairo_ps("PR_status_bvc.eps")
plotMDS(b.counts, col = brewer[b.clinical$PR.Status], method = "bcv",
        pch = 20, main = "PR Status")
legend("bottomleft", legend = c("Negative", "Positive"), pch = 20, 
       col = c("#1B9E77", "#D95F02"))
dev.off()

# Visualize RNAseq breast data by PAM50
cairo_ps("PAM50_bvc.eps")
plotMDS(b.counts, col = brewer[b.clinical$PAM50.mRNA], method = "bcv",
        pch = 20, main = "PAM50 Assignment")
legend("bottomleft", legend = c("Basal-like", "Luminal A", "Luminal B", "HER2-enriched"), 
       pch = 20, col = c("#1B9E77", "#7570B3", "#E7298A", "#D95F02"))
dev.off()

# Read in PAM50 genes for breast cancer subtyping
PAM50 = read.table("pam50_annotation.txt", header = TRUE, colClasses = "character")
PAM50$GeneName[PAM50$GeneName == "ORC6L"] = "ORC6"   #Update to HUGO Gene Symbol

# Create heatmap of PAM50 gene expression
require(gplots)
require(RColorBrewer)
cairo_ps("PAM50_heatmap.eps")
col.pan <- colorpanel(100, "blue", "white", "red")
temp.data = b.logcpm
colnames(temp.data) = b.clinical$PAM50.mRNA
temp.data = temp.data[rownames(temp.data) %in% PAM50$GeneName, ]
colColors = brewer.pal(4, "Dark2")[b.clinical$PAM50.mRNA]   #Color columns by subtype
heatmap.2(temp.data, col = col.pan, Rowv = TRUE, Colv = TRUE, scale = "row",
          trace = "none", dendrogram = "column", cexRow=1, cexCol=0.5,
          margin=c(4,5), lhei = c(1.25,8), lwid = c(1.5,6), ColSideColors = colColors,
          density.info = "histogram", key.par = list(cex=0.5), 
          key.title = "Distribution of Z-scores", key.xlab = "", key.ylab = "Counts")
rm(temp.data, colColors)
dev.off()
#####
#save.image("checkpoint_1.RData") #saves workspace


#############################################
# Step 7 - Filter and impute phospho data
#############################################
#load("checkpoint_1.RData")
#####
# Reorganize phospho data
b.phospho = list(values = b.phospho[1:111],
                 index = b.phospho[112:114])
o.phospho = list(values = o.phospho[1:69],
                 index = o.phospho[70:72])

# Filter out missing phospho data (require quantification in at least 50% of sample)
b.keep = apply(b.phospho$values, 1, function(x) 1 - sum(is.na(x))/length(x))
b.keep = b.keep > 0.5
b.phospho = lapply(b.phospho, function(x) x[b.keep, ])
o.keep = apply(o.phospho$values, 1, function(x) 1 - sum(is.na(x))/length(x))
o.keep = o.keep > 0.5
o.phospho = lapply(o.phospho, function(x) x[o.keep, ])
rm(b.keep, o.keep)

# Impute missing values using k-nearest neighbor (k-NN) imputation
require(pamr)
b.phospho$values = as.data.frame(pamr.knnimpute(list(x = as.matrix(b.phospho$values)))[["x"]])
o.phospho$values = as.data.frame(pamr.knnimpute(list(x = as.matrix(o.phospho$values)))[["x"]])

# Plot distribution of BRCA phospho data
require(ggplot2)
require(tidyr)
cairo_ps("BRCA_phospho_imputed.eps")
b.plot = gather(b.phospho$values, "sample", "logRatio")
ggplot(data = b.plot, aes(b.plot$logRatio)) + 
  geom_histogram() +
  facet_wrap(~ sample) +
  labs(x = "iTRAQ log( fold change )", y = "Count")
rm(b.plot)
dev.off()

# Plot distribution of OV phospho data
require(ggplot2)
require(tidyr)
cairo_ps("OV_phospho_imputed.eps")
o.plot = gather(o.phospho$values, "sample", "logRatio")
ggplot(data = o.plot, aes(o.plot$logRatio)) + 
  geom_histogram() +
  facet_wrap(~ sample) +
  labs(x = "iTRAQ log( fold change )", y = "Count")
rm(o.plot)
dev.off()
#####
#save.image("checkpoint_2.RData")   #Saves workspace


#############################################
# Step 8 - Estimate kinase activities (Kinase Set Enrichment Analysis)
#############################################
#load("checkpoint_2.RData")
#####
# Read in PhosphoSitePlus Kinase-Substrate table
pSitePlus = read.table("PhosphoSitePlus_Kinase_Substrate_Dataset", header = TRUE,
                       colClasses = "character", skip = 3, fill = TRUE, sep = "\t")
pSitePlus = pSitePlus[pSitePlus$KIN_ORGANISM == "human" &
                        pSitePlus$SUB_ORGANISM == "human", ]   #Filter for interactions in human

# Fill in missing substrate gene names in pSitePlus
require(gProfileR)
missing = unique(pSitePlus[pSitePlus$SUB_GENE == "", "SUB_ACC_ID"])
missing = data.frame(genes = missing,
                     cleaned = sub("-.*", "", missing),
                     stringsAsFactors = FALSE)
pConvert = gconvert(missing$cleaned, target = "HGNC", filter_na = FALSE, 
                    mthreshold = 1)[c("alias", "name")]
missing = merge(missing, pConvert, by.x = "cleaned", by.y = "alias")
missing = missing[missing$name != "N/A", ]
missing$name = as.character(missing$name)
for (i in 1:nrow(missing)) {
  pSitePlus[pSitePlus$SUB_ACC_ID == missing[i, "genes"], "SUB_GENE"] = missing[i, "name"]
}
pSitePlus = pSitePlus[pSitePlus$SUB_GENE != "", ]  # Remove rows without substrate gene names
rm(missing, pConvert)

# Define "Site" column as GENE_Psite
pSitePlus$Site = paste0(pSitePlus$SUB_GENE, "_", pSitePlus$SUB_MOD_RSD)   # Contains 3 p-His

# Plot historgram of # of substrates per kinase
cairo_ps("Substrates_per_kinase.eps")
hist(table(pSitePlus$GENE), breaks = 20, main = "Distribution of Number of Substrates per Kinase", 
     xlab = "Number of Substrates")
dev.off()

# Plot top 12 kinases with most substrates
cairo_ps("top10_Substrates_per_kinase_barplot.eps")
bar.data = table(pSitePlus$GENE)[table(pSitePlus$GENE) > 200]
bar.data = bar.data[order(bar.data)]
barplot(bar.data, main = "Top 12 Annotated Kinases", xlab = "Kinase", 
        ylab = "Number of Substrate", las = 2)
dev.off()


# Construct adjacency matrix from pSitePlus data
pMatrix = matrix(0, nrow = length(unique(pSitePlus$Site)), 
                 ncol = length(unique(pSitePlus$GENE)),
                 dimnames = list(row = unique(pSitePlus$Site),
                                 col = unique(pSitePlus$GENE)))
for (i in 1:nrow(pSitePlus)) {
  pMatrix[pSitePlus[i, "Site"], pSitePlus[i, "GENE"]] = 1   #Assign 1 if interaction exists
}

# Kinase activity matrix (m samples by n kinases)
KSEA = function(phospho, pMatrix) {
  # phospho = a list of two elements (phospho data and site information)
  # pMatrix = an adjacency matrix detailing relationship between kinase and substrate
  # Method described in Phosphoproteomics-based Profiling of Kinase Activities in Cancer Cells (BioRxiv 2016)
  store = matrix(NA, nrow = ncol(pMatrix), ncol = ncol(phospho$values),
                 dimnames = list(colnames(pMatrix), colnames(phospho$values)))   #Initialize matrix
  zScore = store   #Initialize matrix
  pValue = store   #Initialize matrix
  
  mP = mean(as.matrix(phospho$values))
  sdP = sd(as.matrix(phospho$values))
  
  for (j in 1:ncol(pMatrix)) {
    # Iterate through each kinase
    sites = rownames(pMatrix)[pMatrix[, j] == 1]
    phospho_id = phospho$values[phospho$index$Site %in% sites, ]
    for (i in 1:ncol(phospho$values)) {
      # Iterate through each sample
      if (nrow(phospho_id) != 0) {
        store[j, i] = mean(phospho_id[, i])   #Calculate KSEA score
        zScore[j, i] = (store[j, i] - mP * sqrt(nrow(phospho_id))) / sdP
      }
    }
  }
  pValue = 2*pnorm(-abs(zScore))   # 2-sided test
  
  return(list(score = store, zScore = zScore, pValue = pValue))
}
b.KSEA = KSEA(b.phospho, pMatrix)
b.KSEA = lapply(b.KSEA, function(x) x[complete.cases(x), ])   #Remove rows with NA

# Plot Heatmap of Kinase Activity Score for BRCA data
require(gplots)
require(RColorBrewer)
cairo_ps("KSEA_BRCA_heatmap.eps")
col.pan <- colorpanel(100, "blue", "white", "red")
##Color columns by subtype
colColors = colnames(b.KSEA$score)
colColors = sub("-.$", "", colColors)
colColors = as.character(b.clinical$PAM50.mRNA[match(colColors, b.clinical$Complete.TCGA.ID)])
colColors[1:3] = "Normal"
temp.data = b.KSEA$score
colnames(temp.data) = colColors
colColors = brewer.pal(5, "Dark2")[factor(colColors)]
heatmap.2(temp.data, col = col.pan, Rowv = TRUE, Colv = TRUE, scale = "none",
          trace = "none", dendrogram = "column", cexRow=0.5, cexCol=0.5,
          margin=c(4,5), lhei = c(1.25,8), lwid = c(1.5,6), ColSideColors = colColors,
          density.info = "histogram", key.par = list(cex=0.5), 
          key.title = "Distribution of Z-scores", key.xlab = "", key.ylab = "Counts")
rm(colColors, temp.data)
dev.off()

# Impletement KSEA with NetworKIN
require(stringr)
b.networKIN = rownames(b.phospho$values)
b.networKIN = data.frame(id = sub("\\..*", "", b.networKIN),
                         residue = unlist(lapply(b.networKIN,
                                                 function(x) str_extract(x, "(?<=:[styh])[0-9]*"))),
                         position = toupper(unlist(lapply(b.networKIN, 
                                                          function(x) str_extract(x, "(?<=:)[styh]")))))
# Write file for input into networKIN (HUMAN - NCBI62 database)
rm_genes = c("NP_001074295", "NP_001092906", "NP_001098549",
             "NP_004049	", "NP_056271", "NP_060686", "NP_112181",
             "NP_443134", "NP_542937", "NP_612206", "NP_699165",
             "NP_851851", "NP_004049")   #Remove proteins that trip up networKIN
write.table(b.networKIN[!(b.networKIN$id %in% rm_genes), ], 
            file = "BRCA_networKIN_input_mod.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# Read in networKIN input

#####





  
  


