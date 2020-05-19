############ :::: ISOFORM DRUG TARGETS PROJECT :::: #############
## ::: Load packages :::
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")

if(!require("data.table")){
  install.packages("data.table")
  library(data.table)
}

if(!require("ggplot2")){
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require("reshape2")){
  install.packages("reshape2")
  library(reshape2)
}

if(!require("RColorBrewer")){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

if(!require("XML")){
  install.packages("XML")
  library(XML)
}

if(!require("Biostrings")){
  biocLite("Biostrings")
  library(Biostrings)
}

if(!require("mygene")){
  biocLite("mygene")
  library(mygene)
}

if(!require("Biobase")){
  biocLite("Biobase")
  library(Biobase)
}

if(!require("PharmacoGx")){
  install.packages("/Users/yjq6148/Desktop/quest/Jerry/Isoform_drug_targets/PharmacoGx", repos = NULL, type="source")
  library(PharmacoGx)
}

if(!require("xtable")){
  install.packages("xtable")
  library(xtable)
  options(xtable.floating = FALSE)
  options(xtable.timestamp = "")
}

if(!require("msa")){
  biocLite("msa")
  library(msa)
}

if(!require("seqinr")){
  biocLite("seqinr")
  library(seqinr)
}

if(!require("doParallel")){
  install.packages("doParallel")
  library(doParallel)
}

if(!require("dplyr")){
  install.packages("dplyr")
  library(dplyr)
}

if(!require("tidyr")){
  install.packages("tidyr")
  library(tidyr)
}

if(!require("tidyverse")){
  install.packages("tidyverse")
  library(tidyverse)
}

registerDoParallel(cores = 5)


## ::::: Load helper function :::::
source("/Users/yjq6148/Desktop/quest/Jerry/Isoform_drug_targets/helpers.R")


## ::::: Load dgidb data and summarize number of transcript variants/protein isoforms :::::
# Directories
DATA <- "/Users/yjq6148/Desktop/quest/Jerry/Isoform_drug_targets/data/"
Anno.dir <- paste0(DATA,"../annotation/")
FDA.approved <- 0
Output <- ifelse(FDA.approved == 0,paste0(DATA,"../results/all/"),paste0(DATA,"../results/FDAapproved/"))
dgidb.dir <- paste0(DATA,"dgidb/")
drug.interactions <- read.table(paste0(dgidb.dir,"interactions.tsv"),sep = '\t',fill = TRUE,quote = "",header = TRUE)
# gene.categories <- read.table(paste0(dgidb.dir,"categories.tsv"),sep = '\t',fill = TRUE,quote = "",header = TRUE)

print(length(unique(drug.interactions$gene_name))) # 2994 unique genes
print(length(unique(drug.interactions$drug_name))) # 6538 unique drugs

# rm(gene.categories)

## Write all druggable genes to a file as input to dgidb webpage (to filter antineoplastic), in total 2994 druggable genes in dgidb
# write.table(unique(drug.interactions$gene_name),file=paste0(dgidb.dir,"all_druggable_genes.txt"),row.names = FALSE,col.names = FALSE, quote = FALSE)
rm(drug.interactions)

## ::: Load resulting dgidb data :::
## 6688 cancer drug interactions (1447 genes ~ 883 drugs) and 3477 FDA approved cancer drug interactions (1122 genes ~ 280 drugs)
cancer.drug.interactions <- read.table(file = paste0(dgidb.dir,"interactions_antineoplastic.tsv"),sep = '\t',fill = TRUE,quote = "",header = TRUE)[,4:5]
print(length(unique(cancer.drug.interactions$gene))) # 1447 unique genes
print(length(unique(cancer.drug.interactions$drug))) # 883 unique drugs
cancer.drug.interactions.FDA <- read.table(file = paste0(dgidb.dir,"interactions_antineoplastic_FDA.tsv"),sep = '\t',fill = TRUE,quote = "",header = TRUE)[,4:5]
print(length(unique(cancer.drug.interactions.FDA$gene))) # 1122 unique genes
print(length(unique(cancer.drug.interactions.FDA$drug))) # 280 unique drugs

cancer.drug.interactions$`FDA_approved` <- ifelse(cancer.drug.interactions$drug %in% cancer.drug.interactions.FDA$drug,'Yes','No')
cancer.drug.interactions$drug <- toupper(cancer.drug.interactions$drug)
rm(cancer.drug.interactions.FDA)

if(FDA.approved  == 1){
  cancer.drug.interactions <<- cancer.drug.interactions[which(cancer.drug.interactions$FDA_approved == "Yes"),]
}

## Summarize drugs for each gene
cancer.drug.interactions <- cancer.drug.interactions %>%
  group_by(gene) %>%
  summarise(drug_list = paste(drug,collapse=", ")) %>%
  ungroup()

## Annotation (GRCh38)
# Mygene --> gene level annotation
mygene.search <- as.data.frame(queryMany(cancer.drug.interactions$gene,scopes="symbol",fields=c("ensembl.gene","entrezgene"),species="human"))
mygene.search$ensembl <- sapply(mygene.search$ensembl,unlist)
mygene.search$ensembl <- sapply(mygene.search$ensembl, function(x) {ifelse(is.null(x), NA, paste(unlist(x),collapse = ','))})
mygene.search <- mygene.search %>%
  group_by(query) %>%
  fill(ensembl,entrezgene) %>%
  fill(ensembl,entrezgene,.direction="up")
mygene.search <- mygene.search %>% 
  mutate(ensembl= strsplit(as.character(ensembl), ",")) %>% 
  unnest(ensembl)
mygene.search <- subset(mygene.search,select=-c(X_id,X_score,notfound))
mygene.search <- mygene.search[ ,order(names(mygene.search))]
colnames(mygene.search) <- c("ensembl_id","entrez_id","gene")
cancer.drug.interactions <- merge(cancer.drug.interactions,mygene.search,by="gene")
cancer.drug.interactions <- na.omit(cancer.drug.interactions) # 1258 unique genes
print(length(unique(cancer.drug.interactions$gene_symbol)))
print(length(unique(cancer.drug.interactions$targeting_drugs)))
colnames(cancer.drug.interactions) <- c("gene_symbol","targeting_drugs","ensembl_id","entrez_id")
rm(mygene.search)

# Transcript
ensembl.transcript <- read.table(paste0(Anno.dir,"ensembl/Ensembl.GRCh37.p13.genes.transcripts.txt"),sep='\t',header=TRUE,quote="",fill=TRUE)
ensembl.transcript.summary <- ensembl.transcript %>%
  group_by(Gene.stable.ID) %>%
  summarise(transcript_variants = n_distinct(Transcript.stable.ID)) %>%
  ungroup()
colnames(ensembl.transcript.summary)[1] <- "ensembl_id"
cancer.drug.interactions <- merge(cancer.drug.interactions,ensembl.transcript.summary,by="ensembl_id") # 1205 unique genes
rm(ensembl.transcript.summary)

# Protein
ensembl.protein <- read.table(paste0(Anno.dir,"ensembl/Ensembl.GRCh37.p13.genes.transcripts.proteins.txt"),sep='\t',header=TRUE,quote="",fill=TRUE)

# Protein sequence
ensembl.protein.seq <- readAAStringSet(filepath = paste0(Anno.dir,"ensembl/Homo_sapiens.GRCh37.pep.all.fa"), format = "fasta")
ensembl.protein.seq <- as.data.frame(ensembl.protein.seq)
colnames(ensembl.protein.seq) <- "peptide_sequence"
ensembl.protein.seq$`Protein.stable.ID` <- gsub("^(ENSP[0-9]*)\\..*","\\1",row.names(ensembl.protein.seq))
ensembl.protein.seq <- ensembl.protein.seq[!duplicated(ensembl.protein.seq),]
ensembl.protein <- merge(ensembl.protein,ensembl.protein.seq,by="Protein.stable.ID")
ensembl.protein.summary <- ensembl.protein %>%
  group_by(Gene.stable.ID) %>%
  summarise(protein_isoforms = n_distinct(peptide_sequence)) %>% # Some IDs have duplicated sequences, so count n_distinct based on sequence
  ungroup()
colnames(ensembl.protein.summary)[1] <- "ensembl_id"
cancer.drug.interactions <- merge(cancer.drug.interactions,ensembl.protein.summary,by="ensembl_id") # 1180 unique genes
rm(ensembl.protein.summary)
rm(ensembl.protein.seq)

transcript.count <- cancer.drug.interactions %>%
  group_by(transcript_variants) %>%
  summarise(freq=n()) 
transcript.count[nrow(transcript.count)+1,] <- c(">5",colSums(transcript.count[5:nrow(transcript.count),"freq"]))
print(sum(cancer.drug.interactions$transcript_variants)/(length(unique(cancer.drug.interactions$gene_symbol)))) # 9.23

protein.count <- cancer.drug.interactions %>%
  group_by(protein_isoforms) %>%
  summarise(freq=n()) 
protein.count[nrow(protein.count)+1,] <- c(">5",colSums(protein.count[5:nrow(protein.count),"freq"]))
print(sum(cancer.drug.interactions$protein_isoforms)/(length(unique(cancer.drug.interactions$gene_symbol)))) # 5.22


# Write output to files
if(FDA.approved  == 1){
  write.csv(cancer.drug.interactions,file=paste0(Output,"cancer.drug.interactions.FDAapproved.csv"))
}else{
  write.csv(cancer.drug.interactions,file=paste0(Output,"cancer.drug.interactions.all.csv"))
}


######### ::::: BIOLIP BINDING POCKETS ::::: ###########
if(file.exists(paste0(Output,"BioLiP.human.pdb.drugbank.chembl.seq.annotated.csv"))){
  BioLiP.human.CHEMBL.seq <- read.csv(file=paste0(Output,"BioLiP.human.pdb.drugbank.chembl.seq.annotated.csv"),header=TRUE,row.names = 1)
}else{
  # Read in all BioLiP entries
  BioLiP <- read.table(paste0(DATA,"biolip/all_summaries_nr/BioLiP_all_nr.txt"),sep="\t",header = FALSE,fill=TRUE)
  BioLiP <- BioLiP[!duplicated(BioLiP),] # 180750 non-redundant entries (Binding site residues identity <= 90%; Receptor sequence identity <= 90%)
  colnames(BioLiP) <- c("PDB_ID",
                        "PDB_chain",
                        "Resolution",
                        "Binding_site_number_code",
                        "Ligand_ID",
                        "Ligand_chain",
                        "Ligand_serial_number",
                        "Binding_site_residues",
                        "Renumbered_binding_site_residues",
                        "Catalytic_site_residues",
                        "Renumbered_catalytic_site_residues",
                        "EC_number",
                        "GO_terms",
                        "Binding_affinity_manual",
                        "Binding_affinity_Binding_MOAD",
                        "Binding_affinity_PDBbind_CN",
                        "Binding_affinity_BindingDB",
                        "UniProt_ID",
                        "PubMed_ID",
                        "Receptor_sequence")
  BioLiP$PDB_ID <- toupper(BioLiP$PDB_ID)
  
  # # # Write to Clipboard for searching UniProt
  # clip <- pipe("pbcopy","w")
  # write.table(paste(unique(BioLiP$UniProt_ID),collapse = " "),file=clip)
  # close(clip)
  
  #### Annotation
  
  # if(anno == "gencode"){
  #   pdb.gencode <- read.table(paste0(Anno.dir,"gencode/gencode.v28lift37.metadata.PDB"),header=FALSE,quote="",fill=TRUE)
  #   colnames(pdb.gencode) <- c("Transcript.stable.ID","PDB_ID")
  #   pdb.gencode$Transcript.stable.ID <- gsub("^(ENST[0-9]*)\\.[0-9]","\\1",pdb.gencode$Transcript.stable.ID)
  #   pdb.gencode <- merge(pdb.gencode,ensembl.transcript[,c("Gene.stable.ID","Transcript.stable.ID")],by="Transcript.stable.ID")
  #   pdb.gencode$Transcript.stable.ID <- NULL
  #   pdb.gencode <- pdb.gencode[!duplicated(pdb.gencode),]
  #   BioLiP.human <<- unique(merge(BioLiP,pdb.gencode,by="PDB_ID")[,c("PDB_ID","PDB_chain","Ligand_ID","Ligand_chain","Binding_site_residues","Renumbered_binding_site_residues","Gene.stable.ID")])
  # }
  # uniprot <- read.table(paste0(Anno.dir,"uniprot/UniProt_to_gene_symbol.txt"),header=TRUE,quote="",fill=TRUE)
  # uniprot$gene_symbol <- toupper(uniprot$gene_symbol)
  # BioLiP1 <- unique(merge(BioLiP,uniprot,by="UniProt_ID")[,c("UniProt_ID","PDB_ID","PDB_chain","Ligand_ID","Ligand_chain","Binding_site_residues","Renumbered_binding_site_residues","gene_symbol")])
  mygene.search <- queryMany(BioLiP$PDB_ID,
                             scopes="pdb",
                             fields=c("ensembl.gene","entrezgene"),
                             species="human",
                             return.as = "DataFrame")
  mygene.search <- mygene.search[which(is.na(mygene.search$notfound)),]
  mygene.search$ensembl <- sapply(mygene.search$ensembl,unlist)
  mygene.search$ensembl <- sapply(mygene.search$ensembl, function(x) {ifelse(is.null(x), NA, paste(unlist(x),collapse = ','))})
  mygene.search <- as.data.frame(mygene.search)
  mygene.search <- mygene.search %>%
    group_by(query) %>%
    fill(ensembl,entrezgene) %>%
    fill(ensembl,entrezgene,.direction="up")
  mygene.search <- mygene.search %>% 
    mutate(ensembl= strsplit(as.character(ensembl), ",")) %>% 
    unnest(ensembl)
  mygene.search <- subset(mygene.search,select=-c(X_id,X_score,notfound,ensembl.gene))
  mygene.search <- mygene.search[ ,order(names(mygene.search))]
  colnames(mygene.search) <- c("ensembl_id","entrez_id","PDB_ID")
  BioLiP.human <- unique(merge(BioLiP,mygene.search,by="PDB_ID")[,c("PDB_ID","PDB_chain","Ligand_ID","Ligand_chain","Binding_site_residues","Renumbered_binding_site_residues","ensembl_id")])
  # 66188 bindings corresponding to 4470 human genes and 8995 ligands
  
  ## Not all ligands are cancer drugs, so filter those ligands belong to cancer drugs
  # Obtain ligand ID annotation from PDB
  # Write to Clipboard for searching PDB
  # clip <- pipe("pbcopy","w")
  # write.table(paste(unique(BioLiP.human$Ligand_ID),collapse = " "),file=clip)
  # close(clip)
  # pdb.drugs <- unique(read.table(paste0(Anno.dir,"pdb/pdb_to_drugbank.csv"),sep=',',header=TRUE,fill=TRUE)[,c(1:3,5)])
  # pdb.drugs <- pdb.drugs %>% 
  #   mutate(Ligand.ID= strsplit(as.character(Ligand.ID), ",")) %>% 
  #   unnest(Ligand.ID)
  
  # Obtain ligand ID annotation from DRUGBANK
  drugbank.drugs <- unique(read.csv(paste0(Anno.dir,"drugbank/drug links.csv"),header=TRUE,fill=TRUE)[,c(1,2,11)])
  drugbank.drugs <- drugbank.drugs[which(drugbank.drugs$HET.ID != ""),]
  colnames(drugbank.drugs) <- c("Drugbank_id","Generic_name","Ligand_ID")
  drugbank.drugs$Generic_name <- toupper(drugbank.drugs$Generic_name)
  
  # Merge with binding data
  BioLiP.human.drugbank <- merge(BioLiP.human,drugbank.drugs,by="Ligand_ID") 
  BioLiP.human.drugbank <- BioLiP.human.drugbank[which(!is.na(BioLiP.human.drugbank$ensembl_id)),] # 15010 bindings corresponding to 2280 human genes and 2298 drugs
  write.csv(BioLiP.human.drugbank,file=paste0(Output,"BioLiP.human.pdb.drugbank.annotated.csv"))
  
  
  
  
  ######### ::::: Merge dgidb (cancer gene/isoform-drug pair, in ChEMBL) and BioLiP (binding data, in DrugBank) ::::: #########
  ## Reading in full DrugBank database from XML
  drugbank <- xmlParse(paste0(Anno.dir,"drugbank/DrugBank.xml"))
  drugbank.list <- xmlToList(drugbank)
  drugbank.list[[11034]] <- NULL # Remove last element
  
  ## Extracting id, name and CHEMBL id
  drugbank.id <- sapply(drugbank.list,function(x){x[["drugbank-id"]][["text"]]})
  drugbank.name <- sapply(drugbank.list,function(x) toupper(x[["name"]]))
  drugbank.identifiers <- lapply(drugbank.list,
                            function(x){
                              x[["external-identifiers"]]
                              })
  drugbank.CHEMBL <- c()
  hasCHEMBL <- 0
  for(drug in drugbank.identifiers){
    for(identifier in drug){
      if(identifier[["resource"]] == "ChEMBL"){
        drugbank.CHEMBL <- c(drugbank.CHEMBL,identifier[["identifier"]])
        hasCHEMBL <- 1
      }
    }
    if(hasCHEMBL != 1){
      drugbank.CHEMBL <- c(drugbank.CHEMBL,NA)
    }
    hasCHEMBL <- 0
  }
  
  drugbank.to.CHEMBL <- cbind.data.frame(drugbank.id,drugbank.name,drugbank.CHEMBL)
  colnames(drugbank.to.CHEMBL) <- c("Drugbank_id","Drugbank_name","ChEMBL_id")
  
  # Merge with binding data
  BioLiP.human.CHEMBL <- merge(BioLiP.human.drugbank,drugbank.to.CHEMBL,by="Drugbank_id")
  BioLiP.human.CHEMBL <- BioLiP.human.CHEMBL[which(!is.na(BioLiP.human.CHEMBL$ChEMBL_id)),] # 10714 bindings (1790 genes ~ 1633 drugs)
  BioLiP.human.CHEMBL$Drugbank_name <- NULL
  write.csv(BioLiP.human.CHEMBL,file=paste0(Output,"BioLiP.human.pdb.drugbank.chembl.annotated.csv")) # Now binding data is ChEMBL annotated
  
  ## Reading in pdb sequence information
  pdb.seq <- readAAStringSet(filepath = paste0(Anno.dir,"pdb/BioLiP.pdb.sequences.fasta"), format = "fasta")
  pdb.seq <- as.data.frame(pdb.seq)
  colnames(pdb.seq) <- "peptide_sequence"
  pdb.seq$`PDB_ID` <- gsub("^(.*):.*","\\1",row.names(pdb.seq))
  pdb.seq$`PDB_chain` <- gsub("^.*:([A-Z])\\|.*","\\1",row.names(pdb.seq))
  pdb.seq <- pdb.seq[!duplicated(pdb.seq),]
  
  BioLiP.human.CHEMBL.seq <- merge(BioLiP.human.CHEMBL,pdb.seq,by=c("PDB_ID","PDB_chain"))
  write.csv(BioLiP.human.CHEMBL.seq,file=paste0(Output,"BioLiP.human.pdb.drugbank.chembl.seq.annotated.csv")) # with seq
  
  ## Make binding pocket sequence
  binding.pocket <- makeBindingPocket(BioLiP.human.CHEMBL.seq)
  BioLiP.human.CHEMBL.seq$Binding_pocket <- binding.pocket
  BioLiP.human.CHEMBL.seq$Binding_pocket <- as.character(BioLiP.human.CHEMBL.seq$Binding_pocket)
  write.csv(BioLiP.human.CHEMBL.seq,file=paste0(Output,"BioLiP.human.pdb.drugbank.chembl.seq.annotated.csv")) # with pocket
}


## Compute complete cancer drug interactions table
drugs <- read.table(paste0(dgidb.dir,"drugs.tsv"),sep = '\t',fill = TRUE,quote = "",header = TRUE)
drugs <- unique(drugs[,c(1:3)])
drugs$drug_name <- toupper(drugs$drug_name)
drugs[which(drugs$drug_name == ""),"drug_name"] <- drugs[which(drugs$drug_name == ""),"drug_claim_name"]
drugs$drug_claim_name <- NULL
drugs$drug_name <- toupper(drugs$drug_name)
drugs <- unique(drugs)

if(FDA.approved  == 1){
  cancer.drug.interactions <- read.table(file = paste0(dgidb.dir,"interactions_antineoplastic_FDA.tsv"),sep = '\t',fill = TRUE,quote = "",header = TRUE)[,4:5]
}else{
  cancer.drug.interactions <- read.table(file = paste0(dgidb.dir,"interactions_antineoplastic.tsv"),sep = '\t',fill = TRUE,quote = "",header = TRUE)[,4:5]
}
cancer.drug.interactions$drug <- toupper(cancer.drug.interactions$drug)
colnames(cancer.drug.interactions) <- c("gene_name","drug_name")

cancer.drug.interactions.complete <- merge(cancer.drug.interactions,drugs,by="drug_name") # 3476/6685 drug interactions

## Annotation (GRCh37.p13)
# Mygene --> gene level annotation
mygene.search <- as.data.frame(queryMany(unique(cancer.drug.interactions.complete$gene_name),scopes="symbol",fields=c("ensembl.gene","entrezgene"),species="human"))
mygene.search$ensembl <- sapply(mygene.search$ensembl,unlist)
mygene.search$ensembl <- sapply(mygene.search$ensembl, function(x) {ifelse(is.null(x), NA, paste(unlist(x),collapse = ','))})
mygene.search <- mygene.search %>%
  group_by(query) %>%
  fill(ensembl,entrezgene) %>%
  fill(ensembl,entrezgene,.direction="up")
mygene.search <- mygene.search %>% 
  mutate(ensembl= strsplit(as.character(ensembl), ",")) %>% 
  unnest(ensembl)
mygene.search <- subset(mygene.search,select=-c(X_id,X_score,notfound))
mygene.search <- unique(mygene.search[ ,order(names(mygene.search))])
colnames(mygene.search) <- c("ensembl_id","entrez_id","gene_name")
cancer.drug.interactions.complete <- merge(cancer.drug.interactions.complete,mygene.search,by="gene_name")
cancer.drug.interactions.complete <- unique(na.omit(cancer.drug.interactions.complete)) # 1110/1432 unique genes (symbol)
cancer.drug.interactions.complete <- cancer.drug.interactions.complete[,1:4]
rm(mygene.search)


## Merge BioLiP with dgidb
binding.data <- unique(BioLiP.human.CHEMBL.seq[,c(8,10,12)])
colnames(binding.data) <- c("ensembl_id","chembl_id","binding_pocket")
cancer.binding.data <- merge(cancer.drug.interactions.complete,binding.data,by=c("ensembl_id","chembl_id"))
length(unique(cancer.binding.data$gene_name)) # 67/112 genes
length(unique(cancer.binding.data$chembl_id)) # 51/137 drugs

## Correct for Ensembl ID of some genes, which have multiple IDs that are not primary assembly
cancer.binding.data[which(cancer.binding.data$gene_name == "DDR1"),"ensembl_id"] <- "ENSG00000204580"
cancer.binding.data[which(cancer.binding.data$gene_name == "PSMB1"),"ensembl_id"] <- "ENSG00000008018"
cancer.binding.data[which(cancer.binding.data$gene_name == "PSMB3"),"ensembl_id"] <- "ENSG00000108294"
cancer.binding.data[which(cancer.binding.data$gene_name == "PSMB8"),"ensembl_id"] <- "ENSG00000204264"
cancer.binding.data <- unique(cancer.binding.data)

## Write output
write.csv(cancer.binding.data,file=paste0(Output,"cancer.binding.data.gene.level.csv")) # Gene level cancer binding data



######### ::::: Multiple alignments ::::: #########
## Merge data frame with protein isoform level sequence information
# Protein
ensembl.protein <- read.table(paste0(Anno.dir,"ensembl/Ensembl.GRCh37.p13.genes.transcripts.proteins.txt"),sep='\t',header=TRUE,quote="",fill=TRUE)

# Protein sequence
ensembl.protein.seq <- readAAStringSet(filepath = paste0(Anno.dir,"ensembl/Homo_sapiens.GRCh37.pep.all.fa"), format = "fasta")
ensembl.protein.seq <- as.data.frame(ensembl.protein.seq)
colnames(ensembl.protein.seq) <- "peptide_sequence"
ensembl.protein.seq$`Protein.stable.ID` <- gsub("^(ENSP[0-9]*)\\..*","\\1",row.names(ensembl.protein.seq))
ensembl.protein.seq <- ensembl.protein.seq[!duplicated(ensembl.protein.seq),]
ensembl.protein <- merge(ensembl.protein,ensembl.protein.seq,by="Protein.stable.ID")
ensembl.protein <- unique(ensembl.protein[,c(1:3,10)])
colnames(ensembl.protein) <- c("protein_id","ensembl_id","gene_name","protein_sequence")
cancer.binding.data.protein <- merge(cancer.binding.data,ensembl.protein,by=c("gene_name","ensembl_id"))
length(unique(cancer.binding.data.protein$gene_name)) # 67 genes
length(unique(cancer.binding.data.protein$chembl_id)) # 51 drugs
write.csv(cancer.binding.data.protein,file=paste0(Output,"cancer.binding.data.protein.level.csv")) # Protein level cancer binding data

## Add binding residues as a separate row, in order to make fasta file in following step
cancer.binding.data.final <- addBindingSeqRow(cancer.binding.data.protein)

## Make fasta file
makeFastaFiles(cancer.binding.data.final,paste0(Output,"fasta/"))


## Output multiple alignments
setwd(paste0(Output,"alignments/"))
files <- list.files(path=paste0(Output,"fasta/"), pattern="*.fasta", full.names=T, recursive=FALSE)
lapply(files, function(file) {
  print(file)
  myseqs <- readAAStringSet(file)
  myalignment <- msa(myseqs,"ClustalOmega",order="input")
  msaPrettyPrint(myalignment, 
                 output="pdf", 
                 file = paste0(Output,"alignments/",gsub("^.*/fasta/(.*)\\..*","\\1",file),".pdf"), 
                 showNames="left",
                 showLogo="top", 
                 logoColors = "structure",
                 shadingMode = "similar",
                 askForOverwrite=FALSE, 
                 verbose=FALSE,
                 furtherCode =c("\\constosingleseq{1}",
                                "\\threshold[50]{0}"))})


######### ::::: Multiple alignments to GRCh38 ::::: #########
## Merge data frame with protein isoform level sequence information
# Protein
ensembl.protein <- read.table(paste0(Anno.dir,"ensembl/Ensembl.GRCh38.p7.genes.transcripts.proteins.txt"),sep='\t',header=TRUE,quote="",fill=TRUE)

# Protein sequence
ensembl.protein.seq <- readAAStringSet(filepath = paste0(Anno.dir,"ensembl/Homo_sapiens.GRCh38.pep.all.fa"), format = "fasta")
ensembl.protein.seq <- as.data.frame(ensembl.protein.seq)
colnames(ensembl.protein.seq) <- "peptide_sequence"
ensembl.protein.seq$`Protein.ID` <- gsub("^(ENSP[0-9]*)\\..*","\\1",row.names(ensembl.protein.seq))
ensembl.protein.seq <- ensembl.protein.seq[!duplicated(ensembl.protein.seq),]
ensembl.protein <- merge(ensembl.protein,ensembl.protein.seq,by="Protein.ID")
ensembl.protein <- unique(ensembl.protein[,c(1:3,5:6)])
colnames(ensembl.protein) <- c("protein_id","ensembl_id","gene_name","protein_name","protein_sequence")
cancer.binding.data.protein <- merge(cancer.binding.data,ensembl.protein,by=c("gene_name","ensembl_id"))
length(unique(cancer.binding.data.protein$gene_name)) # 67 genes
length(unique(cancer.binding.data.protein$chembl_id)) # 51 drugs
write.csv(cancer.binding.data.protein,file=paste0(Output,"/cancer.binding.data.protein.level.GRCh38.csv")) # Protein level cancer binding data

## Add binding residues as a separate row, in order to make fasta file in following step
cancer.binding.data.final <- addBindingSeqRow(cancer.binding.data.protein)

## Make fasta file
makeFastaFiles(cancer.binding.data.final,paste0(Output,"all/fasta.grch38/"))


## Output multiple alignments
setwd(paste0(Output,"/alignments.grch38/"))
files <- list.files(path=paste0(Output,"/fasta.grch38/"), pattern="*.fasta", full.names=T, recursive=FALSE)
lapply(files, function(file) {
  print(file)
  myseqs <- readAAStringSet(file)
  myalignment <- msa(myseqs,"ClustalOmega",order="input")
  msaPrettyPrint(myalignment, 
                 output="pdf", 
                 file = paste0(Output,"/alignments.grch38/",gsub("^.*/fasta.grch38/(.*)\\..*","\\1",gsub(" ","_",file)),".pdf"), 
                 showNames="left",
                 showLogo="top", 
                 logoColors = "structure",
                 shadingMode = "similar",
                 askForOverwrite=FALSE, 
                 verbose=FALSE,
                 furtherCode =c("\\constosingleseq{1}",
                                "\\threshold[50]{0}"))})




######## ::::: DRUG SENSITIVITY AND PERTURBATION MODELING USING PHARMACOGX PACKAGE ::::: ########
# availablePSets()
# GDSC <- downloadPSet("GDSC")
# CCLE <- downloadPSet("CCLE")
load(paste0(DATA,"../PSets/CCLE_hs.RData"), verbose=TRUE) ## Ensembl 87
#downloadPSet("GDSC", saveDir="data")
# load(paste0(DATA,"../PSets/GDSC.RData"), verbose=TRUE)

#### EXAMPLE: ABL1 (3 isoforms)
# feature1 <- fNames(CCLE, "isoforms")[which(featureInfo(CCLE,"isoforms")$TranscriptName == "ABL1-001")]
# drug.sensitivity1 <- drugSensitivitySig(pSet = CCLE,
#                                         mDataType="isoforms",
#                                         drugs=c("Nilotinib"),
#                                         features=feature1,
#                                         sensitivity.measure="auc_published",
#                                         molecular.summary.stat="median",
#                                         sensitivity.summary.stat="median",
#                                         verbose=TRUE)
# 
# feature2 <- fNames(CCLE, "isoforms")[which(featureInfo(CCLE,"isoforms")$TranscriptName == "ABL1-002")]
# drug.sensitivity2 <- drugSensitivitySig(pSet = CCLE,
#                                         mDataType="isoforms",
#                                         drugs=c("Nilotinib"),
#                                         features=feature2,
#                                         sensitivity.measure="auc_published",
#                                         molecular.summary.stat="median",
#                                         sensitivity.summary.stat="median",
#                                         verbose=TRUE)
# 
# feature3 <- fNames(CCLE, "isoforms")[which(featureInfo(CCLE,"isoforms")$TranscriptName == "ABL1-003")]
# drug.sensitivity3 <- drugSensitivitySig(pSet = CCLE,
#                                         mDataType="isoforms",
#                                         drugs=c("Nilotinib"),
#                                         features=feature3,
#                                         sensitivity.measure="auc_published",
#                                         molecular.summary.stat="median",
#                                         sensitivity.summary.stat="median",
#                                         verbose=TRUE)

feature <- fNames(CCLE, "isoforms")[which(grepl("^ABL1",featureInfo(CCLE,"isoforms")$TranscriptName))]
drug.sensitivity <- drugSensitivitySig(pSet = CCLE,
                                       mDataType="isoforms",
                                       drugs=c("Nilotinib"),
                                       features=feature,
                                       tissues="haematopoietic_and_lymphoid_tissue",
                                       sensitivity.measure="auc_published",
                                       molecular.summary.stat="median",
                                       sensitivity.summary.stat="median",
                                       verbose=TRUE)

# sensitivity <- rbind(drug.sensitivity1,drug.sensitivity2,drug.sensitivity3)
sensitivity <- as.data.frame(drug.sensitivity)
rownames(sensitivity) <- c("Nilotinib + ABL1-001",
                           "Nilotinib + ABL1-002",
                           "Nilotinib + ABL1-003")
colnames(sensitivity) <- dimnames(drug.sensitivity)[[3]]
xtable(sensitivity, caption='P Value of Gene-Drug Association')


#### Do this for all examples found
drugs <- unique(cancer.binding.data.final[,c("gene_name","drug_name")])
drugs$gene_name <- as.character(drugs$gene_name)
drugs$drug_name <- as.character(drugs$drug_name)
drugs <- drugs[which(drugs$drug_name %in% c("VANDETANIB", # Approved
                                            "DOVITINIB", # Phase III
                                            "ERLOTINIB", # Approved
                                            "LAPATINIB", # Approved
                                            "NILOTINIB", # Approved
                                            "PD-0325901", # Discontinued
                                            "PALBOCICLIB", # Phase II
                                            "CRIZOTINIB", # Approved
                                            "PLX-4720", # Pre-clinical
                                            "CHIR-265", # Phase I
                                            "SORAFENIB")),] # Approved

# 11 out of 15 unique drugs that have drug response information
# Change names
drugs[which(drugs$drug_name == "DOVITINIB"),"drug_name"] <- "TKI258"
drugs[which(drugs$drug_name == "PALBOCICLIB"),"drug_name"] <- "PD-0332991"
drugs[which(drugs$drug_name == "PLX-4720"),"drug_name"] <- "PLX4720"
drugs[which(drugs$drug_name == "CHIR-265"),"drug_name"] <- "RAF265"
drugs[which(drugs$drug_name == "VANDETANIB"),"drug_name"] <- "Vandetanib"
drugs[which(drugs$drug_name == "ERLOTINIB"),"drug_name"] <- "Erlotinib"
drugs[which(drugs$drug_name == "LAPATINIB"),"drug_name"] <- "lapatinib"
drugs[which(drugs$drug_name == "NILOTINIB"),"drug_name"] <- "Nilotinib"
drugs[which(drugs$drug_name == "CRIZOTINIB"),"drug_name"] <- "Crizotinib"
drugs[which(drugs$drug_name == "SORAFENIB"),"drug_name"] <- "Sorafenib"

# Tissues
tissues <- unique(cellInfo(CCLE)$tissueid)

# Group them
featureinfo <- featureInfo(CCLE,"isoforms")
# for(i in 1:nrow(drugs)){
#   drug <- as.character(drugs[i,"drug_name"])
#   gene <- as.character(drugs[i,"gene_name"])
#   sensitivities <- NULL
#   sensitivity.output <- paste0(Output,"sensitivity/",gene,"_",toupper(drug),"/")
#   dir.create(sensitivity.output)
#   for(tissue in tissues){
#     if(is.na(tissue)) next;
#     print(tissue)
#     feature <- fNames(CCLE, "isoforms")[which(grepl(paste0("^",gene,"-"),featureinfo$TranscriptName))]
#     drug.sensitivity <- drugSensitivitySig(pSet = CCLE,
#                                            mDataType="isoforms",
#                                            drugs=drug,
#                                            features=feature,
#                                            tissues=tissue,
#                                            sensitivity.measure="auc_published",
#                                            molecular.summary.stat="median",
#                                            sensitivity.summary.stat="median",
#                                            nthread = 5,
#                                            verbose=TRUE)
#     sensitivity <- as.data.frame(drug.sensitivity)
#     idtoname <- featureinfo[,c(4,6)]
#     sensitivity <- merge(sensitivity,idtoname,by="row.names")
#     row.names(sensitivity) <- sensitivity$TranscriptName
#     sensitivity$Row.names <- NULL
#     sensitivity$EnsemblTranscriptId <- NULL
#     sensitivity$TranscriptName <- NULL
#     write.csv(sensitivity,file=paste0(sensitivity.output,gene,"_",toupper(drug),"_",tissue,"_auc_published.csv"))
#     
#     # sensitivity[which(sensitivity[,8] > 0.05),1] <- 0
#     sens <- as.data.frame(sensitivity[,1])
#     colnames(sens) <- tissue
#     if(is.null(sensitivities)){
#       sensitivities <- sens
#     }else{
#       sensitivities <- cbind(sensitivities,sens)
#     }
#     row.names(sensitivities) <- row.names(sensitivity)
#   }
#   sensitivities[is.na(sensitivities)] <- 0
#   write.csv(sensitivities,file=paste0(sensitivity.output,gene,"_",toupper(drug),"_all_tissue_auc_published.csv"))
  # 
  # # Heatmap
  # library(gplots)
  # filename = paste0(sensitivity.output,gene,"_",toupper(drug),"_all_tissue_auc_published.png")
  # png(file=filename,width = 600, height = 1600)
  # print(heatmap.2(as.matrix(t(sensitivities)), 
  #           col=bluered(75),
  #           Rowv = F,
  #           Colv = F,
  #           scale="column",
  #           trace="none",
  #           colsep = 1:ncol(sensitivities),
  #           rowsep = 1:nrow(sensitivities),
  #           sepwidth = c(0.001,0.001),
  #           sepcolor = 'white',
  #           # cexRow = 0.5,
  #           # cexCol = 0.8,
  #           dendrogram = "none",
  #           lhei = c(1,8)
  #           # lwid = c(0.5,5)
  #           ))
  # dev.off()
# }








#### Individual examples
drug <- "Crizotinib"
gene <- "MET"
featureinfo <- featureInfo(CCLE,"isoforms")
tissues <- unique(cellInfo(CCLE)$tissueid)
sensitivities <- NULL
pvals <- NULL
sensitivity.output <- paste0(Output,"sensitivity/",gene,"_",toupper(drug),"/")
dir.create(sensitivity.output)
for(tissue in tissues){
  if(is.na(tissue)) next;
  print(tissue)
  feature <- fNames(CCLE, "isoforms")[which(grepl(paste0("^",gene,"-"),featureinfo$TranscriptName))]
  drug.sensitivity <- drugSensitivitySig(pSet = CCLE,
                                         mDataType="isoforms",
                                         drugs=drug,
                                         features=feature,
                                         tissues=tissue,
                                         sensitivity.measure="auc_recomputed",
                                         molecular.summary.stat="median",
                                         sensitivity.summary.stat="median",
                                         nthread = 5,
                                         verbose=TRUE)
  sensitivity <- as.data.frame(drug.sensitivity)
  idtoname <- featureinfo[,c(4,6)]
  sensitivity <- merge(sensitivity,idtoname,by="row.names")
  row.names(sensitivity) <- sensitivity$TranscriptName
  sensitivity$Row.names <- NULL
  sensitivity$EnsemblTranscriptId <- NULL
  sensitivity$TranscriptName <- NULL
  write.csv(sensitivity,file=paste0(sensitivity.output,gene,"_",toupper(drug),"_",tissue,"_auc_recomputed.csv"))
  
  # sensitivity[which(sensitivity[,8] > 0.05),1] <- 0
  sens <- as.data.frame(sensitivity[,1])
  pval <- as.data.frame(sensitivity[,8])
  colnames(sens) <- tissue
  colnames(pval) <- tissue
  if(is.null(sensitivities)){
    sensitivities <- sens
    pvals <- pval
  }else{
    sensitivities <- cbind(sensitivities,sens)
    pvals <- cbind(pvals,pval)
  }
  row.names(sensitivities) <- row.names(sensitivity)
  row.names(pvals) <- row.names(sensitivity)
}
# sensitivities[is.na(sensitivities)] <- 0
write.csv(sensitivities,file=paste0(sensitivity.output,gene,"_",toupper(drug),"_all_tissue_auc_recomputed.csv"))
write.csv(pvals,file=paste0(sensitivity.output,gene,"_",toupper(drug),"_all_tissue_auc_recomputed_pval.csv"))

## Heatmap
sensitivities <- read.csv(file=paste0(sensitivity.output,gene,"_",toupper(drug),"_all_tissue_auc_recomputed.csv"),header=TRUE,row.names=1)
pvals <- read.csv(file=paste0(sensitivity.output,gene,"_",toupper(drug),"_all_tissue_auc_recomputed_pval.csv"),header=TRUE,row.names=1)
# sensitivities[sensitivities == 0] <- NA
sensitivities[pvals > 0.05] <- NA
sensitivities.melted <- melt(as.matrix(sensitivities),na.rm = FALSE)
p <- ggplot(sensitivities.melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  coord_equal() +
  scale_fill_gradientn(
                      # mid = "white",
                      # high = "indianred3",
                      # midpoint = 0,
                      na.value = "grey90",
                      limits = c(-1,1),
                      colours = c("steelblue","white","indianred3"),
                      breaks = c(-1,0,1),
                      labels = format(c(-1,0,1))) +
  ylab('Tissues') +
  xlab('Isoforms') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0(sensitivity.output,gene,"_",toupper(drug),"_all_tissue_auc_recomputed_heatmap.png"),p)