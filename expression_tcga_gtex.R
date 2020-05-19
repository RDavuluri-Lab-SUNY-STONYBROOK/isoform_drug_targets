############ :::: DRUG TARGET EXPRESSION TCGA GTEX :::: #############
library(data.table)
library(stringdist)
library(dplyr)

GTEX_ONLY = TRUE

# data
tcga_gtex <- fread("/projects/b1017/shared/TCGA_data/pan_cancer/TcgaTargetGtex_rsem_isoform_tpm",sep = "\t")
tcga_gtex <- as.data.frame(tcga_gtex)
row.names(tcga_gtex) <- tcga_gtex$sample

# meta
all.meta <- fread("/projects/b1017/shared/TCGA_data/pan_cancer/TcgaTargetGTEX_phenotype.txt",sep = "\t")
# filter recurrent/metastatic
all.meta <- all.meta[which(all.meta$`_study` %in% c("TCGA","GTEX")),]
all.meta <- all.meta[which(all.meta$`_sample_type` %in% c("Primary Tumor",
                                                          "Additional - New Primary",
                                                          "Primary Blood Derived Cancer - Peripheral Blood",
                                                          "Normal Tissue",
                                                          "Solid Tissue Normal",
                                                          "Primary Solid Tumor",
                                                          "Primary Blood Derived Cancer - Bone Marrow"
                                                          )),]
all.meta$id <- ifelse(substring(all.meta$sample,1,4) == "TCGA",
                      substring(all.meta$sample,1,15),
                      all.meta$sample)


# subset
tcga_gtex_tpm <- tcga_gtex[,which(colnames(tcga_gtex) %in% all.meta$id)]

rm(tcga_gtex)


# separate tumor and normal
erroneous_entries <- all.meta[which(all.meta$`_primary_site` == ""),"primary disease or tissue"]
erroneous_entries <- strsplit(erroneous_entries$`primary disease or tissue`," - ")
erroneous_entries <- sapply(erroneous_entries, function(x) x[[1]])
all.meta[which(all.meta$`_primary_site` == ""),"_primary_site"] <- erroneous_entries
meta.n <- all.meta[which(all.meta$`_sample_type` %in% c("Normal Tissue","Solid Tissue Normal"))]
meta.t <- all.meta[which(!(all.meta$sample) %in% meta.n$sample),]

tpm_normal <- tcga_gtex_tpm[,which(colnames(tcga_gtex_tpm) %in% meta.n$id)]
tpm_tumor <- tcga_gtex_tpm[,which(colnames(tcga_gtex_tpm) %in% meta.t$id)]

save(tpm_tumor,tpm_normal, meta.t,meta.n, file = "/projects/b1017/shared/TCGA_data/pan_cancer/TCGA_GTEx_rsem_isoform_tpm.RData")
rm(tcga_gtex_tpm)


# compute mean by each tissue
tpm_tumor_1 <- as.data.frame(t(tpm_tumor))
tpm_tumor_1 <- tpm_tumor_1[order(row.names(tpm_tumor_1)),]
meta.t <- meta.t[order(meta.t$sample),]
tpm_tumor_1$disease <- meta.t$`primary disease or tissue`

tpm_tumor_mean <- setDT(tpm_tumor_1)[, 
                   lapply(
                     .SD, 
                     function(x) 
                       sum(x)/length(x)
                   ), 
                   by=disease]

tpm_tumor_mean <- as.data.frame(tpm_tumor_mean)
row.names(tpm_tumor_mean) <- tpm_tumor_mean$disease
tpm_tumor_mean$disease <- NULL
tpm_tumor_mean_1 <- as.data.frame(t(tpm_tumor_mean))
write.csv(tpm_tumor_mean_1,file = "/projects/b1017/shared/TCGA_data/pan_cancer/TCGA_GTEx_rsem_isoform_tpm_tumor_mean.csv")


if(GTEX_ONLY){
  meta.n <- meta.n[which(meta.n$`_study` == "GTEX"),]
  tpm_normal <- tpm_normal[,which(colnames(tpm_normal) %in% meta.n$id)]
}

tpm_normal_1 <- as.data.frame(t(tpm_normal))
tpm_normal_1 <- tpm_normal_1[order(row.names(tpm_normal_1)),]
meta.n <- meta.n[order(meta.n$sample),]
tpm_normal_1$tissue <- meta.n$`_primary_site`

tpm_normal_mean <- setDT(tpm_normal_1)[, 
                                     lapply(
                                       .SD, 
                                       function(x) 
                                         sum(x)/length(x)
                                     ), 
                                     by=tissue]
tpm_normal_mean <- as.data.frame(tpm_normal_mean)
row.names(tpm_normal_mean) <- tpm_normal_mean$tissue
tpm_normal_mean$tissue <- NULL
tpm_normal_mean_1 <- as.data.frame(t(tpm_normal_mean))
write.csv(tpm_normal_mean_1,file = "/projects/b1017/shared/TCGA_data/pan_cancer/GTEx_rsem_isoform_tpm_normal_mean.csv")

# Some adjustments
colnames(tpm_tumor_mean_1) <- gsub("\\."," ",colnames(tpm_tumor_mean_1))
colnames(tpm_normal_mean_1) <- gsub("\\."," ",colnames(tpm_normal_mean_1))


# Map diseases to tissues
# d <- stringdistmatrix(colnames(tpm_tumor_mean_1),colnames(tpm_normal_mean_1),method = "lcs")
# row.names(d) <- colnames(tpm_tumor_mean_1)
# colnames(d) <- colnames(tpm_normal_mean_1)
# assigned.tissue <- apply(d,1,which.min)
# map <- data.frame(disease=colnames(tpm_tumor_mean_1),tissue=colnames(tpm_normal_mean_1)[assigned.tissue],
#                   stringsAsFactors = FALSE)
# 
# # Manually correct
# map[which(map[,"disease"] == "Glioblastoma Multiforme"),"tissue"] <- "Brain"
# map[which(map[,"disease"] == "Uterine Corpus Endometrioid Carcinoma"),"tissue"] <- "Uterus"
# map[which(map[,"disease"] == "Sarcoma"),"tissue"] <- "Muscle"
# map[which(map[,"disease"] == "Mesothelioma"),"tissue"] <- "Lung"
# map[which(map[,"disease"] == "Cholangiocarcinoma"),"tissue"] <- "Liver"
# map[which(map[,"disease"] == "Acute Myeloid Leukemia"),"tissue"] <- "Blood"
# map[which(map[,"disease"] == "Diffuse Large B Cell Lymphoma"),"tissue"] <- "Blood"
# map[which(map[,"disease"] == "Uveal Melanoma"),"tissue"] <- "Skin"
# # map[which(map[,"disease"] == "Pheochromocytoma Paraganglioma"),"tissue"] <- "Adrenal Gland"
# # Pheochromocytoma Paraganglioma
# # map[which(map[,"disease"] == "Cervical Endocervical Cancer"),"tissue"] <- "Cervix"
# 
# # map <- map[-c(13),]
# # map <- rbind(map,c("Cervical Endocervical Cancer","Cervix"))
# write.csv(map, file = "/projects/b1017/shared/TCGA_data/pan_cancer/TCGA_GTEx_tissue_mapping_GTExonly.csv")
map <- read.csv(file = "/projects/b1017/shared/TCGA_data/pan_cancer/TCGA_GTEx_tissue_mapping_GTExonly.csv")


# Calculate fold change
map <- map[order(map$disease),]
tpm_tumor_mean_1 <- tpm_tumor_mean_1[,order(colnames(tpm_tumor_mean_1))]
tpm_tumor_mean_1 <- tpm_tumor_mean_1[,-which(colnames(tpm_tumor_mean_1) %in% as.character(map[which(map$tissue == ""),"disease"]))]
tpm_normal_mean_2 <- tpm_normal_mean_1[,as.character(map[which(map$tissue != ""),"tissue"])]
tpm_fc <- tpm_tumor_mean_1 - tpm_normal_mean_2
# tpm_fc <- apply(tpm_tumor_mean_1,2, function(x){
#   disease <- colnames(x)
#   tissue <- map[which(map$disease == disease),"tissue"]
#   as.numeric(x) - tpm_normal_mean_1[,which(colnames(tpm_normal_mean_1) == tissue)]
# })
write.csv(tpm_fc,file = "/projects/b1017/shared/TCGA_data/pan_cancer/TCGA_GTEx_rsem_isoform_tpm_fc_GTExonly_corrected.csv")

rm(tpm_tumor,tpm_normal,tpm_tumor_mean,tpm_normal_mean,tpm_tumor_1,tpm_normal_1,tpm_tumor_mean_1,tpm_normal_mean_1,tpm_normal_mean_2)

####################### Count for individual drug targets ######################
load(file = "/projects/b1017/Jerry/Isoform_drug_targets/results/all/isoforms_binding_indicator.RData")
tpm_fc <- read.csv(file = "/projects/b1017/shared/TCGA_data/pan_cancer/TCGA_GTEx_rsem_isoform_tpm_fc_GTExonly_corrected.csv",row.names=1)
row.names(tpm_fc) <- gsub("(.*)\\.[0-9]","\\1",row.names(tpm_fc))

final_res <- as.data.frame(t(tpm_fc[0,]))

for(j in 1:length(res)){
  pair <- res[[j]]
  pairname <- names(res)[j]
  tpm <- tpm_fc[which(row.names(tpm_fc) %in% pair$isoforms),]
  tpm$isoforms <- row.names(tpm)
  tpm <- merge(tpm, pair, by = "isoforms")
  
  # if(strsplit(pairname,"_")[[1]][1] %in% c("ABL1","ABL2")){ # if ABL1 or ABL2, assign type III
  #   final_res[,pairname] <- "Type III"
  # }else{
    if(sum(tpm$ibind) == 1){ # only one isoform binds to the drug
      # look at expression of those isoforms not bound
      # if at least one upregulated, then Type II
      # otherwise specific
      tpm_1 <- tpm[tpm$ibind == 0,] # those not targeted
      for(i in 2:(ncol(tpm)-1)){
        if(sum(tpm_1[,i]>0) >= 1){ # at least 1 upregulated is not targeted
          final_res[i-1, pairname] <- "Type II"
        }else{
          final_res[i-1,pairname] <- "specific" # no missed isoforms
        }
      }
    }else{ # two or more isoforms binds to the drug
      
      tpm_2 <- tpm[tpm$ibind == 1,] # those targeted
      tpm_1 <- tpm[tpm$ibind == 0,] # those not targeted
      for(i in 2:(ncol(tpm)-1)){
        if(sum(tpm_2[,i]<0) >= 1){ # at least 1 downregulated is targeted
          final_res[i-1, pairname] <- "Type I"
          
          if(sum(tpm_1[,i]>0) >= 1){ # at least 1 upregulated is not targeted
            final_res[i-1, pairname] <- "Type III"
          }
        }else if(sum(tpm_1[,i]>0) >= 1){ # at least 1 upregulated is not targeted
          final_res[i-1, pairname] <- "Type II"
        }else{
          final_res[i-1,pairname] <- "specific" # targets all overexpressed isoforms
        }
      }
    }
  # }
}

final_res <- as.data.frame(t(final_res))

sum(final_res == "specific")/(nrow(final_res)*ncol(final_res)) # 0.2417 isoform specific
1- sum(final_res == "specific")/(nrow(final_res)*ncol(final_res)) # 0.7583 not specific

write.csv(final_res,file = "/projects/b1017/Jerry/Isoform_drug_targets/results/all/drug_isoform_specificity_by_cancer_GTExonly_corrected.csv")

## Heatmap
library(ggplot2)
library(reshape2)
library(RColorBrewer)
# library(data.table)
final_res$id <- row.names(final_res)
final_res_melted <- melt(final_res,id.vars="id")
final_res_melted$variable <- gsub("\\."," ",final_res_melted$variable)
p <- ggplot(final_res_melted, aes(x = variable, y = id, fill = value)) +
  geom_tile(alpha = 0.5, color = "white") +
  scale_fill_manual(
    # values=c("palegreen1","firebrick1","darkorange1","mediumpurple2"),
    values=brewer.pal(4,"YlGnBu"),
    aesthetics = c( "fill")) +
  ylab('Drug target interaction') +
  xlab('Disease') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
        axis.text.y = element_text(size=2),
        # axis.ticks=element_line(size=0.4)
        axis.ticks = element_blank()
        )
ggsave(filename = paste0("./results/all/drug_isoform_specificity_by_cancer_GTExonly_corrected.tiff"),p)


## Pie chart
library(dplyr)
library(scales)
final_res_melted %>%
  group_by(value) %>% 
  tally() %>%
  mutate(pct = n/sum(n)) %>%
  ggplot(aes(x="", y = pct, fill=value)) + 
    geom_bar(stat = "identity", position = position_fill(),color = "black", alpha = 0.5) +
    xlab("") + ylab("") +
    scale_fill_manual(
      values=brewer.pal(4,"YlGnBu"),
      aesthetics = c( "fill"))+
    theme_minimal() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_blank()) +
    coord_polar(theta="y") +
    geom_text(aes(label = percent(pct)), position = position_stack(vjust = 0.5))

ggsave(filename = paste0("./results/all/drug_isoform_specificity_by_cancer_GTExonly_pie_corrected.tiff"))


## Venn diagram
library(VennDiagram)

venn.diagram(list(`Type I`=row.names(final_res_melted[which(final_res_melted$value %in% c("Type I","Type III")),]),
                  `Type II`=row.names(final_res_melted[which(final_res_melted$value %in% c("Type II","Type III")),])), 
             fill=c("#A1DAB4", "#41B6C4"), 
             alpha=c(0.5,0.5), 
             cex=2,
             filename="./results/all/drug_isoform_specificity_by_cancer_GTExonly_venn_corrected.tiff")

