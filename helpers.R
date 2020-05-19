############## :::: HELPER FUNCTIONS FOR ISOFORM DRUG TARGETS PROJECT :::: ##################
extractBindingResidues <- function(res){
  # res is a vector containing binding residues information
  # can be from one protein binding or a column of protein bindings
  if(!require("stringr")){
    install.packages("stringr")
    library(stringr)
  }
  res <- strsplit(as.character(res),split=" ")
  res.num <- str_extract_all(res, "[0-9]+")
  res.chr <- str_extract_all(res, "[A-Z]+")
  
  # res.num <- sapply(res.num,as.numeric)
  
  res.list <- list()
  for(i in 1:length(res.num)){
    # print(i)
    res.list[[i]] <- list()
    for(j in 1:length(res.num[[i]])){
      res.list[[i]][[res.num[[i]][[j]]]] <- res.chr[[i]][[j]]
    }
  }
  return(res.list)
}




collapseBindingResidues <- function(object,res,...){
  # Input object needs to be a data frame
  
  ## Collapse binding residues (multiple versions of binding residues available, take union of all the residues)
  object$Binding_site_residues <- res
  object$id <- paste(object$ensembl_id,object$ChEMBL_id,sep="_")
  object$number <- row.names(object)
  object.grouped <- object %>%
    group_by(id) %>%
    split(.,.$id)
  
  new.object <- NULL
  for(gene.drug in object.grouped){
    # if(length(Reduce(intersect,gene.drug$Binding_site_residues)) > 0){ # Do not have conflicting binding residues/numbering
    collapsed.res <- list(Reduce(modifyList,gene.drug$Binding_site_residues))
    gene.drug$Binding_site_residues <- collapsed.res
    new.object <- rbind(new.object,gene.drug)
    # }
  }
  
  ## Sort new object based on original object
  new.object <- new.object[order(match(new.object$number,object$number)),]
  new.res <- new.object$Binding_site_residues
  rm(object.grouped)
  return(new.res)
}

makeBindingPocket <- function(object,...){
  # Input object can be the data frame with binding residues information
  # Or can be individual bindings (rows) within the dataframe.
  
  if(is.data.frame(object)){
    res <- object[,"Binding_site_residues"]
  }else{
    res <- object
  }
  
  res <- extractBindingResidues(res)
  res <- collapseBindingResidues(object,res)
  
  all.bp <- list()
  for(i in 1:length(res)){
    # print(i)
    
    # Sort residues by index, in case some are not sorted and create problems
    res[[i]] <- res[[i]][order(as.numeric(names(unlist(res[[i]]))))]
    res.chr <- unlist(res[[i]])
    res.num <- as.numeric(names(unlist(res[[i]])))
    
    if(anyDuplicated(res.num)){
      # In cases of same residue number but different aa
      # Keep first but remove second
      duplicated.index <- match(TRUE,duplicated(res.num))
      res.num <- res.num[-duplicated.index] 
      res.chr <- res.chr[-duplicated.index] 
    }
    
    gap <- diff(res.num) - 1 # Compute gap between two predicted residues
    gap <- gap[gap >= 0]
    # print(gap)
    bp <- res.chr
    
    if(length(gap) == 0){ # In case of empty row
      next
    }
    
    for(j in 1:length(gap)){
      if(j == 1){
        bp <- append(bp,rep("-",gap[j]),j)
      }else{
        bp <- append(bp,rep("-",gap[j]),res.num[j]-res.num[1]+1)
      }
    }
    all.bp[i] <- paste(bp,collapse="")
  }
  
  return(all.bp)
}




addBindingSeqRow <- function(df, .filter=FALSE)
{
  if(.filter){
    df.grouped <- df %>%
      mutate(id=paste(ensembl_id,drug_name,sep="_")) %>%
      group_by(id) %>%
      filter(n_distinct(protein_sequence) > 1) %>% # Filtered genes that only have 1 protein isoform
      filter(n_distinct(protein_id) > 1) %>% # Filtered genes that only have 1 protein isoform
      split(.,.$id)
  }else{
    df.grouped <- df %>%
      mutate(id=paste(ensembl_id,drug_name,sep="_")) %>%
      group_by(id) %>%
      split(.,.$id)
  }
  
  binding.df <- NULL
  for (binding in df.grouped)
  {
    # print(binding2)
    # print(unique(binding$id))
    binding.seq <- unique(binding$binding_pocket)
    binding$binding_pocket <- NULL
    binding.df1 <- binding[1,]
    binding.df1$protein_sequence <- binding.seq
    binding.df1$protein_id <- binding.df1$drug_name
    binding.df1$protein_name <- binding.df1$drug_name
    binding.df1$transcript_id <- binding.df1$drug_name
    binding <- rbind(binding.df1,binding)
    binding.df <- rbind(binding.df,binding)
  }
  binding.df$id <- NULL
  return(binding.df)
}




makeFastaFiles <- function(df,dir, .filter=FALSE)
{
  df <- as.data.frame(df)
  if(.filter){
    df.grouped <- df %>%
      mutate(id=paste(ensembl_id,drug_name,sep="_")) %>%
      group_by(id) %>%
      filter(n_distinct(protein_sequence) > 1) %>% # Filtered genes that only have 1 protein isoform
      filter(n_distinct(protein_id) > 1) %>% # Filtered genes that only have 1 protein isoform
      split(.,.$id)
  }else{
    df.grouped <- df %>%
      mutate(id=paste(ensembl_id,drug_name,sep="_")) %>%
      group_by(id) %>%
      # filter(n_distinct(protein_id) > 1) %>% # Filtered genes that only have 1 protein isoform
      split(.,.$id)
  }
 
  for (binding in df.grouped)
  {
    require(seqinr)
    write.fasta(sequences = as.list(binding$protein_sequence),
                names = as.list(paste(unique(binding$transcript_id),
                                      # "|",unique(interact2$gene_symbol),
                                      # "|",unique(interact2$Generic_name),
                                      sep="")),
                file.out = paste0(dir,unique(binding$gene_name),"_",unique(binding$drug_name),".fasta"),
                open = "w")
  }
}