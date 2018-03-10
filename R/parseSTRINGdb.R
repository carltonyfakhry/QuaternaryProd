# library(org.Hs.eg.db)
# library(stringr)
# library(readr)
# 
# # Function to parse the StringDB network
# parseStringDB <- function(){
# 
#   # Get the full file name containing the STRINGdb Rels
#   ff <- system.file("extdata", "9606.protein.actions.v10.txt.gz", package="QuaternaryProd")
#   all_rels <- read_tsv(gzfile(ff), col_names = TRUE)
# 
#   # Set new names for columns
#   names(all_rels) <- c("srcuid", "trguid", "mode", "action", "direction","score")
#   Rels <- all_rels[, c("srcuid", "trguid", "mode", "direction")]
# 
#   # Get all rows with causal Rels
#   Rels <- Rels[Rels$mode %in% c("activation", "inhibition","expression"),]
# 
#   # Get causal Rels where direction is not specified, and consider reversed
#   # direction of causality as a valid causal relation
#   Bidirectional <- Rels[Rels$direction == 0 , c("trguid", "srcuid", "mode", "direction")]
#   names(Bidirectional) <- c("srcuid", "trguid", "mode", "direction")
#   Rels <- bind_rows(Rels, Bidirectional)
#   Rels$direction <- NULL
#   Rels <- unique(Rels)
# 
#   # Rename activation as increases, inhibition as decreases, expression
#   # as regulates
#   Rels$mode <- sub("activation", "increases", Rels$mode)
#   Rels$mode <- sub("inhibition", "decreases", Rels$mode)
#   Rels$mode <- sub("expression", "regulates", Rels$mode)
#   Rels <- unique(Rels)
# 
#   # Get all unique protein ensemble ids in the causal network
#   allEns <- unique(c(Rels$srcuid, Rels$trguid))
# 
#   # Get the entities data frame
#   ff2 <- system.file("extdata", "ENSPid_to_Entrezid_v10.csv", package="QuaternaryProd")
#   Ensemble2Entrez = read.table(ff2, sep = ",", header = T, stringsAsFactors = F)
# 
#   Ents <- data.frame(id = Ensemble2Entrez$Entrez_Gene_ID, ensembleid = Ensemble2Entrez$STRING_Locus_ID, stringsAsFactors = F)
#   map = org.Hs.egSYMBOL
#   symbol <- unlist(mget(as.character(Ents$id), map, ifnotfound=NA))
#   symbol[is.na(symbol)] <- "No-Symbol"
#   Ents["symbol"] <- symbol
# 
#   # Keep entities in the network
#   Ents <- Ents[Ents$ensembleid %in% allEns,]
# 
#   # Remove ensemble ids in entities with duplicated entrez id
#   Ents <- Ents[!duplicated(Ents$id),]
# 
#   allEns2 <- allEns[which(!(allEns %in% Ents$ensembleid))]
#   map <- org.Hs.egENSEMBLPROT2EG
#   x <- unlist(mget(str_replace(allEns2, "9606.",""), map, ifnotfound=NA));
#   Ents2 <- data.frame(ensembleid=paste("9606.", names(x), sep=""), id= as.integer(x), stringsAsFactors = F)
#   Ents2$id[is.na(Ents2$id)] <- -1
#   Ents2$symbol <- unlist(mget(as.character(Ents2$id), org.Hs.egSYMBOL, ifnotfound=NA))
#   Ents2 <- Ents2[Ents2$id != -1,]
#   Ents2 <- Ents2[!duplicated(Ents2$id),]
#   Ents2 <- Ents2[,c("id","ensembleid","symbol")]
#   Ents <- rbind(Ents, Ents2)
#   Ents$id[duplicated(Ents$id)] <- "No-EntrezId"
#   Ents$symbol[duplicated(Ents$symbol)] <- "No-Symbol"
# 
#   allEns2 <- allEns[which(!(allEns %in% Ents$ensembleid))]
#   Ents2 <- data.frame(id = "No-EntrezId", ensembleid = allEns2, symbol = "No-Symbol")
#   Ents <- rbind(Ents, Ents2)
#   Ents <- cbind(uid = 1:nrow(Ents), Ents)
#   Ents <- Ents[!duplicated(Ents$ensembleid),]
# 
#   # Get all unique Rels
#   Rels <- Rels[!(Rels$trguid %in% allEns2),]
#   Rels <- Rels[which(Rels$srcuid %in% Ents$ensembleid & Rels$trguid %in% Ents$ensembleid),]
#   Rels2 <- left_join(Rels, Ents, by = c("srcuid" = "ensembleid")) %>% dplyr::select(srcuid = uid, trguid, mode)
#   Rels2 <- left_join(Rels2, Ents, by = c("trguid" = "ensembleid")) %>% dplyr::select(srcuid, trguid = uid, mode)
#   Rels <- data.frame(Rels2)
#   rm(Rels2)
#   Rels$mode[which(Rels$mode == "increases")] <- 1
#   Rels$mode[which(Rels$mode == "decreases")] <- -1
#   Rels$mode[which(Rels$mode == "regulates")] <- 0
#   Rels$mode <- as.integer(Rels$mode)
#   Rels <- Rels[order(Rels$srcuid, Rels$trguid),]
#   Rels <- unique(Rels)
# 
#   # Leave source proteins which contain at least 10 edges
#   sufficientRels <- group_by(Rels, srcuid) %>% summarise(count=n())
#   sufficientRels <- sufficientRels %>% filter(count > 10)
#   Rels <- Rels %>% filter(srcuid %in% sufficientRels$srcuid)
#   Rels <- Rels[order(Rels$srcuid, Rels$trguid),]
# 
#   # Handle ambiguities if source node affects target node in more than one way
#   temp.rels <- Rels[which(Rels$mode %in% c(1,-1)),]
#   temp.rels <- temp.rels[,c("srcuid","trguid")]
#   indices1 <- which(duplicated(temp.rels))
#   for(i in indices1){
#     indices2 <- which(Rels$srcuid == temp.rels$srcuid[i] & Rels$trguid == temp.rels$trguid[i])
#     if(all(c(-1,1) %in% unique(Rels$mode[indices2]))){
#       Rels$mode[indices2] = 0
#     }else if(1 %in% unique(Rels$mode[indices2])){
#       Rels$mode[indices2] = 1
#     }else if(-1 %in% unique(Rels$mode[indices2])){
#       Rels$mode[indices2] = -1
#     }else{
#       Rels$mode[indices2] = 0
#     }
#   }
#   Rels <- Rels[!duplicated(Rels),]
#   Rels <- Rels[which(Rels$srcuid != Rels$trguid),]
# 
#   write.table(Rels, "/home/kolonel/Programs/QuaternaryProd/inst/extdata/StringRels.dat", col.names = T, row.names = F)
#   write.table(Ents, "/home/kolonel/Programs/QuaternaryProd/inst/extdata/StringEnts.dat", col.names = T, row.names = F)
# 
 Rels = read.table("~/Programs/QuaternaryProd/inst/extdata/StringRels.dat", header = T)
 Ents = read.table("~/Programs/QuaternaryProd/inst/extdata/StringEnts.dat", header = T)
 
 # Proteins to be tested
 u.hyps <- unique(Rels$srcuid)

 # For each protein, get its children
 child.uid <- lapply(u.hyps, function(x, Rels) Rels$trguid[which(Rels$srcuid == x)], Rels = Rels)

 # For each protein, Find if its children are upregulated or downregulated
 child.sgn <- lapply(u.hyps, function(x, Rels) ifelse(Rels$mode[which(Rels$srcuid == x)] == 1,
                                                           1, ifelse(Rels$mode[which(Rels$srcuid == x)] == -1, -1, 0)), Rels = Rels)

 write_yaml(u.hyps, file='/home/kolonel/Programs/QuaternaryProd/inst/extdata/u.hyps.yaml')
 write_yaml(child.uid, file='/home/kolonel/Programs/QuaternaryProd/inst/extdata/child.uid.yaml')
 write_yaml(child.sgn, file='/home/kolonel/Programs/QuaternaryProd/inst/extdata/child.sgn.yaml')

# }
