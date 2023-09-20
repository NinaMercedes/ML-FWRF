## Knowledge Graph Annotation
#Load libraries
library(dplyr,quietly = TRUE)
library(visNetwork, quietly = TRUE)
library(shiny,quietly = TRUE)
library(igraph, quietly=TRUE)
set.seed(10)

#Set options
option_list <- list(
  make_option(c("-g", "--genes_of_interest"), type="character", default=FALSE,
    help="all genes of interest (drug-specific)"),
  make_option(c("-k", "--known"), type="character", default=FALSE,
    help="known genes of interest"),
  make_option(c("-c", "--compensatory"), type="character", default=TRUE,
    help="compensatory genes of interest"),
  make_option(c("-d", "--drug_name"), type="character", default=FALSE,
    help="name of drug eg rifampicin"),
  make_option(c("-t", "--drug_type"), type="character", default=TRUE,
    help="eg aminoglycosides"), 
  make_option(c("-h", "--hl"), type="character", default=TRUE,
    help="highly likely threshold"),
 make_option(c("-u", "--ul"), type="character", default=TRUE,
    help="highly unlikely threshold"),
 make_option(c("-x", "--un"), type="character", default=TRUE,
    help="unlikely threshold"),
 make_option(c("-a", "--hgvs"), type="character", default=TRUE,
    help="hgvs annotation"),
 make_option(c("-l", "--lineage"), type="character", default=TRUE,
    help="lineage snp ids (from phylogenetic analysis"),
 make_option(c("-b", "--barcode"), type="character", default=TRUE,
    help="lineage barcode"),
 make_option(c("-p", "--tbdb"), type="character", default=TRUE,
    help="tbdb (TBProfiler database of known drug-resistant mutations)"),
 make_option(c("-s", "--string"), type="character", default=TRUE,
    help="string PPI network"),
 make_option(c("-i", "--interactions"), type="character", default=TRUE,
    help="interactions from MLFWRF"),
    )
#parse options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)

#### Set Parameters ####
genes_of_interest <- opt$genes_of_interest #vector
genes_of_interest_known <- opt$known #vector
genes_of_interest_compensatory <- opt$compensatory #vector
drug_name <- opt$drug
all_drug_names <- opt$drug_type
thresholds_hl <- opt$hl
thresholds_ul <- opt$ul
thresholds_u <- opt$un

#Read in annotation files and select candidate genes and phenotypes of interest
cand_gene <- c("embR","rrs", "fabG1","inhA","rpsA","tlyA","katG","pncA","kasA","panD","embC","embA","embB",NA,"gid", "gyrB","gyrA","rpoB","rpoC","rpsL","ahpC","eis" )
drug_labels <- c("rifampicin", "isoniazid", "ethambutol", "pyrazinamide", "streptomycin", "ofloxacin", "moxifloxacin", "amikacin", "kanamycin", "capreomycin")
hgvs <- read.csv(opt$hgvs
lineage_snp_id <- read.csv(opt$lineage, header=TRUE)
lineage_barcode <- read.csv(opt$barcode, header=TRUE)
tbdb <- read.csv("tbdb.csv" header=TRUE)
string_interactions_short <- read.delim(opt$string, header=FALSE, comment.char="#")
interactions2 <- read.csv(opt$interactions)

#Process files
hgvs <- hgvs %>%   filter(Effect=="MODIFIER" & Gene %in% c("rrs")|!Effect=="MODIFIER" & !Gene  %in% c("rrs"))
hgvs[is.na(hgvs[,c("Protein_Change")]),c("Protein_Change")]  <- hgvs[is.na(hgvs[,c("Protein_Change")]),c("NA_Change")]  
hgvs$ID <- paste0(hgvs$ID,"_",hgvs$Position,"_",hgvs$Alt_Allele)
hgvs <- unique(hgvs)
lineage_barcode$id <- paste0("Chromosome_",lineage_barcode$Pos,"_",lineage_barcode$Alt)
#Read in interactions data frame and annotate with hgvs nomenclature
interactions2 <- interactions2[,-1]
interactions2 <- data.frame(apply(interactions2, 2, function(x) gsub(" ", "", x)))
interactions2 <- unique(interactions2)
interactions2 <- left_join(interactions2, hgvs, by= c("Mutation_1"="ID"),relationship ="many-to-many")
interactions2 <- left_join(interactions2, hgvs, by= c("Mutation_2"="ID"),relationship ="many-to-many")
interactions2[grepl("fs", interactions2$Protein_Change.x),c("Protein_Change.x")] <- interactions2[grepl("fs", interactions2$Protein_Change.x),c("NA_Change.x")]
interactions2[grepl("fs", interactions2$Protein_Change.y),c("Protein_Change.y")] <- interactions2[grepl("fs", interactions2$Protein_Change.y),c("NA_Change.y")]
interactions2[grepl("\\?", interactions2$Protein_Change.x),c("Protein_Change.x")] <- interactions2[grepl("\\?", interactions2$Protein_Change.x),c("NA_Change.x")]
interactions2[grepl("\\?", interactions2$Protein_Change.y),c("Protein_Change.y")] <- interactions2[grepl("\\?", interactions2$Protein_Change.y),c("NA_Change.y")]
interactions2 <- interactions2 %>% filter(Gene.x %in% cand_gene & Gene.y %in% cand_gene)
#Integrate domain knowledge from tbdb (TBProfiler) and calculate knowledge score
interactions2$known_1 <- NA
interactions2$coocurring_1 <- NA
interactions2$known_2 <- NA
interactions2$coocurring_2 <- NA
interactions2$total_effect1 <- NA
interactions2$total_effect2 <- NA
interactions2$candgene1 <- NA
interactions2$candgene2 <- NA
known<- tbdb %>% filter(Drug==paste0(drug_name))
coocurring <- tbdb %>% filter(!Drug %in% all_drug_names)
interactions2$known_1 <- ifelse(interactions2$Protein_Change.x %in% known$Mutation & interactions2$Gene.x %in% c(genes_of_interest), 1,0)
interactions2$known_2 <- ifelse(interactions2$Protein_Change.y %in% known$Mutation & interactions2$Gene.y %in% c(genes_of_interest), 1,0)
interactions2$coocurring_1 <- ifelse(interactions2$Protein_Change.x %in% coocurring$Mutation & !interactions2$Gene.x %in% c(genes_of_interest), 1,0)
interactions2$coocurring_2 <- ifelse(interactions2$Protein_Change.y %in% coocurring$Mutation & !interactions2$Gene.y %in% c(genes_of_interest), 1,0)
interactions2$lineage_1 <- ifelse(interactions2$Mutation_1 %in% lineage_barcode$id, 1,0)
interactions2$lineage_2 <- ifelse(interactions2$Mutation_2 %in% lineage_barcode$id, 1,0)
interactions2$total_effect1 <- ifelse(interactions2$Effect.x %in% c("HIGH", "MODERATE"), 1,0)
interactions2$total_effect2<- ifelse(interactions2$Effect.y %in% c("HIGH", "MODERATE"), 1,0)
interactions2$candgene1 <- ifelse(interactions2$Gene.x %in% c(genes_of_interest), 1,0)
interactions2$candgene2<- ifelse(interactions2$Gene.y %in% c(genes_of_interest), 1,0)
interactions2$score_1 <- NA
interactions2[,25:32] <- apply(interactions2[,25:32],2, function(x) as.numeric(as.character(x)))
interactions2 <- data.frame(interactions2)
interactions2$score_1 <-  interactions2$known_1 + interactions2$known_2 - interactions2$coocurring_1 - interactions2$coocurring_2 + interactions2$total_effect1 + interactions2$total_effect2 + interactions2$candgene1 + interactions2$candgene2 - interactions2$lineage_1 - interactions2$lineage_2
interactions2$score_2 <- interactions2$score_1 * as.numeric(as.character(interactions2$coeff))
interactions2$ks_1 <- NA
#Interpret knowledge score
for (i in 1: nrow(interactions2)){
  if (interactions2$score_2[i]>=thresholds_hl){
    interactions2$ks_1[i] <- "Highly likely"
  }
  if (interactions2$score_2[i]<=thresholds_ul){
    interactions2$ks_1[i] <- "Highly unlikely"
  }
  if (interactions2$score_2[i]>thresholds_ul & interactions2$score_2[i]<thresholds_u){
    interactions2$ks_1[i] <- "Unlikely"
  }
  if (interactions2$score_2[i]>=thresholds_u & interactions2$score_2[i]<thresholds_hl){
    interactions2$ks_1[i] <- "Likely"
  }
}

#Set colour parameters and naming convention
interactions2$ec_1 <- "grey"
interactions2$ec_2 <- "black"
interactions2$gen_int <- 1
interactions2$ks_2 <- interactions2$ks_1
interactions2$name.x <- paste0(interactions2$Gene.x,": ",interactions2$Protein_Change.x)
interactions2$name.y <- paste0(interactions2$Gene.y,": ",interactions2$Protein_Change.y)
#Integrate PPI network
string_interactions_short <- string_interactions_short[,c(1:2,13)]
colnames(string_interactions_short) <- c("P1", "P2", "coeff")
string_interactions_short$candgene1 <- ifelse(string_interactions_short$P1 %in% c(genes_of_interest), 1,0)
string_interactions_short$candgene2<- ifelse(string_interactions_short$P2 %in% c(genes_of_interest), 1,0)
string_interactions_short$known_1 <- ifelse(string_interactions_short$P1 %in% c(genes_of_interest), 1,0)
string_interactions_short$known_2<- ifelse(string_interactions_short$P2 %in% c(genes_of_interest), 1,0)
string_interactions_short$coocurring_1 <- ifelse(string_interactions_short$P1 %in% c(genes_of_interest), 0,1)
string_interactions_short$coocurring_2<- ifelse(string_interactions_short$P2 %in% c(genes_of_interest), 0,1)
string_interactions_short$ks_1 <- "Not Applicable"
string_interactions_short <- unique(string_interactions_short)
string_interactions_short$name.x <- string_interactions_short$P1
string_interactions_short$name.y <- string_interactions_short$P2
string_interactions_short$ec_1 <- "black"
string_interactions_short$group_1 <- paste0(string_interactions_short$P1,"_PPI")
string_interactions_short$group_2 <- paste0(string_interactions_short$P2,"_PPI")

## Let's make nodes data.frame
nodes <- data.frame(id=c(interactions2$Mutation_1, interactions2$Mutation_2, string_interactions_short$P1, string_interactions_short$P2),  
                    gene = c(interactions2$Gene.x, interactions2$Gene.y, string_interactions_short$P1, string_interactions_short$P2), 
                    candgene= c(interactions2$candgene1, interactions2$candgene2, string_interactions_short$candgene1, string_interactions_short$candgene2), 
                    known= c(interactions2$known_1, interactions2$known_2, string_interactions_short$known_1, string_interactions_short$known_2), 
                    coocurring=c(interactions2$coocurring_1, interactions2$coocurring_2, string_interactions_short$coocurring_1, string_interactions_short$coocurring_2),
                    name= c(interactions2$name.x, interactions2$name.y, string_interactions_short$name.x, string_interactions_short$name.y),
                    group=c(interactions2$Gene.x, interactions2$Gene.y, string_interactions_short$group_1, string_interactions_short$group_2))
 nodes <- unique( nodes)
 nodes$title <- paste0("<p><b>",  nodes$name,"</b><br> </p>")
 nodes$type <- NA
for (i in 1: nrow( nodes)){
  if ( nodes$id[i] %in% c(genes_of_interest_known)){
     nodes$type[i] <- "Known"
  }
  if ( nodes$id[i] %in% c(genes_of_interest_compensatory)){
     nodes$type[i] <- "Compensatory"
  }
  if (! nodes$id[i] %in% c(genes_of_interest)){
     nodes$type[i] <- "Co-occurring"
  }
  if ( nodes$known[i]>0){
    if ( nodes$gene[i] %in% c(genes_of_interest_compensatory)){
       nodes$type[i] <- "Compensatory"
    } else {
       nodes$type[i] <- "Known"
    }
  }
  if ( nodes$coocurring[i]>0){
     nodes$type[i] <- "Co-occurring"
  }
  if ( nodes$known[i]<1 &  nodes$coocurring[i]<1){
    if ( nodes$id[i] %in% lineage_barcode$id){
       nodes$type[i] <- "Lineage"
    } else {
       nodes$type[i] <- "Other"
    }
  }
}

 nodes$shape <- NA
 nodes$shadow <- NA
 nodes$border <- NA
 nodes$size <- NA
 nodes$label <- NA
 nodes$value <- NA
for (i in 1:nrow( nodes)){
  if( nodes$id[i] %in% unique( nodes$gene)){
     nodes$shape[i]<- "ellipse"
     nodes$shadow[i]<- TRUE
     nodes$border[i]<- "black"
     nodes$value[i]<- 50000
     nodes$size[i]<- 100
     nodes$label[i]<-  nodes$id[i]
  } else {
     nodes$shape[i]<- "dot"
     nodes$shadow[i]<- FALSE
     nodes$label[i]<- " "
     nodes$size[i]<- 50
  }
}

for (i in 1:nrow( nodes)){ 
  if( nodes$id[i] %in% drug_labels){
     nodes$shape[i]<- "box"
     nodes$shadow[i]<- TRUE
     nodes$border[i]<- "black"
     nodes$value[i]<- 50000
     nodes$size[i]<- 100
     nodes$label[i]<-  nodes$id[i]
     nodes$type[i] <- "Label"
     nodes$group[i] <-  nodes$id[i]
  }
}


#Make edges data frame
interactions2$width1 <- 0.2
interactions2$width2 <- 0.2
interactions2$dash1<- FALSE
interactions2$dash2 <- TRUE
string_interactions_short$width1 <- 0.2
string_interactions_short$ec_1 <- "black"
string_interactions_short$dash <- TRUE
edges <- data.frame(from = c(interactions2$Mutation_1,interactions2$Mutation_1, interactions2$Mutation_2, string_interactions_short$P1), 
                    to = c(interactions2$Mutation_2, interactions2$Gene.x, interactions2$Gene.y, string_interactions_short$P2), 
                    knowledge_score = c(interactions2$ks_1, interactions2$ks_2, interactions2$ks_2, string_interactions_short$ks_1), 
                    coeff= c(interactions2$coeff, interactions2$gen_int, interactions2$gen_int, string_interactions_short$coeff),
                    edge_colour= c(interactions2$ec_1, interactions2$ec_2, interactions2$ec_2, string_interactions_short$ec_1),
                    value= c(interactions2$width1, interactions2$width2,interactions2$width2, string_interactions_short$width1),
                    dashes=c(interactions2$dash1, interactions2$dash2, interactions2$dash2,string_interactions_short$dash))

edges$name <- paste0(edges$to,"_",edges$from)
edges$coeff2 <- NA
edges$coeff2 <- ifelse(edges$coeff>=0.7, "Yes", "No")
edges$value <- NULL
edges <- unique(edges)
interactions2$dash1<- FALSE
interactions2$dash2 <- TRUE
string_interactions_short$width1 <- 0.2
string_interactions_short$ec_1 <- "black"
string_interactions_short$dash <- TRUE
edges_capreomycin <- data.frame(from = c(interactions2$Mutation_1,interactions2$Mutation_1, interactions2$Mutation_2, string_interactions_short$P1), 
                    to = c(interactions2$Mutation_2, interactions2$Gene.x, interactions2$Gene.y, string_interactions_short$P2), 
                    knowledge_score = c(interactions2$ks_1, interactions2$ks_2, interactions2$ks_2, string_interactions_short$ks_1), 
                    coeff= c(interactions2$coeff, interactions2$gen_int, interactions2$gen_int, string_interactions_short$coeff),
                    edge_colour= c(interactions2$ec_1, interactions2$ec_2, interactions2$ec_2, string_interactions_short$ec_1),
                    value= c(interactions2$width1, interactions2$width2,interactions2$width2, string_interactions_short$width1),
                    dashes=c(interactions2$dash1, interactions2$dash2, interactions2$dash2,string_interactions_short$dash))

edges_capreomycin$name <- paste0(edges_capreomycin$to,"_",edges_capreomycin$from)
edges_capreomycin$coeff2 <- NA
edges_capreomycin$coeff2 <- ifelse(edges_capreomycin$coeff>=0.7, "Yes", "No")
edges_capreomycin$value <- NULL
edges_capreomycin <- unique(edges_capreomycin)

#write nodes and edges to file
write.csv(edges, paste0(opt$drug_name,"_edges.csv", row.names=FALSE)
write.csv(nodes, paste0(opt$drug_name,"_edges.csv", row.names=FALSE)