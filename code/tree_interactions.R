library(tidyr)
library(data.table)
set.seed(11)
option_list <- list(
  make_option(c("-d", "--drug"), type="character", default=FALSE,
    help="drug name")
    )
#parse options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)
# Calculate the depth of each node in a single tree obtained from a forest with ranger::treeInfo
calculate_tree_depth_ranger <- function(frame){
  if(!all(c("rightChild", "leftChild") %in% names(frame))){
    stop("The data frame has to contain columns called 'rightChild' and 'leftChild'!
         It should be a product of the function ranger::treeInfo().")
  }
  frame$depth <- NA
  frame$depth[1] <- 0
  for(i in 2:nrow(frame)){
    frame[i, "depth"] <-
      frame[(!is.na(frame[, "leftChild"]) & frame[, "leftChild"] == frame[i, "nodeID"]) |
              (!is.na(frame[, "rightChild"]) & frame[, "rightChild"] == frame[i, "nodeID"]), "depth"] + 1
  }
  return(frame)
}
get_interactions <- function(ranger_obj, n_trees){
  summaries <- list()
  for(tree_n in 1:n_trees){
    tree1 <- treeInfo(ranger_obj, tree=tree_n)
    tree1 <- calculate_tree_depth_ranger(tree1)
    Dt <- max(tree1[,10])
    tree1$leftChild <- tree1$leftChild + 1
    tree1$rightChild <- tree1$rightChild + 1
    splitvars <- tree1$splitvarName
    splitvars <- splitvars[!is.na(splitvars)]
    split_summary <- list()
    for (splitvar in unique(splitvars)){
      tree <- tree1 %>% filter(splitvarName==splitvar)
      subtrees <- list()
      for (i in tree$nodeID){
        A <- i + 1
        child_parent <- matrix(ncol=4, nrow=nrow(tree1)*nrow(tree1))
        for (j in A:nrow(tree1)){
          parent = tree1[j,5]
          child1_n = tree1[j,2]
          child2_n= tree1[j,3]
          child1 = tree1[child1_n,5]
          child2 = tree1[child2_n,5]
          depth_initial = tree1[j,10]
          depth1 = tree1[child1_n, 10] 
          depth2 = tree1[child2_n, 10] 
          row_1 = 2*j-1
          row_2 = 2*j
          child_parent[row_1,1] = parent
          child_parent[row_1,2] = child1
          child_parent[row_1,3] = depth_initial
          child_parent[row_1,4] = depth1
          child_parent[row_2,1] = parent
          child_parent[row_2,2] = child2
          child_parent[row_2,3] = depth_initial
          child_parent[row_2,4] = depth2
          }
        child_parent<- data.frame(child_parent, stringsAsFactors = FALSE)
        colnames(child_parent) <- c("parent","child","depth_A", "depth_B")
        min_depth <- min(as.numeric(as.character(child_parent$depth_A)), na.rm=TRUE)
        child_parent$depth_AB <- as.numeric(as.character(child_parent$depth_B)) - min_depth
        child_parent$depth_A <- min_depth
        child_parent$parent <- splitvar
        child_parent<- child_parent[!is.na(child_parent$child),]
        subtrees[[A]] <- child_parent        
        }
      subtrees <- rbindlist(subtrees)
      subtrees <- data.frame(subtrees)
      subtrees <- subtrees[!is.na(subtrees),]
      subtrees_final <- subtrees %>% group_by(parent, child) %>% slice(which.min(depth_A))
      subtrees_final$OAB <- 1
      subtrees_final$Dt <- Dt
      subtrees_final$Tree <- tree_n
      split_summary[[splitvar]] <- subtrees_final
      }
    split_summary_df <- rbindlist(split_summary)
    summaries[[tree_n]] <- split_summary_df
    }
    return(rbindlist(summaries))
  }
#Calculate interaction coefficients
get_coef <- function(interactions_frame){ 
  OB <- data.frame(interactions_frame[,2],interactions_frame[,8])
  OB <- unique(OB)
  OB$num <- 1
  OB <- OB %>% group_by(child) %>% summarise(sum_OB = sum(num))
  OB <- data.frame(OB)
  interactions_frame <- data.frame(interactions_frame)
  interactions_frame$sum_cmd <- (interactions_frame$Dt - interactions_frame$depth_AB- 1)/interactions_frame$Dt
  interactions_frame2 <- interactions_frame %>% group_by(parent, child) %>% summarise(sum_DAB = sum(sum_cmd), sum_OAB = sum(OAB))
   interactions_frame <- left_join(interactions_frame2, OB)
  interactions_frame$coeff <- (interactions_frame$sum_OAB/interactions_frame$sum_OB)*(interactions_frame$sum_DAB/interactions_frame$sum_OAB)
  interactions_frame <- interactions_frame %>%
    rowwise() %>%      # for each row
    mutate(interactions = paste(sort(c(parent, child)), collapse = " - ")) %>%  # sort the teams alphabetically and then combine them separating with -
    ungroup() 
  interactions_frame <- aggregate(
    cbind(sum_DAB, sum_OAB, sum_OB, coeff) ~  interactions, 
    data = interactions_frame,
    FUN = mean)   
  interactions_frame <- interactions_frame %>% separate(interactions, c("Mutation_1", "Mutation_2"), sep="-")
  return(interactions_frame)
  }
model <- readRDS("model",opt$drug,".rds"))
interactions <- get_interactions(model[[opt$drug]]$finalModel, 1000)
interactions2 <- get_coef(interactions)
write.csv(rif_interactions2, paste0(opt$drug,"_interactions2.csv"))


