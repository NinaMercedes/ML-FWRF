#Load libraries
library(MLeval)
library(caret)
library(randomForestExplainer)
library(ranger)
library(dplyr)
set.seed(10)
option_list <- list(
  make_option(c("-x", "--train"), type="character", default=FALSE,
    help="csv file for train dataset"),
  make_option(c("-y", "--test"), type="character", default=FALSE,
    help="csv file for test dataset"),
  make_option(c("-p", "--phenotypes"), type="character", default=TRUE,
    help="phenotypes"),
  make_option(c("-l", "--pheno_order"), type="character", default=FALSE,
    help="order of labels to use in chain"),
  make_option(c("-w", "--fw_file"), type="character", default=TRUE,
    help="file with feature weights"), 
  make_option(c("-o", "--output_id"), type="character", default=TRUE,
    help="id to add to file name eg order1")
    )
#parse options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)
#Read in files
train_dataset <- read.csv(opt$train, header=TRUE)
rownames(train_dataset) <- train_dataset$id
train_dataset <- train_dataset[,-1]
fw_file <- feature_weights_file
test_dataset <- read.csv(opt$test, header=TRUE)
rownames(test_dataset) <- test_dataset$id
test_dataset <- test_dataset[,-1]
#Functions for normalisation
norm <- function(x) {
  return ((x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)))
}

norm_fitch <- function(fitch_file, train_geno){
  fitch_df <- read.csv(fitch_file, header=TRUE)
  fitch_df <- fitch_df %>% filter(Mutation %in% colnames(train_geno[,-ncol(train_geno)]))
  fitch_df$norm <- norm(fitch_df[,2])
  return(fitch_df$norm)
}
#Read in weights
sw <- norm_fitch(opt$fw_file, train_data)
#tunegrid
tunegrid <- expand.grid(
  mtry = c(sqrt(10000)), #default is to use root
  min.node.size = c(1), # default for classification
  splitrule = c("extratrees","gini"))
#order of phenotypes
order_label = opt$pheno_order
phenotypes = opt$phenotypes
#model function
train_model <- function(train_data, phenotypes, order_label, select_weights, tunegrid, cv_number=5, num_trees=1000, max_depth=10, seed=10, metric="ROC", importance="impurity"){
  dataset = train_data
  dataset2 = train_data[,!colnames(train_data) %in% phenotypes]
  models <- list()
  prediction <- list()
  for (label in order_label){
        labelled_data <- dataset[, label]
        dataset2$phenotype<-  labelled_data
        dataset_complete <- dataset2[!is.na(dataset2$phenotype),]
        model_weights <- ifelse(dataset_complete$phenotype == "0",(1/table(dataset_complete$phenotype)[1]) * 0.5,(1/table(dataset_complete$phenotype)[2]) * 0.5)
        dataset_complete$phenotype <- factor(dataset_complete$phenotype, levels = c("0","1"), labels = c("S", "R"))
        if (label == order_label[1]){
          mod = train(phenotype~., 
            data = dataset_complete, 
            method = "ranger", 
            metric= metric, 
            trControl = trainControl(method="cv",  number = cv_number, allowParallel = TRUE, verbose = TRUE ,savePredictions = TRUE, classProbs = TRUE, summaryFunction=twoClassSummary),
            tuneGrid = tunegrid,
            num.trees = num_trees, 
            weights=model_weights, 
            importance = importance, 
            max.depth= max_depth, 
            seed= seed, 
            split.select.weights=select_weights)
          prediction[[label]] = predict(mod, dataset2)
          prediction[[label]] = as.character(prediction[[label]])
          prediction[[label]] <- ifelse(prediction[[label]]=="S",0,1)
         dataset2 <- cbind(dataset2, prediction[[label]])  
         colnames(dataset2)[ncol(dataset2)]<- label
          }  else {
            select_weights <- c(select_weights, 1)
            mod= train(phenotype~., 
              data = dataset_complete, 
              method = "ranger", 
              metric= metric,
              trControl = trainControl(method="cv",  number = cv_number, allowParallel = TRUE, verbose = TRUE ,savePredictions = TRUE, classProbs = TRUE, summaryFunction=twoClassSummary), 
              tuneGrid = tunegrid, 
              num.trees = num_trees, 
              weights=model_weights, 
              importance = importance, 
              max.depth= max_depth, 
              seed=seed,
              split.select.weights=select_weights)
            prediction[[label]] = predict(mod, dataset2)
            prediction[[label]] = as.character(prediction[[label]])
            prediction[[label]] <- ifelse(prediction[[label]]=="S",0,1)
            dataset2 <- cbind(dataset2, prediction[[label]])
            colnames(dataset2)[ncol(dataset2)]<- label
          }
        models[[label]] <- mod
      }
      return(models)
    }
#Function to predict    
model.predict <- function(models, test_data, phenotypes, order_label){
  y_dataset = test_data
  y_dataset2 = test_data[,!colnames(test_data) %in% phenotypes]
  y_prediction <- list()
  for (label in order_label){
        labelled_data2 <- y_dataset[, label]
        y_dataset2$phenotype<-  labelled_data2
        y_dataset2$phenotype <- factor(y_dataset2$phenotype, levels = c("0","1"), labels = c("S", "R"))
        if (label == order_label[1]){
          y_prediction[[label]] = predict(models[[label]], y_dataset2)
          y_prediction[[label]] = as.character(y_prediction[[label]])
          y_prediction[[label]] <- ifelse(y_prediction[[label]]=="S",0,1)
         y_dataset2 <- cbind(y_dataset2, y_prediction[[label]])  
         colnames(y_dataset2)[ncol(y_dataset2)]<- label
          }  else {
            y_prediction[[label]] = predict(models[[label]], y_dataset2)
            y_prediction[[label]] = as.character(y_prediction[[label]])
            y_prediction[[label]] <- ifelse(y_prediction[[label]]=="S",0,1)
            y_dataset2 <- cbind(y_dataset2, y_prediction[[label]])
            colnames(y_dataset2)[ncol(y_dataset2)]<- label
          }
      }
      return(y_prediction)
  }
#Function to predict probability
model.predict.proba <- function(models, test_data, phenotypes, order_label){
  y_dataset = test_data
  y_dataset2 = test_data[,!colnames(test_data) %in% phenotypes]
  y_prediction <- list()
  for (label in order_label){
        labelled_data2 <- y_dataset[, label]
        y_dataset2$phenotype<-  labelled_data2
        if (label == order_label[1]){
          y_prediction[[label]] = predict(models[[label]], y_dataset2, type="prob")
          y = predict(models[[label]], y_dataset2)
          y = as.character(y)
          y <- ifelse(y=="S",0,1)
         y_dataset2 <- cbind(y_dataset2, y)  
         colnames(y_dataset2)[ncol(y_dataset2)]<- label
          }  else {
            y_prediction[[label]] = predict(models[[label]], y_dataset2, type="prob")
            y = predict(models[[label]], y_dataset2)
            y = as.character(y)
            y <- ifelse(y=="S",0,1)
            y_dataset2 <- cbind(y_dataset2, y)
            colnames(y_dataset2)[ncol(y_dataset2)]<- label
          }
      }
      return(y_prediction)
  }
#perform training and prediction
MLFWRF <- model.train(train_data=train_dataset, phenotypes=phenotypes, order_label=order_label, select_weights=sw, tunegrid=tunegrid)
pred <- model.predict(models=MLFWRF, test_data= test_dataset, phenotypes=phenotypes, order_label=order_label)
pred2 <- model.predict.proba(models=MLFWRF, test_data= test_dataset, phenotypes=phenotypes, order_label=order_label)
#Function to return results
get_results <- function(test_data, phenotypes, predictions, predictions_proba, model){
  results_1 <- list()
  for (label in phenotypes){
    pred_frame <- data.frame(test_data[,label], predictions[[label]])
    pred_frame2 <- data.frame(pred2[[label]], test_data[,label] )
    colnames(pred_frame) <- c("True_Label", "Pred_Label")
    pred_frame <- pred_frame[!is.na(pred_frame$True_Label),]
    pred_frame2 <- pred_frame2[!is.na(pred_frame2[,3]),]
    pred_frame2[,3] <- factor(pred_frame2[,3], levels = c("0","1"), labels = c("S", "R"))
    pred_frame$Pred_Label<- factor(pred_frame$Pred_Label, levels = c("0","1"), labels = c("S", "R"))
    pred_frame$True_Label<- factor(pred_frame$True_Label, levels = c("0","1"), labels = c("S", "R"))
    conf_matrix <- confusionMatrix(pred_frame$Pred_Label, pred_frame$True_Label, positive="R")
    evaluate <- evalm(model[[label]], showplots=FALSE)
    evaluate <- evaluate$stdres
    evaluate2 <- evalm(pred_frame2, showplots=FALSE)
    evaluate2 <- evaluate2$stdres
    results_1[[label]] <- conf_matrix
    results_1[[label]]$eval <- evaluate
    results_1[[label]]$eval2 <- evaluate2
    }
  return(results_1)  
}
#Get results
results1 <- get_results(test_dataset, phenotypes, pred, pred2, MLFWRF)
accuracy_results <- matrix(nrow=7,ncol=length(phenotypes))
colnames(accuracy_results) <- phenotypes
sensitivity_results <- matrix(nrow=11,ncol=length(phenotypes))
colnames(sensitivity_results) <- phenotypes
crossval_results_AUC <- matrix(nrow=15,ncol=length(phenotypes))
colnames(crossval_results_AUC) <- phenotypes
crossval_results_CI <- matrix(nrow=15,ncol=length(phenotypes))
colnames(crossval_results_CI) <- phenotypes
holdout_results_AUC <- matrix(nrow=15,ncol=length(phenotypes))
colnames(holdout_results_AUC) <- phenotypes
holdout_results_CI <- matrix(nrow=15,ncol=length(phenotypes))
colnames(holdout_results_CI) <- phenotypes
rownames(accuracy_results) <- names(results1[[1]][[3]])
rownames(sensitivity_results) <- names(results1[[1]][[4]])
rownames(crossval_results_AUC) <- rownames(results1[[1]][[7]][[1]])
rownames(crossval_results_CI) <- rownames(results1[[1]][[7]][[1]])
rownames(holdout_results_AUC) <- rownames(results1[[1]][[8]][[1]])
rownames(holdout_results_CI) <- rownames(results1[[1]][[8]][[1]])
for (label in phenotypes){
  accuracy_results[,label] <- results1[[label]][[3]]
  sensitivity_results[,label] <- results1[[label]][[4]]
  crossval_results_AUC[,label] <- results1[[label]][7][[1]][[1]][[1]]
  crossval_results_CI[,label] <- results1[[label]][7][[1]][[1]][[2]]
  holdout_results_AUC[,label] <- results1[[label]][8][[1]][[1]][[1]]
  holdout_results_CI[,label] <- results1[[label]][8][[1]][[1]][[2]]
}
#Write results to csv
write.csv(accuracy_results, paste0("RF_multi/multilabel_accuracy_",opt$output_id,".csv"), row.names=TRUE)
write.csv(sensitivity_results, paste0("RF_multi/multilabel_sensitivity_",opt$output_id,".csv"), row.names=TRUE)
write.csv(crossval_results_AUC, paste0("RF_multi/multilabel__crossval_AUC_",opt$output_id,".csv"), row.names=TRUE)
write.csv(crossval_results_CI, paste0("RF_multi/multilabel_crossval_CI_",opt$output_id,".csv"), row.names=TRUE)
write.csv(holdout_results_AUC, paste0("RF_multi/multilabel__holdout_AUC_",opt$output_id,".csv"), row.names=TRUE)
write.csv(holdout_results_CI, paste0("RF_multi/multilabel_weighted_holdout_CI_",opt$output_id,".csv"), row.names=TRUE)
#Get feature importances
for (label in phenotypes){
 imp <- varImp(mod[[label]])
 write.csv(imp$importance, paste0("RF_multi/multilabel_weighted_importance_",label,".csv"), row.names=TRUE)
}
#Get feature interactions data frame
for (label in phenotypes){
  interactions_frame <- min_depth_interactions(MLFWRF[[label]]$finalModel)
  write.csv(interactions_frame, paste0("RF_multi/multilabel_weighted_interactions",label,".csv"), row.names=TRUE)
}
#Save models
for (label in phenotypes){
  saveRDS(MLFWRF[[label]], paste0("model",label,".rds"))
  }

 

