# Requirements
if (!require(discretization)){
  install.packages("discretization")
  require(discretization)
} else{
  require(discretization)
}

# I assume that we run the script from the scripts directory
setwd("/home/fcalvet/Documents/replicate_yale_analysis/scripts")
setwd("..")
home <- getwd()

ref_home <- paste(home, "/data/reference_sets/", sep = "")
input_home <- paste(home, "/data/input/", sep = "")
day_date = Sys.Date()
# day_date = "2020-05-08"


analysis_home <- paste(home,day_date, sep = "/analysis_")
system(paste("mkdir -p", analysis_home))

analysis_home_data <- paste(analysis_home,"data", sep = "/")
system(paste("mkdir -p", analysis_home_data))

  
# read_data

uORFs_candidates <- paste(input_home, "allENSG_15-5-2020_at_17-12_complete_uORFs_with_exp_homATG.tsv", sep = "")
# uORFs_candidates <- paste(input_home, "complete_uORFs_with_exp_homATG_phyloCSF.tsv", sep = "")
# data <- read.delim(uORFs_candidates)



check_reference <- function(x){
  # this function checks if a value is in the set of reference uORFs
  return (x %in% positive_location)
}



a <- c("heart_science","tools_heart","tools_science")
b <- c("heart_science_intersection","tools_heart_intersection","tools_science_intersection")
c <- c(a, b)

for (name_ref in b){
#for (name_ref in c("tools_science")){
  data <- read.delim(uORFs_candidates)
  
  # choose the highest phyloCSF
  data$phylocsf_real <- data$phylocsf_3
  
  data$mean_expression <- rowMeans(data[40:70])
  
  # name_ref <-"heart_science"
  # Read the data from the known positive uORFs
  positive_cases <- paste(ref_home, name_ref, sep = "")
  positive <- read.delim(positive_cases, header = F)
  # Paste the columns to make the intersection with location column
  positive_location <- paste(paste(paste(positive$V1, positive$V2, sep = ":"), positive$V3,
                                   sep = "-"), positive$V4, sep = ":")
  
  
  data$label <- sapply(data$location, check_reference)
  data$label <- as.factor(as.integer(data$label))
  
  
  cleared_dataframe <- data[,c(9,13,18:32,39:74)]
  write.table(x = cleared_dataframe, file = paste(analysis_home, "/", name_ref, "_uORFs_pre_discretize.tsv", sep = ""),
              sep = "\t", row.names = F)
  
  
  # Combine the positive dataset and unlabeled datasets and discretize them
  data <- read.delim(paste(analysis_home, "/", name_ref, "_uORFs_pre_discretize.tsv", sep = ""))
  IDs <- data$uORF_id
  data$uORF_id <- c()
  data$exon_initiating <- c()
  data$stop_initiating <- c()
  colnames(data)[dim(data)[2]]  <- "class.label"
  
  print("Reached Discretization step")
  # install.packages("discretization")
  library(discretization)
  
  
  # The mdlp discretization algorithm uses the last column of the matrix as the class variables
  # pos.d = cbind(positive_data, rep(1, dim(positive_data)[1]))
  # unl.d = cbind(unlabeled_data, rep(0, dim(unlabeled_data)[1]))
  # neu.d = cbind(neutral_data, rep(2, dim(neutral_data)[1]))
  # colnames(pos.d)[dim(pos.d)[2]] <- "class.label"
  # colnames(unl.d)[dim(unl.d)[2]] <- "class.label"
  # colnames(neu.d)[dim(unl.d)[2]] <- "class.label"
  # all.data = rbind(pos.d, unl.d, neu.d)
  
  
  
  discretized_data <- data.frame(matrix(nrow=nrow(data), ncol=0))
  max_cols <- ncol(data)
  
  # # These are the columns with NA values
  # colnames(data)[colSums(is.na(data)) > 0]
  # # I am replacing NA with -1 in all the columns related with expression
  # # these were the only columns with NA values, and -1 is not a possible
  # # value obtained by experimental means
  data[is.na(data)] = -1
  
  print(paste0("Running mdlp algorithm: ", Sys.time()))
  
  
  for (i in 1:10){
    initial <- i*5-4
    if (initial+4 >= max_cols){
      print(colnames(data)[initial:(max_cols-1)])
      reduced_data <- data[,c(initial:(max_cols-1),max_cols)]    
    } else {
      print(colnames(data)[initial:(initial+4)])
      reduced_data <- data[,c(initial:(initial+4),max_cols)]
    }
  #  print(head(reduced_data))
    i.disc = mdlp(reduced_data)$Disc.data
    discretized_data <- cbind(discretized_data, i.disc)
    print("DONE!")
  }
  
  print(paste0("Discretized the data: ", Sys.time()))
  
  # save.image("./workspace_discretized.RData")
  # write.table(x = discretized_data, file = paste(home, "discretized_data.tsv", sep = ""),
  #             sep = "\t", row.names = F)
  
  c.lab = discretized_data["class.label"]
  
  labels_cols <- which(names(discretized_data) == "class.label")
  all.disc.corrected = subset(discretized_data, select = -c(labels_cols))
  all.disc.corrected = cbind(all.disc.corrected, c.lab)
  colnames(all.disc.corrected)[dim(all.disc.corrected)[2]] = "class.label"
  
  print(paste0("Restored label variable: ", Sys.time()))
  
  write.table(all.disc.corrected, file=paste(analysis_home_data, "/", name_ref, "_discretized_uORF_table_no_row_names.tsv", sep = ""),
              quote=FALSE, sep="\t", col.names=TRUE, row.names = F)
  
  print(paste0("Wrote the discretized data in a file: ", Sys.time()))

}




# # Partitioned protocol
# 
# 
# all_columns <- read.table("discretized_uORF_table_all_2", header=TRUE ,sep="\t",row.names=1,comment.char="", quote="")
# all_columns <- all_columns[,c(1), drop=FALSE]
# 
# for (i in 3:89) {
#   column_current <- read.table(paste("discretized_uORF_table_all_", i, sep=""), header=TRUE ,sep="\t",row.names=1,comment.char="", quote="")
#   column_current <- column_current[,c(1), drop=FALSE]
#   all_columns <- cbind(all_columns, column_current)
#   print(paste("Binding column #", i, "...", sep=""))
# }
# 
# final_column <- read.table(paste("discretized_uORF_table_all_", i, sep=""), header=TRUE ,sep="\t",row.names=1,comment.char="", quote="")
# final_column <- final_column[,c(2), drop=FALSE]
# 
# print("Completing input reading...")
# 
# disc.data <- cbind(all_columns, final_column)
# 
# write.table(disc.data, file=paste("./discretized_data_all_1.txt", sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
