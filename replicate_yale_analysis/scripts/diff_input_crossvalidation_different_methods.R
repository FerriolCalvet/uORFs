# Requirements
if (!require(discretization)){
  install.packages("discretization")
  require(discretization)
} else{
  require(discretization)
}

if (!require(ROCR)){
  install.packages("ROCR")
  require(ROCR)
} else{
  require(ROCR)
}

##################################################
########## DEFINE PATHS FOR THE ANALYSIS #########
##################################################

setwd("/home/fcalvet/Documents/replicate_yale_analysis/scripts")
setwd("..")
home <- getwd()

# day_date = Sys.Date()
day_date = "2020-05-25"
# day_date = "2020-05-18"
# day_date = "2020-05-08"

analysis_home <- paste(home,day_date, sep = "/analysis_")
system(paste("mkdir -p", analysis_home))


ref_home <- paste(home, "/data/reference_sets/", sep = "")
input_home <- paste(home, "/data/input/", sep = "")

analysis_home_data <- paste(analysis_home,"data", sep = "/")
system(paste("mkdir -p", analysis_home_data))


analysis_home_plots <- paste(analysis_home,"plots", sep = "/")
analysis_home_plots_histograms <- paste(analysis_home_plots,"var_histogram", sep = "/")
system(paste("mkdir -p", analysis_home_plots))
system(paste("mkdir -p", analysis_home_plots_histograms))


analysis_home_out_files <- paste(analysis_home,"output", sep = "/")
system(paste("mkdir -p", analysis_home_out_files))

analysis_home_prob_matrices <- paste(analysis_home_out_files, "prob_matrices", sep = "/")
system(paste("mkdir -p", analysis_home_prob_matrices))

analysis_home_out_files_test <- paste(analysis_home_out_files, "test", sep = "/")
system(paste("mkdir -p", analysis_home_out_files_test))


setwd(analysis_home)




##################################################
############ READ AND MODIFY THE DATA ############    
##################################################

name_ref <- "heart_science_intersection"
disc.data <- read.delim(paste(analysis_home_data, "/", name_ref, "_discretized_uORF_table_no_row_names.tsv", sep = ""))
# 1 length
# 7 riboseq coverage
# 6 sequence conservation
# 1 homology
# 31 tissue expression
# 1 tissue entropy
# 1 phylocsf real
# 1 mean expression
# 1 class.label
names(disc.data)[c(16:46)]

mod.disc.data <- disc.data[,-c(16:46)]
mod.disc.data <- mod.disc.data[, c(1:14, rep(15, 5), rep(16, 5), rep(17, 5), rep(18, 5), 19)]
# View(mod.disc.data)
# names(mod.disc.data)
# 1 length
# 7 riboseq coverage
# 6 sequence conservation
# (1 homology) x5
# (1 tissue entropy) x5
# (1 phylocsf real) x5
# (1 mean expression) x5
# 1 class.label
no_label_data <- mod.disc.data[,-ncol(mod.disc.data)]





#the dimensions of the discretized data are measured

d.M = dim(mod.disc.data)[1]
d.N = dim(mod.disc.data)[2]

#the positive and unlabelled data sets are separated based on the class column, the class column is removed.

all_pos.disc = (mod.disc.data[mod.disc.data[,d.N] == 1,])[, -d.N]	# shaves off the class column
all_unl.disc = (mod.disc.data[mod.disc.data[,d.N] == 0,])[, -d.N]	# shaves off the class column



#retaining a number of the examples from the positive set, to test the quality of the algorithm:

d.p.M = dim(all_pos.disc)[1] #number of rows in the positive example.
d.p.N = dim(all_pos.disc)[2] #number of columns in the positive example.
tenpercentsamplesize = trunc(d.p.M/10)
tenpercentremainder = d.p.M %% 10

# In the following section of code, looping over the data set is accomplished, to make multiple retained sets,
# and multiple products.

#First we randomly sample the positive data:

random_max = length(all_pos.disc[,1])
nonrandom_vector = seq(1,random_max,1)
random_vector = sample(nonrandom_vector)
all_pos.disc.mixed = all_pos.disc[random_vector,] # we shuffle the positive data

retained = list()
for (i in 0:9) {
  retained[[(i+1)]] = all_pos.disc.mixed[(i*tenpercentsamplesize+1):((i+1)*tenpercentsamplesize),1:d.p.N]
}

if (tenpercentremainder == 0) {
  positivestail = list()
  positivestail[[10]] = list()
  for (i in 1:9) {
    positivestail[[(i)]] = all_pos.disc.mixed[(i*tenpercentsamplesize+1):d.p.M,1:d.p.N]
  }
} else {
  positivestail = list()
  for (i in 1:10) {
    positivestail[[(i)]] = all_pos.disc.mixed[(i*tenpercentsamplesize+1):d.p.M,1:d.p.N]
  }
}

positiveshead = list()
positiveshead[[(1)]] = list()
for (i in 2:10) {
  positiveshead[[(i)]] = all_pos.disc.mixed[1:((i-1)*tenpercentsamplesize),1:d.p.N]
}

positives = list()
for (i in 1:10) {
  positives[[(i)]] = rbind(positiveshead[[i]],positivestail[[i]])
}


#now add the positive test group, back to the unlabelled set, for later retrieval
bound = list()
for (i in 1:10) {
  bound[[i]] = rbind(retained[[i]],all_unl.disc)
}

#-----------------------------------------------
# The following will be looped:

#The dimensions of both the unlabelled, and positive, matrices of discretized values are measured.
setwd(analysis_home_out_files_test)

#The maximum number of values for a variable is measured.
max.len = 1
for (k in 1:d.N) {
  len = length(table(mod.disc.data[,k]))
  if (len > max.len) max.len <- len
}
print(max.len)














#############################################

#-------------------------------------
# Function Definitions

# This function allows us to remove Inf values from a vector
#   so that we can compute the mean of the values
check_inf <- function(x){
  max_not_inf <- max(x[x != Inf])
  x[x == Inf] = max_not_inf
  
  min_not_inf <- min(x[x != -Inf])
  x[x == -Inf] = min_not_inf
  
  return(x)
}


# this function is for checking that all columns have some variation
screen = function(M) {
  m = dim(M)[1]
  n = dim(M)[2]
  bad.c = c()
  for (k in 1:n) {
    fr = mean(M[,k] == 1)
    if (fr == 0.0 || fr == 1.00) {bad.c = c(bad.c, k)}
  }
  if (length(bad.c) > 0){
    return(M[, -c(bad.c)])  
  } else {
    print("No variable was homogeneous")
    return(M)
  }
  
}


percent_pos = function(training, pred, class, freq=F) {
  # If both are 1 in the same place, then we have percent of positives
  n_class = sum(training == class)
  n_pos = sum(((training == class) + (pred == 1)) == 2)
  if (freq) return(n_pos / n_class)
  else return(n_pos)
}



get.perf.measures = function(pred, labels) {
  p_pos_pos = percent_pos(labels, pred, T, freq = T)
  p_unl_pos = percent_pos(labels, pred, F, freq = T)
  
  recall = p_pos_pos			# These are the same thing
  
  # Precision calculation
  total_length = length(pred)
  pos_pos = percent_pos(labels, pred, T)
  unl_pos = percent_pos(labels, pred, F)
  # prec = pos_pos / (pos_pos + unl_pos)
  # F.stat = 2 * (prec*recall) / (prec + recall)
  prob_positive = (pos_pos + unl_pos) / total_length
  
  F.stat_estimate = (recall * recall) / prob_positive
  
  print(F.stat_estimate)
  
  return(list(p_pos_pos, p_unl_pos, prob_positive, recall, F.stat_estimate))
  
  # recall = true-positives / (true-positives + false-negatives)
  # Precision = true-positives / (true-positives + false-positives)
  #     Precision cannot be computed in PU learning due to the unknown number of real positives
  #     we will use the probability of being predicted as positive to compute an alternative
  #     score with the same behaviour.
  
  # F statistic is the geometric mean of Precision and Recall:
  # F = 2 * (precision * recall) / (precision + recall)
  #     however as we cannot compute precision properly we are using an alternative score with the same behaviour
  # F.alt = (recall * recall) / P(y_prediction=1)
}

n.Bayes = function(all.data, u.prb, p.prb, prior) {
  M = dim(all.data)[1]
  N = dim(all.data)[2]
  data  = all.data[,-N]
  labels = all.data[,N]
  pred = logical(M)
  
  for (i in 1:M) {
    pred[i] = pred.Bayes(as.numeric(data[i,]), u.prb, p.prb, prior)
  }
  
  return(get.perf.measures(pred, labels))
}


# p is the probability of positive (prior probability)
pred.Bayes = function(x, u.prb, p.prb, p, ratio=FALSE) {
  
  p.unl = log10((1-p)) # we define the probability of being unlabelled
  n = dim(u.prb)[1] # n is the total number of variables, and is retrieved as the number of rows in the probability matrix
  
  # we iterate for each of the variables and compute the probabilities of being positive or unlabelled
  for (j in 1:n) {
    # we compute the probability of being unlabelled by retrieving the value from the matrix
    # and we add it to the total probability to add up all the features using the logaritmic value
    p.unl = p.unl + log10(u.prb[j, x[j]])
  }
  # print(p.prb)
  # print(u.prb)
  
  # we initialize the positive probability with the logarithm of the prior probability
  p.pos = log10(p)
  for (j in 1:n) {
    # we compute the probability of being unlabelled by retrieving the value from the matrix
    # and we add it to the total probability to add up all the features using the logaritmic value
    p.pos = p.pos + log10(p.prb[j, x[j]])
  }
  # print(p.pos)
  # print(p.unl)
  if (ratio == TRUE) {
    return(p.pos - p.unl)
  }
  if (p.pos > p.unl) return(TRUE)
  else return(FALSE)
}



# this function is the first function called by the write predictions function
# recevies as input the same as the write.predictions function except the filename
n.Bayes.score = function(data, u.prb, p.prb, prior) {
  scores = data[,1]	# To retain the row names, identifiers for the samples(/uORFs)
  M = dim(data)[1]  # get the number of rows as we will iterate once per each

  for (i in 1:M) {
    # for each row/sample, we pass the data of that sample to the function pred.Bayes
    # as well as the probability matrices of the positive and unlaballed samples, and the prior
    scores[i] = pred.Bayes(as.numeric(data[i,]), u.prb, p.prb, prior, ratio=T)
  }
  return(scores)
}



sample.mat = function(M, p) {
  n = dim(M)[1]
  m = dim(M)[2]
  f = 100000
  p.data = as.integer(p * f)
  log = c(rep(T, p.data), rep(F, f - p.data))
  log.s = sample(log, n, replace=TRUE)
  return(M[log.s,])
}



write.predictions = function(data, u.prb, p.prb, prior, filename="uORF_function_unlabelled_predictions_sorted_") {
  # receives as input:
  # data: a matrix with the values for each of the variables
  # u.prb: probability matrix of the unlabelled samples
  # p.prb: probability matrix of the positive samples
  # prior: prior probability of a uORF beign positive
  # filename: name of the file to store the output
  
  # compute the scores for each uORF, in fact they are probability ratios
  # but we use them as scores
  prob.ratios = n.Bayes.score(data, u.prb, p.prb, prior)
  
  # store the scores into a dataframe and sort them in the proper order
  uORFs.ratio = data.frame(prob.ratios)
  rownames(uORFs.ratio) <- rownames(data)
  uORFs.ratio[,2] <- prob.ratios				# To keep rownames when I sort the data frame
  
  # sort the uORFs according to the scores, and put them into a dataframe
  b = uORFs.ratio[order(-uORFs.ratio[,1]),]
  sorted.scores = data.frame(b[,1])
  rownames(sorted.scores) <- rownames(b)
  
  # write the predictions into the desired file
  output = paste(filename, sep="")
  write.table(sorted.scores, file=output, quote=FALSE, sep="\t", col.names=FALSE)
  
  # return the filename where we saved the data
  return(output)
}





compute.prob.matrix <- function(compute_prob_unlabelled, compute_prob_positive, max.len){
  # compute_prob_unlabelled = it_unlabelled
  # compute_prob_positive = it_positive
  u.N = dim(compute_prob_unlabelled)[2]
  
  #probability matrices are constructed (empty at this stage)
  u.prb_int = matrix(0, nrow=u.N, ncol=max.len)
  p.prb_int = matrix(0, nrow=u.N, ncol=max.len)
  dimnames(u.prb_int)[1] = dimnames(compute_prob_unlabelled)[2]
  dimnames(p.prb_int)[1] = dimnames(compute_prob_positive)[2]
  
  #numbers are entered into the probability matrices, indicating the probability a given variable will have a given value.
  for (k in 1:u.N) {
    for (j in 1:max.len) {
      u.prb_int[k, j] = mean(compute_prob_unlabelled[,k] == j) # as the output from the equality is either True or False, computing the mean is the probability
      p.prb_int[k, j] = mean(compute_prob_positive[,k] == j)	
    }
  }
  
  return(list(u.prb_int, p.prb_int))
}




compute.single.prob.matrix <- function(compute_prob, max.len){
  # compute_prob_unlabelled = it_unlabelled
  # compute_prob_positive = it_positive
  N = dim(compute_prob)[2]
  
  #probability matrices are constructed (empty at this stage)
  prb_int = matrix(0, nrow=N, ncol=max.len)
  dimnames(prb_int)[1] = dimnames(compute_prob)[2]
  
  #numbers are entered into the probability matrices, indicating the probability a given variable will have a given value.
  for (k in 1:N) {
    for (j in 1:max.len) {
      prb_int[k, j] = mean(compute_prob[,k] == j) # as the output from the equality is either True or False, computing the mean is the probability
    }
  }
  
  return(prb_int)
}






train.matrices.updating.unlabelled = function(data, unlabelled, positive, max.len, prior,
                                              it_min = 1,
                                              it_max = 15,
                                              quantile_threshold = 0.97,
                                              update_positive_score_threshold = 5,
                                              skip_differences = F,
                                              filename=paste(analysis_home_prob_matrices, "probability_matrix_updating_unlabelled_", sep = "")) {
  # receives as input:
  # data: a matrix with the samples in rows and variables in columns
  #       at each cell there is the value for that sample and variable
  #       IT MUST BE THE SAME AS THE UNLABELLED DATASET
  # unlabelled: dataset of unlabelled samples that will be used for computing the probability matrix
  # positive: dataset of positive samples that will be used for computing the probability matrix
  # prior: prior probability of a uORF being positive
  # filename: name of the file to store the output
  
  # data = all_unl.disc
  # unlabelled = all_unl.disc
  # positive = all_pos.disc
  # max.len = max.len
  # we use this dataframe to store the stats
  
  internal_internal_stats_df <- data.frame(matrix(ncol = 8, nrow = 0))
  
  names_unlabelled <- rownames(unlabelled)
  names_unlabelled_chosen_positives <- c()
  it_unlabelled <- unlabelled
  
  
  diff_between_matrices <- 0
  
  
  p.prb = compute.single.prob.matrix(positive, max.len)
  
  #for (i in 1:10){
  for (i in 1:it_max){
    cat(paste("\n\nStarting iteration", i))
    print("")
    u.prb = compute.single.prob.matrix(it_unlabelled, max.len = max.len)
    #prob_matrices = compute.prob.matrix(it_unlabelled, positive, max.len = max.len)
    
    
    
    print("Matrix computed")
    
    # compute the scores for each uORF, in fact they are probability ratios
    # but we use them as scores
    prob.ratios = n.Bayes.score(it_unlabelled, u.prb, p.prb, prior)
    
    print("Scores computed")
    
    # store the scores into a dataframe and sort them in the proper order
    uORFs.ratio = data.frame(prob.ratios)
    rownames(uORFs.ratio) <- rownames(it_unlabelled) # these are the rownames from the initial data matrix
    #                                         where labelled and unlabelled were together
    uORFs.ratio[,2] <- prob.ratios				# To keep rownames when I sort the data frame
    
    # sort the uORFs according to the scores, and put them into a dataframe
    b = uORFs.ratio[order(-uORFs.ratio[,1]),]
    sorted.scores = data.frame(b[,1])
    rownames(sorted.scores) <- rownames(b)
    
    print("uORFs sorted")
    
    positive_scored <- b[b[,1] > 0,]
    print(paste("Number of positives", nrow(positive_scored), nrow(positive_scored) + length(names_unlabelled_chosen_positives)))
    # print(summary(positive_scored))
    
    threshold_positive <- quantile(positive_scored[,1], quantile_threshold) # I should adjust this threshold
    print(paste("Threshold quantile", threshold_positive))
    
    
    selected_from_unlabelled <- rownames(sorted.scores)[sorted.scores[,1] > threshold_positive] # get the positions of the uORFs above the threshold
    
    names_unlabelled_chosen_positives <- c(names_unlabelled_chosen_positives,
                                           selected_from_unlabelled)
    
    print(paste("chosen_positive", length(selected_from_unlabelled), length(names_unlabelled_chosen_positives)))
    names_unlabelled <- setdiff(names_unlabelled, names_unlabelled_chosen_positives)
    
    print(paste("Unlabelled", length(names_unlabelled)))
    
    it_unlabelled <- data[names_unlabelled,]
    print("Unlabelled filtered")
    
    
    new_difference <- sum(abs(u.prb - p.prb))
    print(paste("The new difference is", new_difference))
    
    
    mean_prob_sum_int <- mean(abs(check_inf(b[,1])))
    print(paste("Mean probability of all remaining unlabelled (absolute value)", mean_prob_sum_int))
    
    
    # before finishing each iteration, we store some stats
    internal_internal_stats_df <- rbind(internal_internal_stats_df, c(i,
                                                                      nrow(positive_scored) + length(names_unlabelled_chosen_positives),
                                                                      threshold_positive,
                                                                      length(names_unlabelled),
                                                                      0,
                                                                      0,
                                                                      new_difference,
                                                                      mean_prob_sum_int))
    
    
    
    if (!skip_differences){
      if (new_difference < diff_between_matrices + 0.1){
        break
      } else {
        diff_between_matrices <- new_difference
      }
    }
    
    if (threshold_positive < update_positive_score_threshold ){ break }
    
    
    
  }
  filename_u = paste(filename, "u.prb", sep = "_")
  filename_p = paste(filename, "p.prb", sep = "_")
  write.table(u.prb, file=filename_u, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
  write.table(p.prb, file=filename_p, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)  
  
  # return the two probability matrices
  return(list(u.prb, p.prb, internal_internal_stats_df))
}






train.matrices.updating.unlabelled.positive = function(data, unlabelled, positive, max.len, prior,
                                                       it_min = 1,
                                                       it_max = 15,
                                                       quantile_threshold = 0.97,
                                                       update_positive_score_threshold = 5,
                                                       skip_differences = F,
                                                       filename=paste(analysis_home_prob_matrices,"probability_matrix_updating_unlabelled_positive_")) {
  # receives as input:
  # data: a matrix with the samples in rows and variables in columns
  #       at each cell there is the value for that sample and variable
  #       IT MUST BE THE SAME AS THE UNLABELLED DATASET
  # unlabelled: dataset of unlabelled samples that will be used for computing the probability matrix
  # positive: dataset of positive samples that will be used for computing the probability matrix
  # prior: prior probability of a uORF being positive
  # filename: name of the file to store the output
  
  # data = all_unl.disc
  # unlabelled = all_unl.disc
  # positive = all_pos.disc
  # max.len = max.len
  
  # we use this dataframe to store the stats
  internal_internal_stats_df <- data.frame(matrix(ncol = 8, nrow = 0))
  
  names_unlabelled <- rownames(unlabelled)
  names_unlabelled_chosen_positives <- c()
  it_unlabelled <- unlabelled
  it_positive <- positive
  
  diff_between_matrices <- 0
  
  #for (i in 1:10){
  for (i in 1:it_max){
    cat(paste("\n\nStarting iteration", i))
    print("")
    prob_matrices = compute.prob.matrix(it_unlabelled, it_positive, max.len = max.len)
    #prob_matrices = compute.prob.matrix(it_unlabelled, positive, max.len = max.len)
    
    u.prb = prob_matrices[[1]]
    p.prb = prob_matrices[[2]]
    
    print("Matrices computed")
    
    # compute the scores for each uORF, in fact they are probability ratios
    # but we use them as scores
    prob.ratios = n.Bayes.score(it_unlabelled, u.prb, p.prb, prior)
    
    print("Scores computed")
    
    # store the scores into a dataframe and sort them in the proper order
    uORFs.ratio = data.frame(prob.ratios)
    rownames(uORFs.ratio) <- rownames(it_unlabelled) # these are the rownames from the initial data matrix
    #                                         where labelled and unlabelled were together
    uORFs.ratio[,2] <- prob.ratios				# To keep rownames when I sort the data frame
    
    # sort the uORFs according to the scores, and put them into a dataframe
    b = uORFs.ratio[order(-uORFs.ratio[,1]),]
    sorted.scores = data.frame(b[,1])
    rownames(sorted.scores) <- rownames(b)
    
    # print("uORFs sorted")
    
    positive_scored <- b[b[,1] > 0,]
    print(paste("Number of positives", nrow(positive_scored), nrow(positive_scored) + length(names_unlabelled_chosen_positives)))
    # print(summary(positive_scored))
    
    threshold_positive <- quantile(positive_scored[,1], quantile_threshold) 
    print(paste("Threshold quantile", threshold_positive))
    
    
    selected_from_unlabelled <- rownames(sorted.scores)[sorted.scores[,1] > threshold_positive] # get the positions of the uORFs above the threshold
    
    names_unlabelled_chosen_positives <- c(names_unlabelled_chosen_positives,
                                           selected_from_unlabelled)
    
    print(paste("chosen_positive", length(selected_from_unlabelled), length(names_unlabelled_chosen_positives)))
    
    names_unlabelled <- setdiff(names_unlabelled, names_unlabelled_chosen_positives)
    
    print(paste("Unlabelled", length(names_unlabelled)))
    
    it_unlabelled <- data[names_unlabelled,]
    print("Unlabelled filtered")
    
    rows_from_unlabelled_to_positive <- data[as.integer(selected_from_unlabelled),]
    
    it_positive <- rbind(it_positive, rows_from_unlabelled_to_positive)
    print("Positive updated")
    
    
    new_difference <- sum(abs(u.prb - p.prb))
    print(paste("The new difference is", new_difference))
    
    
    mean_prob_sum_int <- mean(abs(check_inf(b[,1])))
    print(paste("Mean probability of all remaining unlabelled (absolute value)", mean_prob_sum_int))
    
    
    # before finishing each iteration, we store some stats
    internal_internal_stats_df <- rbind(internal_internal_stats_df, c(i,
                                                                      nrow(positive_scored) + length(names_unlabelled_chosen_positives),
                                                                      threshold_positive,
                                                                      length(names_unlabelled),
                                                                      0,
                                                                      0,
                                                                      new_difference,
                                                                      mean_prob_sum_int))
    
    
    
    #if (i > it_min){
      
      if (!skip_differences){
        if (new_difference < diff_between_matrices + 0.1){
          break
        } else {
          diff_between_matrices <- new_difference
        }
      }
      
      if (threshold_positive < update_positive_score_threshold){ break }
    #}
    
  }
  filename_u = paste(filename, "u.prb", sep = "_")
  filename_p = paste(filename, "p.prb", sep = "_")
  write.table(u.prb, file=filename_u, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
  write.table(p.prb, file=filename_p, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)  

  # return the two probability matrices
  return(list(u.prb, p.prb, internal_internal_stats_df))
}



train.matrices.updating.negative = function(data, unlabelled, positive, max.len, prior,
                                            it_min = 1,
                                            it_max = 50,
                                            quantile_threshold_positive = 0.97,
                                            update_positive_score_threshold = 5,
                                            quantile_threshold_negative = 0.97,
                                            update_negative_score_threshold = 5,
                                            skip_differences = F,
                                            filename=paste(analysis_home_prob_matrices,"probability_matrix_updating_negative_")) {
  # receives as input:
  # unlabelled: dataset of unlabelled samples that will be used for computing the probability matrix
  # positive: dataset of positive samples that will be used for computing the probability matrix
  # prior: prior probability of a uORF being positive
  # filename: name of the file to store the output
  
  # data = all_unl.disc
  # unlabelled = all_unl.disc
  # positive = all_pos.disc
  # max.len = max.len
  
  # we use this dataframe to store the stats
  internal_internal_stats_df <- data.frame(matrix(ncol = 8, nrow = 0))
  
  names_unlabelled <- rownames(unlabelled)
  names_unlabelled_chosen_positives <- c()
  names_unlabelled_chosen_negatives <- c()
  it_positive <- positive
  it_negative <- list()
  # it_negative <- unlabelled # starts with all unlabelled but then contains the uORFs that have been selected as negatives
  it_unlabelled <- unlabelled  # it will contain the uORFs that have not been selected as negative nor positive
  
  diff_between_matrices <- 0
  
  p.prb = compute.single.prob.matrix(positive, max.len = max.len)
  
  #for (i in 1:10){
  for (i in 1:it_max){
    cat(paste("\n\nStarting iteration", i))
    print("")
    if (i == 1) {
      n.prb = compute.single.prob.matrix(it_unlabelled, max.len = max.len)
    } else{
      n.prb = compute.single.prob.matrix(it_negative, max.len = max.len)
    }
    #prob_matrices = compute.prob.matrix(it_unlabelled, positive, max.len = max.len)
    #print(head(n.prb))
    print("Matrix computed")
    
    # compute the scores for each uORF, in fact they are probability ratios
    # but we use them as scores
    prob.ratios = n.Bayes.score(it_unlabelled, n.prb, p.prb, prior)
    
    print("Scores computed")
    
    # store the scores into a dataframe and sort them in the proper order
    uORFs.ratio = data.frame(prob.ratios)
    rownames(uORFs.ratio) <- rownames(it_unlabelled) # these are the rownames from the initial data matrix
    #                                         where labelled and unlabelled were together
    uORFs.ratio[,2] <- prob.ratios				# To keep rownames when I sort the data frame
    
    # sort the uORFs according to the scores, and put them into a dataframe
    b = uORFs.ratio[order(-uORFs.ratio[,1]),]
    sorted.scores = data.frame(b[,1])
    rownames(sorted.scores) <- rownames(b)
    
    # print("uORFs sorted")
    
    # Select top positives
    positive_scored <- b[b[,1] > 0,]
    print(paste("Number of positives", nrow(positive_scored), nrow(positive_scored) + length(names_unlabelled_chosen_positives)))
    
    threshold_positive <- quantile(check_inf(positive_scored[,1]), quantile_threshold_positive) 
    print(paste("Threshold positive quantile", threshold_positive))
    
    positive_selected_from_unlabelled <- rownames(sorted.scores)[sorted.scores[,1] > threshold_positive] # get the positions of the uORFs above the threshold
    
    names_unlabelled_chosen_positives <- c(names_unlabelled_chosen_positives,
                                           positive_selected_from_unlabelled)
    
    print(paste("chosen_positive", length(positive_selected_from_unlabelled), length(names_unlabelled_chosen_positives)))
    
    names_unlabelled <- setdiff(names_unlabelled, names_unlabelled_chosen_positives)
    
    
    
    # Select most negatives
    negative_scored <- b[b[,1] < 0,]
    print(paste("Number of negatives", nrow(negative_scored), nrow(negative_scored) + length(names_unlabelled_chosen_negatives)))
    
    # neg_inf_labels <- rownames(negative_scored)[negative_scored[,1] == -Inf]
    # print(length(neg_inf_labels))
    # inf_removed <- negative_scored[-c(as.integer(neg_inf_labels)),]
    # print(nrow(inf_removed))
    
    inf_removed <- negative_scored[negative_scored[,1] != -Inf,]
    
    #if (i == 1 ){
    #  threshold_negative <- quantile(check_inf(inf_removed[,1]), 0.7) 
    #} else {
      threshold_negative <- quantile(check_inf(inf_removed[,1]), 1-quantile_threshold_negative) 
    #}
    
    print(paste("Threshold negative quantile", threshold_negative))
    
    negative_selected_from_unlabelled <- rownames(sorted.scores)[sorted.scores[,1] < threshold_negative] # get the positions of the uORFs above the threshold
    
    names_unlabelled_chosen_negatives <- c(names_unlabelled_chosen_negatives,
                                           negative_selected_from_unlabelled)
    
    print(paste("chosen_negative", length(negative_selected_from_unlabelled), length(names_unlabelled_chosen_negatives)))
    
    names_unlabelled <- setdiff(names_unlabelled, names_unlabelled_chosen_negatives)
    
    
    
    print(paste("Unlabelled", length(names_unlabelled)))
    
    it_unlabelled <- data[names_unlabelled,]
    print("Unlabelled filtered")
    
    rows_from_unlabelled_to_positive <- data[as.integer(positive_selected_from_unlabelled),]
    
    it_positive <- rbind(it_positive, rows_from_unlabelled_to_positive)
    print("Positive updated")
    
    
    rows_from_unlabelled_to_negative <- data[as.integer(negative_selected_from_unlabelled),]
    # print(head(rows_from_unlabelled_to_negative))
    it_negative <- rbind(it_negative, rows_from_unlabelled_to_negative)
    # print(mean(it_negative[,1]))
    print(paste("Negative updated", nrow(it_negative)))
    
    
    new_difference <- sum(abs(n.prb - p.prb))
    print(paste("The new difference is", new_difference))
    
    
    mean_prob_sum_int <- mean(abs(check_inf(b[,1])))
    print(paste("Mean probability of all remaining unlabelled (absolute value)", mean_prob_sum_int))
    
    
    # before finishing each iteration, we store some stats
    internal_internal_stats_df <- rbind(internal_internal_stats_df, c(i,
                                                                      nrow(positive_scored) + length(names_unlabelled_chosen_positives),
                                                                      threshold_positive,
                                                                      length(names_unlabelled),
                                                                      nrow(negative_scored) + length(names_unlabelled_chosen_negatives),
                                                                      threshold_negative,
                                                                      new_difference,
                                                                      mean_prob_sum_int))
    
    
    #if (i > it_min){
      
    if (!skip_differences){
      if (new_difference < diff_between_matrices + 0.1){
        break
      } else {
        diff_between_matrices <- new_difference
      }
    }
      
      
    #if (threshold_positive < update_positive_score_threshold){ break }
    if (threshold_negative > update_negative_score_threshold){ break }
    #}
    
    
    
  }
  # filename_u = paste(filename, "u.prb", sep = "_")
  # filename_p = paste(filename, "p.prb", sep = "_")
  # write.table(u.prb, file=filename_u, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
  # write.table(p.prb, file=filename_p, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)  
  
  # return the two probability matrices
  return(list(n.prb, p.prb, internal_internal_stats_df))
}
# train.matrices.updating.negative(disc.data, bound[[1]], positives[[1]], max.len, prior, it_min = 4, it_max = 50, quantile_threshold_positive = 0.99, update_positive_score_threshold = 5, quantile_threshold_negative = 0.95, update_negative_score_threshold = -10)





train.matrices.updating.negative.positive = function(data, unlabelled, positive, max.len, prior,
                                                     it_min = 1,
                                                     it_max = 15,
                                                     quantile_threshold_positive = 0.97,
                                                     update_positive_score_threshold = 5,
                                                     quantile_threshold_negative = 0.97,
                                                     update_negative_score_threshold = 5,
                                                     skip_differences = F,
                                                     filename=paste(analysis_home_prob_matrices,"probability_matrix_updating_negative_positive_")) {
  # receives as input:
  # unlabelled: dataset of unlabelled samples that will be used for computing the probability matrix
  # positive: dataset of positive samples that will be used for computing the probability matrix
  # prior: prior probability of a uORF being positive
  # filename: name of the file to store the output
  
  # data = all_unl.disc
  # unlabelled = all_unl.disc
  # positive = all_pos.disc
  # max.len = max.len
  
  # we use this dataframe to store the stats
  internal_internal_stats_df <- data.frame(matrix(ncol = 8, nrow = 0))
  
  names_unlabelled <- rownames(unlabelled)
  names_unlabelled_chosen_positives <- c()
  names_unlabelled_chosen_negatives <- c()
  it_positive <- positive
  it_negative <- list()
  # it_negative <- unlabelled # starts with all unlabelled but then contains the uORFs that have been selected as negatives
  it_unlabelled <- unlabelled  # it will contain the uORFs that have not been selected as negative nor positive
  
  diff_between_matrices <- 0
  
  
  
  #for (i in 1:10){
  for (i in 1:it_max){
    cat(paste("\n\nStarting iteration", i))
    print("")
    if (i == 1) {
      n.prb = compute.single.prob.matrix(it_unlabelled, max.len = max.len)
      p.prb = compute.single.prob.matrix(positive, max.len = max.len)
    } else{
      list.prb = compute.prob.matrix(compute_prob_positive = it_positive, compute_prob_unlabelled = it_negative, max.len = max.len)
      n.prb = list.prb[[1]]
      p.prb = list.prb[[2]]
    }
    print("Matrix computed")
    
    # compute the scores for each uORF, in fact they are probability ratios
    # but we use them as scores
    prob.ratios = n.Bayes.score(it_unlabelled, n.prb, p.prb, prior)
    print("Scores computed")
    
    
    # store the scores into a dataframe and sort them in the proper order
    uORFs.ratio = data.frame(prob.ratios)
    rownames(uORFs.ratio) <- rownames(it_unlabelled) # these are the rownames from the initial data matrix
    #                                         where labelled and unlabelled were together
    uORFs.ratio[,2] <- prob.ratios				# To keep rownames when I sort the data frame
    
    # sort the uORFs according to the scores, and put them into a dataframe
    b = uORFs.ratio[order(-uORFs.ratio[,1]),]
    sorted.scores = data.frame(b[,1])
    rownames(sorted.scores) <- rownames(b)
    # print("uORFs sorted")
    
    
    # Select most positives for updating the positive set
    positive_scored <- b[b[,1] > 0,]
    print(paste("Number of positives", nrow(positive_scored), nrow(positive_scored) + length(names_unlabelled_chosen_positives)))
    
    inf_removed_pos <- positive_scored[positive_scored[,1] != Inf,]
    threshold_positive <- quantile(check_inf(inf_removed_pos[,1]), quantile_threshold_positive) 
    print(paste("Threshold positive quantile", threshold_positive))
    
    # get the positions/labels of the uORFs above the threshold
    positive_selected_from_unlabelled <- rownames(sorted.scores)[sorted.scores[,1] > threshold_positive] # get the positions of the uORFs above the threshold
    names_unlabelled_chosen_positives <- c(names_unlabelled_chosen_positives,
                                           positive_selected_from_unlabelled)
    
    print(paste("chosen_positive", length(positive_selected_from_unlabelled), length(names_unlabelled_chosen_positives)))
    
    # update unlabelled names
    names_unlabelled <- setdiff(names_unlabelled, names_unlabelled_chosen_positives)
    
    
    
    # Select most negatives for updating the negative set
    negative_scored <- b[b[,1] < 0,]
    print(paste("Number of negatives", nrow(negative_scored), nrow(negative_scored) + length(names_unlabelled_chosen_negatives)))
    
    inf_removed_neg <- negative_scored[negative_scored[,1] != -Inf,]
    threshold_negative <- quantile(check_inf(inf_removed_neg[,1]), 1-quantile_threshold_negative) 
    print(paste("Threshold negative quantile", threshold_negative))
    
    # get the positions/labels of the uORFs below the threshold
    negative_selected_from_unlabelled <- rownames(sorted.scores)[sorted.scores[,1] < threshold_negative]
    names_unlabelled_chosen_negatives <- c(names_unlabelled_chosen_negatives,
                                           negative_selected_from_unlabelled)
    print(paste("chosen_negative", length(negative_selected_from_unlabelled), length(names_unlabelled_chosen_negatives)))
    
    # update unlabelled names
    names_unlabelled <- setdiff(names_unlabelled, names_unlabelled_chosen_negatives)
    
    # update unlabelled dataset
    it_unlabelled <- data[names_unlabelled,]
    print(paste("Unlabelled updated", nrow(it_unlabelled)))
    
    
    rows_from_unlabelled_to_positive <- data[as.integer(positive_selected_from_unlabelled),]
    it_positive <- rbind(it_positive, rows_from_unlabelled_to_positive)
    print(paste("Positive updated", nrow(it_positive)))
    
    
    rows_from_unlabelled_to_negative <- data[as.integer(negative_selected_from_unlabelled),]
    it_negative <- rbind(it_negative, rows_from_unlabelled_to_negative)
    print(paste("Negative updated", nrow(it_negative)))
    
    
    new_difference <- sum(abs(n.prb - p.prb))
    print(paste("The new difference is", new_difference))
    
    
    mean_prob_sum_int <- mean(abs(check_inf(b[,1])))
    print(paste("Mean probability of all remaining unlabelled (absolute value)", mean_prob_sum_int))
    
    
    # before finishing each iteration, we store some stats
    internal_internal_stats_df <- rbind(internal_internal_stats_df, c(i,
                                                                      nrow(positive_scored) + length(names_unlabelled_chosen_positives),
                                                                      threshold_positive,
                                                                      length(names_unlabelled),
                                                                      nrow(negative_scored) + length(names_unlabelled_chosen_negatives),
                                                                      threshold_negative,
                                                                      new_difference,
                                                                      mean_prob_sum_int))
    
    
    #if (i > it_min){
    
    if (!skip_differences){
      if (new_difference < diff_between_matrices + 0.1){
        break
      } else {
        diff_between_matrices <- new_difference
      }
    }
    
    
    #if (threshold_positive < update_positive_score_threshold){ break }
    if (threshold_negative > update_negative_score_threshold){ break }
    #}
    
    
    
  }
  # filename_u = paste(filename, "u.prb", sep = "_")
  # filename_p = paste(filename, "p.prb", sep = "_")
  # write.table(u.prb, file=filename_u, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
  # write.table(p.prb, file=filename_p, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)  
  
  # return the two probability matrices
  return(list(n.prb, p.prb, internal_internal_stats_df))
}




hist.looping = function(pred, xmin, xmax, breaks=100, main="No main provided", xlab="Positive/unlabelled probability ratio\n prior positive distribution is 0.5", col="lightgreen") {
  hist(pred[pred > xmin & pred < xmax], breaks=breaks, main=main, xlab=xlab, col=col)
}



#The number of variables is measured, this is common for all the methods
d.N = dim(mod.disc.data)[2]
max.len = 1
for (k in 1:d.N) {
  len = length(table(mod.disc.data[,k]))
  if (len > max.len) max.len <- len
}


# store the performance measures at the end of each iteration and model
performance_measures_df <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(performance_measures_df) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                      "recall", "f1_score_like_estimate", "number_of_unlabelled_positive")


stats_internal_iterations_df <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(stats_internal_iterations_df) = c("method","iteration", "internal_iteration", "number_of_positives", "threshold_positive",
                                           "number_of_unlabelled",
                                           "number_of_negatives", "threshold_negative",
                                           "difference_between_matrices",
                                           "mean_absolute_probability")


internal_internal_stats_df <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(stats_internal_iterations_df) = c("internal_iteration", "number_of_positives", "threshold_positive",
                                           "number_of_unlabelled",
                                           "number_of_negatives", "threshold_negative",
                                           "difference_between_matrices",
                                           "mean_absolute_probability")


prior = 0.5

method = 0
cat.unl = list()
cat.pos = list()
cat.ret = list()
for (i in 1:8) {
  
  unl.looping.disc = bound[[i]]
  pos.looping.disc = positives[[i]]
  ret.looping.disc = retained[[i]]
  
  prob_matrices_list0 = compute.prob.matrix(unl.looping.disc, pos.looping.disc, max.len)
  print("Matrices 0 computed!")

  u.prb0 = prob_matrices_list0[[1]]
  p.prb0 = prob_matrices_list0[[2]]

  
  #write.predictions
  unl.pred.file0 = write.predictions(unl.looping.disc, u.prb0, p.prb0, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  pos.pred.file0 = write.predictions(pos.looping.disc, u.prb0, p.prb0, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  ret.pred.file0 = write.predictions(ret.looping.disc, u.prb0, p.prb0, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  
  
  catfile.unl = paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.unl.add = read.table(catfile.unl, header=FALSE,sep="\t",comment.char="", quote="")
  cat.unl = rbind(cat.unl,cat.unl.add)
  catfile.pos = paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.pos.add = read.table(catfile.pos, header=FALSE,sep="\t",comment.char="", quote="")
  cat.pos = rbind(cat.pos,cat.pos.add)
  catfile.ret = paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.ret.add = read.table(catfile.ret, header=FALSE,sep="\t",comment.char="", quote="")
  cat.ret = rbind(cat.ret,cat.ret.add)
  
  
  
  df_for_perf_measures = rbind(cat.unl.add, cat.pos.add, cat.ret.add)
  df_for_perf_measures = cbind(df_for_perf_measures, mod.disc.data$class.label[df_for_perf_measures[,1]])
  df_for_perf_measures[,2] = as.integer(df_for_perf_measures[,2] > 0)
  
  
  
  perf_measures_fragment <- c(method, i, get.perf.measures(df_for_perf_measures[,2], df_for_perf_measures[,3]))
  names(perf_measures_fragment) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                    "recall", "f1_score_like_estimate")
  
  colnames(performance_measures_df) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                        "recall", "f1_score_like_estimate")
  performance_measures_df <- rbind(performance_measures_df,
                                   perf_measures_fragment)
  
  colnames(performance_measures_df) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                        "recall", "f1_score_like_estimate")
  print("performance measures added")
  
  
  print(paste(i, "DONE!"))
  print(cat("\n\n\n\n"))
  
}

svg(paste(analysis_home_plots, '/ret_pos_unl_distro_method_corrected_', method,'.svg', sep = ""))
par(mfrow = c(3,1))
# xlim = c(1.00, 1000)
xlim = c(-100, 100)
hist.looping(cat.pos[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for positive uORFs", col="lightblue")
hist.looping(cat.unl[,2], xlim[1], xlim[2], breaks=40, main="Class probability ratio distribution for unlabelled uORFs", col="lightgreen")
hist.looping(cat.ret[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for retained uORFs")
dev.off()
  
# Some statistics, resulting from the looping algorithm:
quantile(cat.unl[,2], c(.10, .25, .50, .75, .90))
quantile(cat.pos[,2], c(.10, .25, .50, .75, .90))
quantile(cat.ret[,2], c(.10, .25, .50, .75, .90))








method = 1
cat.unl = list()
cat.pos = list()
cat.ret = list()

for (i in 1:8) {
  
  unl.looping.disc = bound[[i]]
  pos.looping.disc = positives[[i]]
  ret.looping.disc = retained[[i]]
  
  quantile_threshold1 = 0.97
  update_positive_score_threshold1 = 10
  prob_matrices_list1 = train.matrices.updating.unlabelled(no_label_data, unl.looping.disc, pos.looping.disc,
                                                           max.len, prior,
                                                           #it_min = 1,
                                                           it_max = 50,
                                                           quantile_threshold = quantile_threshold1,
                                                           update_positive_score_threshold = update_positive_score_threshold1,
                                                           filename=paste(paste(analysis_home_prob_matrices, "/uORF_probability_matrices_updating_unlabelled", sep = ""), i, prior, method, quantile_threshold1, update_positive_score_threshold1, sep='_'))
  print("Matrices 1 computed!")
  
  u.prb1 = prob_matrices_list1[[1]]
  p.prb1 = prob_matrices_list1[[2]]
  iter_stats_fragment = cbind(method, i, prob_matrices_list1[[3]])
  colnames(iter_stats_fragment) <- c("method","iteration","internal_iteration", "number_of_positives", "threshold_positive",
                                     "number_of_unlabelled",
                                     "number_of_negatives", "threshold_negative",
                                     "difference_between_matrices",
                                     "mean_absolute_probability")
  
  colnames(stats_internal_iterations_df) <- c("method","iteration","internal_iteration", "number_of_positives", "threshold_positive",
                                              "number_of_unlabelled",
                                              "number_of_negatives", "threshold_negative",
                                              "difference_between_matrices",
                                              "mean_absolute_probability")
  stats_internal_iterations_df = rbind(stats_internal_iterations_df, iter_stats_fragment)
  colnames(stats_internal_iterations_df) <- c("method","iteration","internal_iteration", "number_of_positives", "threshold_positive",
                                              "number_of_unlabelled",
                                              "number_of_negatives", "threshold_negative",
                                              "difference_between_matrices",
                                              "mean_absolute_probability")
  print("stats added")
  
  #write.predictions
  unl.pred.file1 = write.predictions(unl.looping.disc, u.prb1, p.prb1, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  pos.pred.file1 = write.predictions(pos.looping.disc, u.prb1, p.prb1, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  ret.pred.file1 = write.predictions(ret.looping.disc, u.prb1, p.prb1, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  
  
  catfile.unl = paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.unl.add = read.table(catfile.unl, header=FALSE,sep="\t",comment.char="", quote="")
  cat.unl = rbind(cat.unl,cat.unl.add)
  catfile.pos = paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.pos.add = read.table(catfile.pos, header=FALSE,sep="\t",comment.char="", quote="")
  cat.pos = rbind(cat.pos,cat.pos.add)
  catfile.ret = paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.ret.add = read.table(catfile.ret, header=FALSE,sep="\t",comment.char="", quote="")
  cat.ret = rbind(cat.ret,cat.ret.add)
  
  df_for_perf_measures = rbind(cat.unl.add, cat.pos.add, cat.ret.add)
  df_for_perf_measures = cbind(df_for_perf_measures, mod.disc.data$class.label[df_for_perf_measures[,1]])
  df_for_perf_measures[,2] = as.integer(df_for_perf_measures[,2] > 0)
  
  
  perf_measures_fragment <- c(method, i, get.perf.measures(df_for_perf_measures[,2], df_for_perf_measures[,3]))
  names(perf_measures_fragment) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                    "recall", "f1_score_like_estimate")
  
  colnames(performance_measures_df) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                        "recall", "f1_score_like_estimate")
  performance_measures_df <- rbind(performance_measures_df,
                                   perf_measures_fragment)
  
  colnames(performance_measures_df) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                        "recall", "f1_score_like_estimate")
  
  print("performance measures added")
  
  print(paste(i, "DONE!"))
  print(cat("\n\n\n\n"))
  
}


svg(paste(analysis_home_plots, '/ret_pos_unl_distro_method_corrected_', method,'.svg', sep = ""))
par(mfrow = c(3,1))
# xlim = c(1.00, 1000)
xlim = c(-100, 100)
hist.looping(cat.pos[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for positive uORFs", col="lightblue")
hist.looping(cat.unl[,2], xlim[1], xlim[2], breaks=40, main="Class probability ratio distribution for unlabelled uORFs", col="lightgreen")
hist.looping(cat.ret[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for retained uORFs")
dev.off()

# Some statistics, resulting from the looping algorithm:

quantile(cat.unl[,2], c(.10, .25, .50, .75, .90))
quantile(cat.pos[,2], c(.10, .25, .50, .75, .90))
quantile(cat.ret[,2], c(.10, .25, .50, .75, .90))






method = 2
cat.unl = list()
cat.pos = list()
cat.ret = list()
for (i in 1:8) {
  
  unl.looping.disc = bound[[i]]
  pos.looping.disc = positives[[i]]
  ret.looping.disc = retained[[i]]
  
  
  quantile_threshold2 = 0.97
  update_positive_score_threshold2 = 10
  prob_matrices_list2 = train.matrices.updating.unlabelled.positive(no_label_data,unl.looping.disc, pos.looping.disc,
                                                                    max.len, prior,
                                                                    #it_min = 5,
                                                                    it_max = 50,
                                                                    quantile_threshold = quantile_threshold2,
                                                                    update_positive_score_threshold = update_positive_score_threshold2,
                                                                    #skip_differences = T,
                                                                    filename=paste(paste(analysis_home_prob_matrices, "/uORF_probability_matrices_updating_unlabelled_and_positive", sep = ""), i, prior,  method, quantile_threshold2, update_positive_score_threshold2, sep='_'))
  print("Matrices 2 computed!")
  
  
  u.prb2 = prob_matrices_list2[[1]]
  p.prb2 = prob_matrices_list2[[2]]
  iter_stats_fragment = cbind(method, i, prob_matrices_list2[[3]])
  colnames(iter_stats_fragment) <- c("method","iteration","internal_iteration", "number_of_positives", "threshold_positive",
                                     "number_of_unlabelled",
                                     "number_of_negatives", "threshold_negative",
                                     "difference_between_matrices",
                                     "mean_absolute_probability")
  
  colnames(stats_internal_iterations_df) <- c("method","iteration","internal_iteration", "number_of_positives", "threshold_positive",
                                              "number_of_unlabelled",
                                              "number_of_negatives", "threshold_negative",
                                              "difference_between_matrices",
                                              "mean_absolute_probability")
  stats_internal_iterations_df = rbind(stats_internal_iterations_df, iter_stats_fragment)
  colnames(stats_internal_iterations_df) <- c("method","iteration","internal_iteration", "number_of_positives", "threshold_positive",
                                              "number_of_unlabelled",
                                              "number_of_negatives", "threshold_negative",
                                              "difference_between_matrices",
                                              "mean_absolute_probability")
  print("stats added")
  
  #write.predictions
  unl.pred.file2 = write.predictions(unl.looping.disc, u.prb2, p.prb2, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  pos.pred.file2 = write.predictions(pos.looping.disc, u.prb2, p.prb2, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  ret.pred.file2 = write.predictions(ret.looping.disc, u.prb2, p.prb2, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))    
  
  
  catfile.unl = paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.unl.add = read.table(catfile.unl, header=FALSE,sep="\t",comment.char="", quote="")
  cat.unl = rbind(cat.unl,cat.unl.add)
  catfile.pos = paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.pos.add = read.table(catfile.pos, header=FALSE,sep="\t",comment.char="", quote="")
  cat.pos = rbind(cat.pos,cat.pos.add)
  catfile.ret = paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.ret.add = read.table(catfile.ret, header=FALSE,sep="\t",comment.char="", quote="")
  cat.ret = rbind(cat.ret,cat.ret.add)
  
  df_for_perf_measures = rbind(cat.unl.add, cat.pos.add, cat.ret.add)
  df_for_perf_measures = cbind(df_for_perf_measures, mod.disc.data$class.label[df_for_perf_measures[,1]])
  df_for_perf_measures[,2] = as.integer(df_for_perf_measures[,2] > 0)
  
  
  perf_measures_fragment <- c(method, i, get.perf.measures(df_for_perf_measures[,2], df_for_perf_measures[,3]))
  names(perf_measures_fragment) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                    "recall", "f1_score_like_estimate")
  
  colnames(performance_measures_df) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                        "recall", "f1_score_like_estimate")
  performance_measures_df <- rbind(performance_measures_df,
                                   perf_measures_fragment)
  
  colnames(performance_measures_df) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                        "recall", "f1_score_like_estimate")
  
  print("performance measures added")
  
  
  
  
  print(paste(i, "DONE!"))
  print(cat("\n\n\n\n"))
  
  
  
}


svg(paste(analysis_home_plots, '/ret_pos_unl_distro_method_corrected_', method,'.svg', sep = ""))
par(mfrow = c(3,1))
# xlim = c(1.00, 1000)
xlim = c(-100, 100)
hist.looping(cat.pos[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for positive uORFs", col="lightblue")
hist.looping(cat.unl[,2], xlim[1], xlim[2], breaks=40, main="Class probability ratio distribution for unlabelled uORFs", col="lightgreen")
hist.looping(cat.ret[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for retained uORFs")
dev.off()

# Some statistics, resulting from the looping algorithm:

quantile(cat.unl[,2], c(.10, .25, .50, .75, .90))
quantile(cat.pos[,2], c(.10, .25, .50, .75, .90))
quantile(cat.ret[,2], c(.10, .25, .50, .75, .90))






method = 3
cat.unl = list()
cat.pos = list()
cat.ret = list()
for (i in 1:8) {
  
  unl.looping.disc = bound[[i]]
  pos.looping.disc = positives[[i]]
  ret.looping.disc = retained[[i]]
  
  quantile_threshold3 = 0.97
  update_positive_score_threshold3 = 10
  update_negative_score_threshold3 = -update_positive_score_threshold3
  
  prob_matrices_list3 = train.matrices.updating.negative(no_label_data, unl.looping.disc, pos.looping.disc,
                                                         max.len, prior,
                                                         #it_min = 4,
                                                         it_max = 50,
                                                         quantile_threshold_positive = quantile_threshold3,
                                                         update_positive_score_threshold = update_positive_score_threshold3,
                                                         quantile_threshold_negative = quantile_threshold3,
                                                         update_negative_score_threshold = update_negative_score_threshold3,
                                                         filename=paste(paste(analysis_home_prob_matrices,
                                                                              "/uORF_probability_matrices_updating_negative", sep = ""),
                                                                        i, prior, method, quantile_threshold3,
                                                                        update_positive_score_threshold3, sep='_'))
  print("Matrices 3 computed!")
  
  
  u.prb3 = prob_matrices_list3[[1]]
  p.prb3 = prob_matrices_list3[[2]]
  iter_stats_fragment = cbind(method, i, prob_matrices_list3[[3]])
  colnames(iter_stats_fragment) <- c("method","iteration","internal_iteration", "number_of_positives", "threshold_positive",
                                     "number_of_unlabelled",
                                     "number_of_negatives", "threshold_negative",
                                     "difference_between_matrices",
                                     "mean_absolute_probability")
  
  colnames(stats_internal_iterations_df) <- c("method","iteration","internal_iteration", "number_of_positives", "threshold_positive",
                                              "number_of_unlabelled",
                                              "number_of_negatives", "threshold_negative",
                                              "difference_between_matrices",
                                              "mean_absolute_probability")
  stats_internal_iterations_df = rbind(stats_internal_iterations_df, iter_stats_fragment)
  colnames(stats_internal_iterations_df) <- c("method","iteration","internal_iteration", "number_of_positives", "threshold_positive",
                                              "number_of_unlabelled",
                                              "number_of_negatives", "threshold_negative",
                                              "difference_between_matrices",
                                              "mean_absolute_probability")
  print("stats added")
  
  
  
  #write.predictions
  unl.pred.file3 = write.predictions(unl.looping.disc, u.prb3, p.prb3, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  pos.pred.file3 = write.predictions(pos.looping.disc, u.prb3, p.prb3, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  ret.pred.file3 = write.predictions(ret.looping.disc, u.prb3, p.prb3, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  
  
  catfile.unl = paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.unl.add = read.table(catfile.unl, header=FALSE,sep="\t",comment.char="", quote="")
  cat.unl = rbind(cat.unl,cat.unl.add)
  catfile.pos = paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.pos.add = read.table(catfile.pos, header=FALSE,sep="\t",comment.char="", quote="")
  cat.pos = rbind(cat.pos,cat.pos.add)
  catfile.ret = paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.ret.add = read.table(catfile.ret, header=FALSE,sep="\t",comment.char="", quote="")
  cat.ret = rbind(cat.ret,cat.ret.add)
  
  
  df_for_perf_measures = rbind(cat.unl.add, cat.pos.add, cat.ret.add)
  df_for_perf_measures = cbind(df_for_perf_measures, mod.disc.data$class.label[df_for_perf_measures[,1]])
  df_for_perf_measures[,2] = as.integer(df_for_perf_measures[,2] > 0)
  
  
  perf_measures_fragment <- c(method, i, get.perf.measures(df_for_perf_measures[,2], df_for_perf_measures[,3]))
  names(perf_measures_fragment) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                    "recall", "f1_score_like_estimate")
  
  colnames(performance_measures_df) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                        "recall", "f1_score_like_estimate")
  performance_measures_df <- rbind(performance_measures_df,
                                   perf_measures_fragment)
  
  colnames(performance_measures_df) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                        "recall", "f1_score_like_estimate")
  
  print("performance measures added")
  
  
  print(paste(i, "DONE!"))
  print(cat("\n\n\n\n"))
  
  
}


svg(paste(analysis_home_plots, '/ret_pos_unl_distro_method_corrected_', method,'.svg', sep = ""))
par(mfrow = c(3,1))
# xlim = c(1.00, 1000)
xlim = c(-100, 100)
hist.looping(cat.pos[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for positive uORFs", col="lightblue")
hist.looping(cat.unl[,2], xlim[1], xlim[2], breaks=40, main="Class probability ratio distribution for unlabelled uORFs", col="lightgreen")
hist.looping(cat.ret[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for retained uORFs")
dev.off()

# Some statistics, resulting from the looping algorithm:

quantile(cat.unl[,2], c(.10, .25, .50, .75, .90))
quantile(cat.pos[,2], c(.10, .25, .50, .75, .90))
quantile(cat.ret[,2], c(.10, .25, .50, .75, .90))





method = 4
cat.unl = list()
cat.pos = list()
cat.ret = list()
for (i in 1:8) {
  
  unl.looping.disc = bound[[i]]
  pos.looping.disc = positives[[i]]
  ret.looping.disc = retained[[i]]
  
  quantile_threshold4 = 0.99
  update_positive_score_threshold4 = 10
  
  prob_matrices_list4 = train.matrices.updating.negative.positive(no_label_data, unl.looping.disc, pos.looping.disc,
                                                                  max.len, prior,
                                                                  #it_min = 4,
                                                                  it_max = 50,
                                                                  quantile_threshold_positive = 0.99, # 0.99
                                                                  update_positive_score_threshold = 5,
                                                                  quantile_threshold_negative = 0.95, # 0.95
                                                                  update_negative_score_threshold = -10,
                                                                  skip_differences = T,
                                                                  filename=paste(paste(analysis_home_prob_matrices,
                                                                                       "/uORF_probability_matrices_updating_negative_and_positive", sep = ""),
                                                                                 i, prior, method, quantile_threshold4, update_positive_score_threshold4, sep='_'))
  print("Matrices 4 computed!")
  
  
  u.prb4 = prob_matrices_list4[[1]]
  p.prb4 = prob_matrices_list4[[2]]
  iter_stats_fragment = cbind(method, i, prob_matrices_list4[[3]])
  colnames(iter_stats_fragment) <- c("method","iteration","internal_iteration", "number_of_positives", "threshold_positive",
                                     "number_of_unlabelled",
                                     "number_of_negatives", "threshold_negative",
                                     "difference_between_matrices",
                                     "mean_absolute_probability")
  
  colnames(stats_internal_iterations_df) <- c("method","iteration","internal_iteration", "number_of_positives", "threshold_positive",
                                              "number_of_unlabelled",
                                              "number_of_negatives", "threshold_negative",
                                              "difference_between_matrices",
                                              "mean_absolute_probability")
  stats_internal_iterations_df = rbind(stats_internal_iterations_df, iter_stats_fragment)
  colnames(stats_internal_iterations_df) <- c("method","iteration","internal_iteration", "number_of_positives", "threshold_positive",
                                              "number_of_unlabelled",
                                              "number_of_negatives", "threshold_negative",
                                              "difference_between_matrices",
                                              "mean_absolute_probability")
  print("stats added")
  
  
  
  #write.predictions
  unl.pred.file4 = write.predictions(unl.looping.disc, u.prb4, p.prb4, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  pos.pred.file4 = write.predictions(pos.looping.disc, u.prb4, p.prb4, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  ret.pred.file4 = write.predictions(ret.looping.disc, u.prb4, p.prb4, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_'))
  
  
  catfile.unl = paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.unl.add = read.table(catfile.unl, header=FALSE,sep="\t",comment.char="", quote="")
  cat.unl = rbind(cat.unl,cat.unl.add)
  catfile.pos = paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.pos.add = read.table(catfile.pos, header=FALSE,sep="\t",comment.char="", quote="")
  cat.pos = rbind(cat.pos,cat.pos.add)
  catfile.ret = paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
  cat.ret.add = read.table(catfile.ret, header=FALSE,sep="\t",comment.char="", quote="")
  cat.ret = rbind(cat.ret,cat.ret.add)
  
  
  df_for_perf_measures = rbind(cat.unl.add, cat.pos.add, cat.ret.add)
  df_for_perf_measures = cbind(df_for_perf_measures, mod.disc.data$class.label[df_for_perf_measures[,1]])
  df_for_perf_measures[,2] = as.integer(df_for_perf_measures[,2] > 0)
  
  
  perf_measures_fragment <- c(method, i, get.perf.measures(df_for_perf_measures[,2], df_for_perf_measures[,3]))
  names(perf_measures_fragment) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                    "recall", "f1_score_like_estimate")
  
  colnames(performance_measures_df) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                        "recall", "f1_score_like_estimate")
  performance_measures_df <- rbind(performance_measures_df,
                                   perf_measures_fragment)
  
  colnames(performance_measures_df) = c("method","iteration","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                        "recall", "f1_score_like_estimate")
  
  print("performance measures added")
  
  print(paste(i, "DONE!"))
  print(cat("\n\n\n\n"))
  
  
}


svg(paste(analysis_home_plots, '/ret_pos_unl_distro_method_corrected_', method,'.svg', sep = ""))
par(mfrow = c(3,1))
# xlim = c(1.00, 1000)
xlim = c(-100, 100)
hist.looping(cat.pos[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for positive uORFs", col="lightblue")
hist.looping(cat.unl[,2], xlim[1], xlim[2], breaks=40, main="Class probability ratio distribution for unlabelled uORFs", col="lightgreen")
hist.looping(cat.ret[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for retained uORFs")
dev.off()

# Some statistics, resulting from the looping algorithm:

quantile(cat.unl[,2], c(.10, .25, .50, .75, .90))
quantile(cat.pos[,2], c(.10, .25, .50, .75, .90))
quantile(cat.ret[,2], c(.10, .25, .50, .75, .90))










# Summarize the stats in plots
svg(paste(analysis_home_plots, '/stats_internal_iterations_by_method.svg', sep = ""), width = 18, height = 16)
par(mfrow = c(1,1))
plot(stats_internal_iterations_df, col = as.factor(stats_internal_iterations_df$method))
dev.off()

svg(paste(analysis_home_plots, '/stats_performance_by_method.svg', sep = ""), width = 18, height = 16)
par(mfrow = c(1,1))
plot(performance_measures_df, col = as.factor(performance_measures_df$method), pch = 19)
dev.off()
  

#}

# # #Concatenate the files together:
# 
# # setwd("~/Sandbox/")
# 
# 
# hist.looping = function(pred, xmin, xmax, breaks=100, main="No main provided", xlab="Positive/unlabelled probability ratio\n prior positive distribution is 0.5", col="lightgreen") {
#   hist(pred[pred > xmin & pred < xmax], breaks=breaks, main=main, xlab=xlab, col=col)
# }
# 
# for (method in c(0,1,2)){
#   
#   cat.unl = list()
#   cat.pos = list()
#   cat.ret = list()
#   
#   for (i in 1:8) {
#     catfile.unl = paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
#     cat.unl.add = read.table(catfile.unl, header=FALSE,sep="\t",comment.char="", quote="")
#     cat.unl = rbind(cat.unl,cat.unl.add)
#     catfile.pos = paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
#     cat.pos.add = read.table(catfile.pos, header=FALSE,sep="\t",comment.char="", quote="")
#     cat.pos = rbind(cat.pos,cat.pos.add)
#     catfile.ret = paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
#     cat.ret.add = read.table(catfile.ret, header=FALSE,sep="\t",comment.char="", quote="")
#     cat.ret = rbind(cat.ret,cat.ret.add)
#     
#     
#   }
#   
#   
#   svg(paste(analysis_home_plots, '/ret_pos_unl_distro_method_corrected_', method,'.svg', sep = ""))
#   par(mfrow = c(3,1))
#   # xlim = c(1.00, 1000)
#   xlim = c(-100, 100)
#   hist.looping(cat.pos[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for positive uORFs", col="lightblue")
#   hist.looping(cat.unl[,2], xlim[1], xlim[2], breaks=40, main="Class probability ratio distribution for unlabelled uORFs", col="lightgreen")
#   hist.looping(cat.ret[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for retained uORFs")
#   dev.off()
#   
#   # Some statistics, resulting from the looping algorithm:
#   
#   quantile(cat.unl[,2], c(.10, .25, .50, .75, .90))
#   quantile(cat.pos[,2], c(.10, .25, .50, .75, .90))
#   quantile(cat.ret[,2], c(.10, .25, .50, .75, .90))
#   
#   
# }
# 





evaluate.diff.priors <- function(data_with_labels, unlabelled_prob, positive_prob){
  scores_diff_prior = list()
  for (pri in c(seq(0.05, 0.95, 0.1))){
    scores_diff_prior = rbind(scores_diff_prior, c(pri,n.Bayes(data_with_labels, u.prb = unlabelled_prob, p.prb = positive_prob, prior = pri)))
  }
  return(scores_diff_prior)
}


evaluate.diff.methods <- function(data_with_labels, unlabelled_prob, positive_prob){
  scores_diff_method = list()
  compute.prob.matrix(compute_prob_unlabelled = , compute_prob_positive = , max.len = max.len)
  
  scores_diff_method = rbind(scores_diff_method, c(pri,n.Bayes(data_with_labels, u.prb = unlabelled_prob, p.prb = positive_prob, prior = pri)))
  
  return(scores_diff_prior)
}






hist.ratio = function(pred, xmin, xmax, breaks=100) {
  hist(pred[pred > xmin & pred < xmax], breaks=breaks)
  
}


hist.trial = function(pred, xmin, xmax, breaks=30, main="No main provided", xlab="value", col="lightgreen") {
  hist(pred[pred > xmin & pred < xmax], breaks=breaks, main=main, xlab=xlab, col=col)
  
}






# Imaging the table ------------------------------------

#The names of the features before removal of homogenous features are written in ./featuresbefore.txt

allfeaturenamesbefore = colnames(mod.disc.data)
write.table(allfeaturenamesbefore, file=paste("./output/featuresbefore.txt", sep=""), quote=FALSE, sep="\t",
            col.names=FALSE, row.names=FALSE)

#The data from each of the discretized features (both unl and pos combined) is plotted next, on an individual basis.

xlim = c(0, 10000)
xlab = paste("trial", sep="")
tracker = 1

for (aname in allfeaturenamesbefore) {
  
  #samplename = paste("./plots/", allfeaturenames[k], ".pdf", sep="")
  pdf(paste("./plots/var_histogram/", aname, ".pdf", sep=""))
  hist.trial(mod.disc.data[,tracker], xlim[1], xlim[2], breaks=30, main=aname, xlab, col="lightblue")
  dev.off()
  tracker <- tracker + 1
  
}

# Imaging the table-------------------------------------

print("Removing homogenous features...")

# this is not working, but we can skip it because there is variation in all variables
mod.disc.data = screen(mod.disc.data)			# Throw out homogenous features (all 1's or all 2's)

write.table(mod.disc.data, file=paste("./output/discretized_data_all_1_screened.txt", sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)

#The names of each feature, retained following discretization, are written to the file ./featuresafter.txt

allfeaturenamesafter = colnames(mod.disc.data)
write.table(allfeaturenamesafter, file=paste("./output/featuresafter.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

#the names of the rejected features (homogenous features) are written to the file ./rejected.txt

rejected = setdiff(allfeaturenamesbefore, allfeaturenamesafter)
write.table(rejected, file=paste("./output/rejected.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

print("Homogenous features rejected.")

#the dimensions of the discretized data are measured

d.M = dim(mod.disc.data)[1]
d.N = dim(mod.disc.data)[2]

#the positive and unlabelled data sets are separated based on the class column, the class column is removed.

all_pos.disc = (mod.disc.data[mod.disc.data[,d.N] == 1,])[, -d.N]	# shaves off the class column
all_unl.disc = (mod.disc.data[mod.disc.data[,d.N] == 0,])[, -d.N]	# shaves off the class column



#retaining a number of the examples from the positive set, to test the quality of the algorithm:

d.p.M = dim(all_pos.disc)[1] #number of rows in the positive example.
d.p.N = dim(all_pos.disc)[2] #number of columns in the positive example.
tenpercentsamplesize = trunc(d.p.M/10)
tenpercentremainder = d.p.M %% 10

# In the following section of code, looping over the data set is accomplished, to make multiple retained sets,
# and multiple products.

#First we randomly sample the positive data:

random_max = length(all_pos.disc[,1])
nonrandom_vector = seq(1,random_max,1)
random_vector = sample(nonrandom_vector)
all_pos.disc.mixed = all_pos.disc[random_vector,] # we shuffle the positive data

retained = list()
for (i in 0:9) {
  retained[[(i+1)]] = all_pos.disc.mixed[(i*tenpercentsamplesize+1):((i+1)*tenpercentsamplesize),1:d.p.N]
}

if (tenpercentremainder == 0) {
  positivestail = list()
  positivestail[[10]] = list()
  for (i in 1:9) {
    positivestail[[(i)]] = all_pos.disc.mixed[(i*tenpercentsamplesize+1):d.p.M,1:d.p.N]
  }
} else {
  positivestail = list()
  for (i in 1:10) {
    positivestail[[(i)]] = all_pos.disc.mixed[(i*tenpercentsamplesize+1):d.p.M,1:d.p.N]
  }
}

positiveshead = list()
positiveshead[[(1)]] = list()
for (i in 2:10) {
  positiveshead[[(i)]] = all_pos.disc.mixed[1:((i-1)*tenpercentsamplesize),1:d.p.N]
}

positives = list()
for (i in 1:10) {
  positives[[(i)]] = rbind(positiveshead[[i]],positivestail[[i]])
}


#now add the positive test group, back to the unlabelled set, for later retrieval
bound = list()
for (i in 1:10) {
  bound[[i]] = rbind(retained[[i]],all_unl.disc)
}

# #-----------------------------------------------
# # The following will be looped:
# 
# #The dimensions of both the unlabelled, and positive, matrices of discretized values are measured.
# setwd(analysis_home_out_files_test)
# 
# #The maximum number of values for a variable is measured.
# max.len = 1
# for (k in 1:d.N) {
#   len = length(table(mod.disc.data[,k]))
#   if (len > max.len) max.len <- len
# }
# print(max.len)
# 
# #for (prior in c(seq(0.05, 0.95, 0.05))){
# #for (prior in c(seq(0.25, 0.95, 0.1))){
#   prior = 0.5
#   for (i in 1:8) {
# 
#     unl.looping.disc = bound[[i]]
#     pos.looping.disc = positives[[i]]
#     ret.looping.disc = retained[[i]]
#     
#     #The number of variables is measured.
#     d.N = dim(disc.data)[2]
#     max.len = 1
#     for (k in 1:d.N) {
#       len = length(table(disc.data[,k]))
#       if (len > max.len) max.len <- len
#     }
#     # print(max.len)
#     
#     
#     prob_matrices_list0 = compute.prob.matrix(unl.looping.disc, pos.looping.disc, max.len)
#     print("Matrices 0 computed!")
#     
#     quantile_threshold1 = 0.95
#     update_positive_score_threshold1 = 10
#     prob_matrices_list1 = train.matrices.updating.unlabelled(unl.looping.disc, pos.looping.disc,
#                                                             max.len, prior,
#                                                             it_min = 4,
#                                                             it_max = 50,
#                                                             quantile_threshold = quantile_threshold1,
#                                                             update_positive_score_threshold = update_positive_score_threshold1,
#                                                             filename=paste(paste(analysis_home_prob_matrices, "/uORF_probability_matrices_updating_unlabelled", sep = ""), i, prior, "1", quantile_threshold1, update_positive_score_threshold1, sep='_'))
#     print("Matrices 1 computed!")
#     
# 
#     quantile_threshold2 = 0.99
#     update_positive_score_threshold2 = 10
#     
#     prob_matrices_list2 = train.matrices.updating.unlabelled.positive(unl.looping.disc, pos.looping.disc,
#                                                                      max.len, prior,
#                                                                      it_min = 5,
#                                                                      it_max = 50,
#                                                                      quantile_threshold = quantile_threshold2,
#                                                                      update_positive_score_threshold = update_positive_score_threshold2,
#                                                                      skip_differences = F,
#                                                                      filename=paste(paste(analysis_home_prob_matrices, "/uORF_probability_matrices_updating_unlabelled_and_positive", sep = ""), i, prior, "2", quantile_threshold2, update_positive_score_threshold2, sep='_'))
#     print("Matrices 2 computed!")
#     
# 
#     u.prb0 = prob_matrices_list0[[1]]
#     p.prb0 = prob_matrices_list0[[2]]
#     
#     u.prb1 = prob_matrices_list1[[1]]
#     p.prb1 = prob_matrices_list1[[2]]
#     
#     u.prb2 = prob_matrices_list2[[1]]
#     p.prb2 = prob_matrices_list2[[2]]
#     
#     
#     #write.predictions
#     unl.pred.file0 = write.predictions(unl.looping.disc, u.prb0, p.prb0, prior,
#                                       filename=paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, "0", sep='_'))
#     pos.pred.file0 = write.predictions(pos.looping.disc, u.prb0, p.prb0, prior,
#                                       filename=paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior,"0", sep='_'))
#     ret.pred.file0 = write.predictions(ret.looping.disc, u.prb0, p.prb0, prior,
#                                       filename=paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior,"0", sep='_'))
#     
#     
#     unl.pred.file1 = write.predictions(unl.looping.disc, u.prb1, p.prb1, prior,
#                                       filename=paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, "1", sep='_'))
#     pos.pred.file1 = write.predictions(pos.looping.disc, u.prb1, p.prb1, prior,
#                                       filename=paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior,"1", sep='_'))
#     ret.pred.file1 = write.predictions(ret.looping.disc, u.prb1, p.prb1, prior,
#                                       filename=paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior,"1", sep='_'))
#     
#     unl.pred.file2 = write.predictions(unl.looping.disc, u.prb2, p.prb2, prior,
#                                       filename=paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, "2", sep='_'))
#     pos.pred.file2 = write.predictions(pos.looping.disc, u.prb2, p.prb2, prior,
#                                       filename=paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior,"2", sep='_'))
#     ret.pred.file2 = write.predictions(ret.looping.disc, u.prb2, p.prb2, prior,
#                                       filename=paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior,"2", sep='_'))    
# 
#     print(paste(i, "DONE!"))
#     print(cat("\n\n\n\n"))
#   }
# #}
#   
# # #Concatenate the files together:
# 
# # setwd("~/Sandbox/")
# 
#   
# hist.looping = function(pred, xmin, xmax, breaks=100, main="No main provided", xlab="Positive/unlabelled probability ratio\n prior positive distribution is 0.5", col="lightgreen") {
#   hist(pred[pred > xmin & pred < xmax], breaks=breaks, main=main, xlab=xlab, col=col)
# }
#   
# for (method in c(0,1,2)){
# 
#   cat.unl = list()
#   cat.pos = list()
#   cat.ret = list()
#   
#   for (i in 1:8) {
#     catfile.unl = paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
#     cat.unl.add = read.table(catfile.unl, header=FALSE,sep="\t",comment.char="", quote="")
#     cat.unl = rbind(cat.unl,cat.unl.add)
#     catfile.pos = paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
#     cat.pos.add = read.table(catfile.pos, header=FALSE,sep="\t",comment.char="", quote="")
#     cat.pos = rbind(cat.pos,cat.pos.add)
#     catfile.ret = paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
#     cat.ret.add = read.table(catfile.ret, header=FALSE,sep="\t",comment.char="", quote="")
#     cat.ret = rbind(cat.ret,cat.ret.add)
#   }
#   
#   
#   svg(paste(analysis_home_plots, '/ret_pos_unl_distro_method_second_round_', method,'.svg', sep = ""))
#   par(mfrow = c(3,1))
#   # xlim = c(1.00, 1000)
#   xlim = c(-100, 100)
#   hist.looping(cat.pos[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for positive uORFs", col="lightblue")
#   hist.looping(cat.unl[,2], xlim[1], xlim[2], breaks=40, main="Class probability ratio distribution for unlabelled uORFs", col="lightgreen")
#   hist.looping(cat.ret[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for retained uORFs")
#   dev.off()
#   
#   # Some statistics, resulting from the looping algorithm:
#   
#   quantile(cat.unl[,2], c(.10, .25, .50, .75, .90))
#   quantile(cat.pos[,2], c(.10, .25, .50, .75, .90))
#   quantile(cat.ret[,2], c(.10, .25, .50, .75, .90))
#   
#   
# }
# 
# 
# 
# # Compare distributions after the reassigning of classes
# 
# ddd <- rbind(cat.pos, cat.unl)
# positive_names <- ddd$V1[which(ddd$V2 > 0)]
# negative_names <- ddd$V1[which(ddd$V2 <= 0)]
# pos.data.after <- no_label_data[positive_names,]
# neg.data.after <- no_label_data[negative_names,]
# 
# 
# 
# #A matrix is constructed, allowing for measurement of the t-test value for each variable.
# allfeaturenamesafter <- names(mod.disc.data)
# 
# # all_pos.disc <- all.disc.corrected[,-49][all.disc.corrected[,49] == 1,]
# # all_unl.disc <- all.disc.corrected[,-49][all.disc.corrected[,49] == 0,]
# 
# ttestmat <- mat.or.vec((length(allfeaturenamesafter)-1),5)
# dimnames(ttestmat)[1] = dimnames(pos.data.after)[2]
# kstestmat <- mat.or.vec((length(allfeaturenamesafter)-1),5)
# dimnames(kstestmat)[1] = dimnames(pos.data.after)[2]
# tracker=1
# for (name in allfeaturenamesafter[1:(length(allfeaturenamesafter)-1)]) {
#   print(name)
#   ttest = t.test(pos.data.after[,tracker],neg.data.after[,tracker])
#   kstest = ks.test(pos.data.after[,tracker],neg.data.after[,tracker])
#   ttestmat[tracker,1] <- ttest$statistic
#   ttestmat[tracker,2] <- ttest$parameter
#   ttestmat[tracker,3] <- ttest$p.value
#   ttestmat[tracker,4] <- ttest$conf.int[1]
#   ttestmat[tracker,5] <- ttest$conf.int[2]
#   
#   kstestmat[tracker,1] <- kstest$statistic
#   kstestmat[tracker,2] <- kstest$p.value
#   kstestmat[tracker,3] <- kstest$alternative
#   kstestmat[tracker,4] <- kstest$method
#   kstestmat[tracker,5] <- kstest$data.name
#   tracker <- tracker + 1
# }
# 
# #Order the t-test matrix according to the magnitude of the t-test value, print to ./ttestmatordered.txt
# 
# ttestmatordered=ttestmat[order(abs(ttestmat[,1])),]
# kstestmatordered=kstestmat[order(abs(as.numeric(kstestmat[,1]))),]
# system("mkdir -p ./output/after/")
# write.table(ttestmatordered, file=paste("./output/after/ttestmatordered.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
# write.table(kstestmatordered, file=paste("./output/after/kstestmatordered.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
# 
# #All variables t-test plot
# system("mkdir -p ./plots/after/ttest/")
# pdf(paste("./plots/after/ttest/", "ttestordered", ".pdf", sep=""))
# barplot(abs(ttestmatordered[,1]), main="t-test", horiz=TRUE, names.arg=c(1:dim(ttestmatordered)[1]), las=1)
# dev.off()
# 
# #All variable ks-test plot
# system("mkdir -p ./plots/after/kstest/")
# pdf(paste("./plots/after/kstest/", "kstestordered_original", ".pdf", sep=""))
# barplot(abs(as.numeric(kstestmatordered[,1])), main="ks-test", horiz=TRUE, names.arg=c(1:dim(kstestmatordered)[1]), las=1)
# dev.off()
# 
# #Top 10 variables ks-test plot
# pdf(paste("./plots/after/kstest/", "kstestordered_top10", ".pdf", sep=""))
# par(mar=c(5.1, 13 ,4.1 ,2.1))
# barplot(abs(as.numeric(kstestmatordered[1:10,1])), main="ks-test", horiz=TRUE, names.arg=tail(dimnames(kstestmatordered)[[1]], 10), las=1)
# dev.off()
# 
# 
# 
# 
# #Legend for the t-test table is constructed, so that the numbers on the table can be associated with variables.
# ttestlegend = ttestmatordered[,1:2]
# ttestlegend[,2] = ttestmatordered[,1]
# ttestlegend[,1] = 1:dim(ttestmatordered)[1]
# write.table(ttestlegend, file=paste("./output/after/ttestlegend.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
# write.table(ttestlegend, file=paste("./plots/after/ttestlegend.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
# #---------------------------------------------------
# 
# 
# 
# 
# #The dimensions of both the unlabelled, and positive, matrices of discretized values are measured.
# 
# u.M = dim(neg.data.after)[1]
# print(u.M)
# u.N = dim(neg.data.after)[2]
# print(u.N)
# p.M = dim(pos.data.after)[1]
# print(p.M)
# p.N = dim(pos.data.after)[2]
# print(p.N)
# print(d.N)
# 
# #The maximum number of values for a variable is measured.
# max.len = 1
# for (k in 1:d.N) {
#   len = length(table(mod.disc.data[,k]))
#   if (len > max.len) max.len <- len
# }
# print(max.len)
# 
# 
# #probability matrices are constructed (empty at this stage)
# u.prb = matrix(0, nrow=u.N, ncol=max.len)
# p.prb = matrix(0, nrow=u.N, ncol=max.len)
# dimnames(u.prb)[1] = dimnames(neg.data.after)[2]
# dimnames(p.prb)[1] = dimnames(pos.data.after)[2]
# 
# 
# #numbers are entered into the probability matrices, indicating the probability a given variable will have a given value.
# for (k in 1:u.N) 
#   for (j in 1:max.len) {
#     u.prb[k, j] = mean(neg.data.after[,k] == j) # as the output from the equality is either True or False, computing the mean is the probability
#     p.prb[k, j] = mean(pos.data.after[,k] == j)	
#   }
# 
# 
# #A table with the number of instances of each value of each variable is constructed (both positive and unlabelled)
# u.count = matrix(0, nrow=u.N, ncol=max.len)
# p.count = matrix(0, nrow=u.N, ncol=max.len)
# dimnames(u.count)[1] = dimnames(neg.data.after)[2]
# dimnames(p.count)[1] = dimnames(pos.data.after)[2]
# 
# 
# #A table with the frequency (percentwise) of each value of for each variable is constructed (both positive and unlabelled)
# u.cumulfreq = matrix(0, nrow=u.N, ncol=max.len)
# p.cumulfreq = matrix(0, nrow=u.N, ncol=max.len)
# dimnames(u.cumulfreq)[1] = dimnames(neg.data.after)[2]
# dimnames(p.cumulfreq)[1] = dimnames(pos.data.after)[2]
# 
# # We fill the last two sets of tables we created, the counts and the cumulative frequency
# for (k in 1:u.N) {
#   varsp=table(pos.data.after[,k])
#   print(varsp)
#   varsu=table(neg.data.after[,k])
#   print(varsu)
#   varsall=table(c(neg.data.after[,k],pos.data.after[,k]))
#   print(varsall)
#   u.count[k,1:length(varsu)] = table(neg.data.after[,k])
#   p.count[k,1:length(varsp)] = table(pos.data.after[,k])
#   vars=table(neg.data.after[,k])/length(neg.data.after[,k])
#   #print(vars)
#   u.cumulfreq[k,1:length(varsu)] = (table(neg.data.after[,k])/length(neg.data.after[,k]))
#   p.cumulfreq[k,1:length(varsp)] = (table(pos.data.after[,k])/length(pos.data.after[,k]))
# }
# 
# print(p.count)
# 
# print(p.cumulfreq)
# print(u.cumulfreq)
# 
# print(p.cumulfreq[1,1:3])
# print(u.cumulfreq[1,1:3])
# 
# 
# #Plot with differential scores between categories ALL DIFFERENT PLOTS
# tracker = 1
# anames = rownames(p.cumulfreq)
# system("mkdir -p ./plots/after/retained/")
# 
# for (k in 1:dim(p.cumulfreq)[1]) {
#   
#   cumulfreqmatvect = c(p.cumulfreq[k,1:dim(p.cumulfreq)[2]], u.cumulfreq[k,1:dim(u.cumulfreq)[2]])
#   
#   print(cumulfreqmatvect)
#   colnumber = length(cumulfreqmatvect)/2
#   
#   trialmat = matrix(cumulfreqmatvect, nrow=2, ncol=colnumber, byrow = TRUE)
#   
#   trialmatnormalized = trialmat
#   for (index in 1:dim(trialmat)[2]) {
#     trialmatnormalized[,index] = trialmat[,index]/sum(trialmat[,index])
#   }
#   
#   pdf(paste("./plots/after/retained/", anames[k], ".pdf", sep=""))
#   barplot(trialmatnormalized, main=anames[k], xlab="value (L to R, 1 to 2 (or 3))", col=c("darkblue","red"), legend = c("positive", "unlabelled"))
#   dev.off()
#   tracker <- tracker + 1
#   
# }
# 
# 
# 
# # #ALL ON THE SAME FILE
# 
# pdf(paste("./plots/after/retained/", "0ALLPLOT_after", ".pdf", sep=""))
# par(mfrow=c(2,2))
# for (k in 1:dim(p.cumulfreq)[1]) {
#   
#   cumulfreqmatvect = c(p.cumulfreq[k,1:dim(p.cumulfreq)[2]], u.cumulfreq[k,1:dim(u.cumulfreq)[2]])
#   colnumber = length(cumulfreqmatvect)/2
#   trialmat = matrix(cumulfreqmatvect, nrow=2, ncol=colnumber, byrow = TRUE)
#   
#   trialmatnormalized = trialmat
#   for (index in 1:dim(trialmat)[2]) {
#     trialmatnormalized[,index] = trialmat[,index]/sum(trialmat[,index])
#   }
#   
#   barplot(trialmatnormalized, main=anames[k], xlab="value", col=c("darkblue","red"), legend = c("positive", "unlabelled"))
#   tracker <- tracker + 1
#   
# }
# dev.off()
# 



# #-------------------------------------












#A matrix is constructed, allowing for measurement of the t-test value for each variable.
allfeaturenamesafter <- names(mod.disc.data)

# all_pos.disc <- all.disc.corrected[,-49][all.disc.corrected[,49] == 1,]
# all_unl.disc <- all.disc.corrected[,-49][all.disc.corrected[,49] == 0,]

ttestmat <- mat.or.vec((length(allfeaturenamesafter)-1),5)
dimnames(ttestmat)[1] = dimnames(all_pos.disc)[2]
kstestmat <- mat.or.vec((length(allfeaturenamesafter)-1),5)
dimnames(kstestmat)[1] = dimnames(all_pos.disc)[2]
tracker=1
for (name in allfeaturenamesafter[1:(length(allfeaturenamesafter)-1)]) {
  ttest = t.test(all_pos.disc[,tracker],all_unl.disc[,tracker])
  kstest = ks.test(all_pos.disc[,tracker],all_unl.disc[,tracker])
  ttestmat[tracker,1] <- ttest$statistic
  ttestmat[tracker,2] <- ttest$parameter
  ttestmat[tracker,3] <- ttest$p.value
  ttestmat[tracker,4] <- ttest$conf.int[1]
  ttestmat[tracker,5] <- ttest$conf.int[2]
  
  kstestmat[tracker,1] <- kstest$statistic
  kstestmat[tracker,2] <- kstest$p.value
  kstestmat[tracker,3] <- kstest$alternative
  kstestmat[tracker,4] <- kstest$method
  kstestmat[tracker,5] <- kstest$data.name
  tracker <- tracker + 1
}

#Order the t-test matrix according to the magnitude of the t-test value, print to ./ttestmatordered.txt

ttestmatordered=ttestmat[order(abs(ttestmat[,1])),]
kstestmatordered=kstestmat[order(abs(as.numeric(kstestmat[,1]))),]
write.table(ttestmatordered, file=paste("./output/ttestmatordered.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
write.table(kstestmatordered, file=paste("./output/kstestmatordered.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)


#All variables t-test plot
system("mkdir -p ./plots/ttest/")
pdf(paste("./plots/ttest/", "ttestordered", ".pdf", sep=""))
barplot(abs(ttestmatordered[,1]), main="t-test", horiz=TRUE, names.arg=c(1:dim(ttestmatordered)[1]), las=1)
dev.off()

#All variable ks-test plot
system("mkdir -p ./plots/kstest/")
pdf(paste("./plots/kstest/", "kstestordered_original", ".pdf", sep=""))
barplot(abs(as.numeric(kstestmatordered[,1])), main="ks-test", horiz=TRUE, names.arg=c(1:dim(kstestmatordered)[1]), las=1)
dev.off()

#Top 10 variables ks-test plot
pdf(paste("./plots/kstest/", "kstestordered_top10", ".pdf", sep=""))
par(mar=c(5.1, 13 ,4.1 ,2.1))
barplot(abs(as.numeric(kstestmatordered[(nrow(kstestmatordered)-9):nrow(kstestmatordered),1])), main="ks-test", horiz=TRUE, names.arg=tail(dimnames(kstestmatordered)[[1]], 10), las=1)
dev.off()




#Legend for the t-test table is constructed, so that the numbers on the table can be associated with variables.
ttestlegend = ttestmatordered[,1:2]
ttestlegend[,2] = ttestmatordered[,1]
ttestlegend[,1] = 1:dim(ttestmatordered)[1]
write.table(ttestlegend, file=paste("./output/ttestlegend.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
write.table(ttestlegend, file=paste("./plots/ttestlegend.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
#---------------------------------------------------




#The dimensions of both the unlabelled, and positive, matrices of discretized values are measured.

u.M = dim(all_unl.disc)[1]
print(u.M)
u.N = dim(all_unl.disc)[2]
print(u.N)
p.M = dim(all_pos.disc)[1]
print(p.M)
p.N = dim(all_pos.disc)[2]
print(p.N)
print(d.N)

#The maximum number of values for a variable is measured.
max.len = 1
for (k in 1:d.N) {
  len = length(table(mod.disc.data[,k]))
  if (len > max.len) max.len <- len
}
print(max.len)


#probability matrices are constructed (empty at this stage)
u.prb = matrix(0, nrow=u.N, ncol=max.len)
p.prb = matrix(0, nrow=u.N, ncol=max.len)
dimnames(u.prb)[1] = dimnames(all_unl.disc)[2]
dimnames(p.prb)[1] = dimnames(all_pos.disc)[2]


#numbers are entered into the probability matrices, indicating the probability a given variable will have a given value.
for (k in 1:u.N) 
  for (j in 1:max.len) {
    u.prb[k, j] = mean(all_unl.disc[,k] == j) # as the output from the equality is either True or False, computing the mean is the probability
    p.prb[k, j] = mean(all_pos.disc[,k] == j)	
  }


#A table with the number of instances of each value of each variable is constructed (both positive and unlabelled)
u.count = matrix(0, nrow=u.N, ncol=max.len)
p.count = matrix(0, nrow=u.N, ncol=max.len)
dimnames(u.count)[1] = dimnames(all_unl.disc)[2]
dimnames(p.count)[1] = dimnames(all_pos.disc)[2]


#A table with the frequency (percentwise) of each value of for each variable is constructed (both positive and unlabelled)
u.cumulfreq = matrix(0, nrow=u.N, ncol=max.len)
p.cumulfreq = matrix(0, nrow=u.N, ncol=max.len)
dimnames(u.cumulfreq)[1] = dimnames(all_unl.disc)[2]
dimnames(p.cumulfreq)[1] = dimnames(all_pos.disc)[2]

# We fill the last two sets of tables we created, the counts and the cumulative frequency
for (k in 1:u.N) {
  varsp=table(all_pos.disc[,k])
  print(varsp)
  varsu=table(all_unl.disc[,k])
  print(varsu)
  varsall=table(c(all_unl.disc[,k],all_pos.disc[,k]))
  print(varsall)
  u.count[k,1:length(varsu)] = table(all_unl.disc[,k])
  p.count[k,1:length(varsp)] = table(all_pos.disc[,k])
  vars=table(all_unl.disc[,k])/length(all_unl.disc[,k])
  #print(vars)
  u.cumulfreq[k,1:length(varsu)] = (table(all_unl.disc[,k])/length(all_unl.disc[,k]))
  p.cumulfreq[k,1:length(varsp)] = (table(all_pos.disc[,k])/length(all_pos.disc[,k]))
}

print(p.count)

print(p.cumulfreq)
print(u.cumulfreq)

print(p.cumulfreq[1,1:3])
print(u.cumulfreq[1,1:3])


#Plot with differential scores between categories ALL DIFFERENT PLOTS
tracker = 1
anames = rownames(p.cumulfreq)
system("mkdir -p ./plots/retained/")

for (k in 1:dim(p.cumulfreq)[1]) {
  
  cumulfreqmatvect = c(p.cumulfreq[k,1:dim(p.cumulfreq)[2]], u.cumulfreq[k,1:dim(u.cumulfreq)[2]])
  
  print(cumulfreqmatvect)
  colnumber = length(cumulfreqmatvect)/2
  
  trialmat = matrix(cumulfreqmatvect, nrow=2, ncol=colnumber, byrow = TRUE)
  
  trialmatnormalized = trialmat
  for (index in 1:dim(trialmat)[2]) {
    trialmatnormalized[,index] = trialmat[,index]/sum(trialmat[,index])
  }
  
  pdf(paste("./plots/retained/", anames[k], ".pdf", sep=""))
  barplot(trialmatnormalized, main=anames[k], xlab="value (L to R, 1 to 2 (or 3))", col=c("darkblue","red"), legend = c("positive", "unlabelled"))
  dev.off()
  tracker <- tracker + 1
  
}



# #ALL ON THE SAME FILE

pdf(paste("./plots/retained/", "0ALLPLOT", ".pdf", sep=""))
par(mfrow=c(2,2))
for (k in 1:dim(p.cumulfreq)[1]) {
  
  cumulfreqmatvect = c(p.cumulfreq[k,1:dim(p.cumulfreq)[2]], u.cumulfreq[k,1:dim(u.cumulfreq)[2]])
  colnumber = length(cumulfreqmatvect)/2
  trialmat = matrix(cumulfreqmatvect, nrow=2, ncol=colnumber, byrow = TRUE)
  
  trialmatnormalized = trialmat
  for (index in 1:dim(trialmat)[2]) {
    trialmatnormalized[,index] = trialmat[,index]/sum(trialmat[,index])
  }
  
  barplot(trialmatnormalized, main=anames[k], xlab="value", col=c("darkblue","red"), legend = c("positive", "unlabelled"))
  tracker <- tracker + 1
  
}
dev.off()













# Make ttest, kstest and everything after the classification and the reassignment of classes


for (method in 0:4){
  prior = 0.5
  cat.unl = list()
  cat.pos = list()
  cat.ret = list()

  for (i in 1:8) {
    catfile.unl = paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
    cat.unl.add = read.table(catfile.unl, header=FALSE,sep="\t",comment.char="", quote="")
    cat.unl = rbind(cat.unl,cat.unl.add)
    catfile.pos = paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
    cat.pos.add = read.table(catfile.pos, header=FALSE,sep="\t",comment.char="", quote="")
    cat.pos = rbind(cat.pos,cat.pos.add)
    catfile.ret = paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior, method, sep='_')
    cat.ret.add = read.table(catfile.ret, header=FALSE,sep="\t",comment.char="", quote="")
    cat.ret = rbind(cat.ret,cat.ret.add)
  }
  
  
  
  # Compare distributions after the reassigning of classes
  ddd <- rbind(cat.pos, cat.unl) # we leave the unlabelled ones out so that we don't have them repeated in both datasets
  positive_names <- ddd$V1[which(ddd$V2 > 0)]
  negative_names <- ddd$V1[which(ddd$V2 <= 0)]
  pos.data.after <- no_label_data[positive_names,]
  neg.data.after <- no_label_data[negative_names,]
  
  
  
  #A matrix is constructed, allowing for measurement of the t-test value for each variable.
  allfeaturenamesafter <- names(mod.disc.data)
  
  # all_pos.disc <- all.disc.corrected[,-49][all.disc.corrected[,49] == 1,]
  # all_unl.disc <- all.disc.corrected[,-49][all.disc.corrected[,49] == 0,]
  
  ttestmat <- mat.or.vec((length(allfeaturenamesafter)-1),5)
  dimnames(ttestmat)[1] = dimnames(pos.data.after)[2]
  kstestmat <- mat.or.vec((length(allfeaturenamesafter)-1),5)
  dimnames(kstestmat)[1] = dimnames(pos.data.after)[2]
  tracker=1
  for (name in allfeaturenamesafter[1:(length(allfeaturenamesafter)-1)]) {
    print(name)
    ttest = t.test(pos.data.after[,tracker],neg.data.after[,tracker])
    kstest = ks.test(pos.data.after[,tracker],neg.data.after[,tracker])
    ttestmat[tracker,1] <- ttest$statistic
    ttestmat[tracker,2] <- ttest$parameter
    ttestmat[tracker,3] <- ttest$p.value
    ttestmat[tracker,4] <- ttest$conf.int[1]
    ttestmat[tracker,5] <- ttest$conf.int[2]
    
    kstestmat[tracker,1] <- kstest$statistic
    kstestmat[tracker,2] <- kstest$p.value
    kstestmat[tracker,3] <- kstest$alternative
    kstestmat[tracker,4] <- kstest$method
    kstestmat[tracker,5] <- kstest$data.name
    tracker <- tracker + 1
  }
  
  #Order the t-test matrix according to the magnitude of the t-test value, print to ./ttestmatordered.txt
  
  ttestmatordered=ttestmat[order(abs(ttestmat[,1])),]
  kstestmatordered=kstestmat[order(abs(as.numeric(kstestmat[,1]))),]
  system(paste("mkdir -p ./output/method", method, sep = ""))
  write.table(ttestmatordered, file=paste("./output/method", method, "/ttestmatordered.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
  write.table(kstestmatordered, file=paste("./output/method", method, "/kstestmatordered.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
  
  #All variables t-test plot
  system(paste("mkdir -p ./plots/method", method, "/", sep = ""))
  system(paste("mkdir -p ./plots/method", method, "/ttest/", sep = ""))
  pdf(paste("./plots/method", method, "/ttest/", "ttestordered", ".pdf", sep=""))
  barplot(abs(ttestmatordered[,1]), main="t-test", horiz=TRUE, names.arg=c(1:dim(ttestmatordered)[1]), las=1)
  dev.off()
  
  #All variable ks-test plot
  system(paste("mkdir -p ./plots/method", method, "/kstest/", sep = ""))
  pdf(paste("./plots/method", method, "/kstest/", "kstestordered_original", ".pdf", sep=""))
  barplot(abs(as.numeric(kstestmatordered[,1])), main="ks-test", horiz=TRUE, names.arg=c(1:dim(kstestmatordered)[1]), las=1)
  dev.off()
  
  #Top 10 variables ks-test plot
  pdf(paste("./plots/method", method, "/kstest/", "kstestordered_top10", ".pdf", sep=""))
  par(mar=c(5.1, 13 ,4.1 ,2.1))
  barplot(abs(as.numeric(kstestmatordered[(nrow(kstestmatordered)-9):nrow(kstestmatordered),1])), main="ks-test", horiz=TRUE, names.arg=tail(dimnames(kstestmatordered)[[1]], 10), las=1)
  dev.off()
  
  
  
  
  #Legend for the t-test table is constructed, so that the numbers on the table can be associated with variables.
  ttestlegend = ttestmatordered[,1:2]
  ttestlegend[,2] = ttestmatordered[,1]
  ttestlegend[,1] = 1:dim(ttestmatordered)[1]
  write.table(ttestlegend, file=paste("./output/method", method,"/ttestlegend.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
  write.table(ttestlegend, file=paste("./plots/method", method, "/ttestlegend.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
  #---------------------------------------------------
  
  
  
  
  #The dimensions of both the unlabelled, and positive, matrices of discretized values are measured.
  
  u.M = dim(neg.data.after)[1]
  print(u.M)
  u.N = dim(neg.data.after)[2]
  print(u.N)
  p.M = dim(pos.data.after)[1]
  print(p.M)
  p.N = dim(pos.data.after)[2]
  print(p.N)
  print(d.N)
  
  #The maximum number of values for a variable is measured.
  max.len = 1
  for (k in 1:d.N) {
    len = length(table(mod.disc.data[,k]))
    if (len > max.len) max.len <- len
  }
  print(max.len)
  
  
  #probability matrices are constructed (empty at this stage)
  u.prb = matrix(0, nrow=u.N, ncol=max.len)
  p.prb = matrix(0, nrow=u.N, ncol=max.len)
  dimnames(u.prb)[1] = dimnames(neg.data.after)[2]
  dimnames(p.prb)[1] = dimnames(pos.data.after)[2]
  
  
  #numbers are entered into the probability matrices, indicating the probability a given variable will have a given value.
  for (k in 1:u.N) 
    for (j in 1:max.len) {
      u.prb[k, j] = mean(neg.data.after[,k] == j) # as the output from the equality is either True or False, computing the mean is the probability
      p.prb[k, j] = mean(pos.data.after[,k] == j)	
    }
  
  
  #A table with the number of instances of each value of each variable is constructed (both positive and unlabelled)
  u.count = matrix(0, nrow=u.N, ncol=max.len)
  p.count = matrix(0, nrow=u.N, ncol=max.len)
  dimnames(u.count)[1] = dimnames(neg.data.after)[2]
  dimnames(p.count)[1] = dimnames(pos.data.after)[2]
  
  
  #A table with the frequency (percentwise) of each value of for each variable is constructed (both positive and unlabelled)
  u.cumulfreq = matrix(0, nrow=u.N, ncol=max.len)
  p.cumulfreq = matrix(0, nrow=u.N, ncol=max.len)
  dimnames(u.cumulfreq)[1] = dimnames(neg.data.after)[2]
  dimnames(p.cumulfreq)[1] = dimnames(pos.data.after)[2]
  
  # We fill the last two sets of tables we created, the counts and the cumulative frequency
  for (k in 1:u.N) {
    varsp=table(pos.data.after[,k])
    print(varsp)
    varsu=table(neg.data.after[,k])
    print(varsu)
    varsall=table(c(neg.data.after[,k],pos.data.after[,k]))
    print(varsall)
    u.count[k,1:length(varsu)] = table(neg.data.after[,k])
    p.count[k,1:length(varsp)] = table(pos.data.after[,k])
    vars=table(neg.data.after[,k])/length(neg.data.after[,k])
    #print(vars)
    u.cumulfreq[k,1:length(varsu)] = (table(neg.data.after[,k])/length(neg.data.after[,k]))
    p.cumulfreq[k,1:length(varsp)] = (table(pos.data.after[,k])/length(pos.data.after[,k]))
  }
  
  print(p.count)
  
  print(p.cumulfreq)
  print(u.cumulfreq)
  
  print(p.cumulfreq[1,1:3])
  print(u.cumulfreq[1,1:3])
  
  
  #Plot with differential scores between categories ALL DIFFERENT PLOTS
  tracker = 1
  anames = rownames(p.cumulfreq)
  system(paste("mkdir -p ./plots/method", method, "/retained", sep = ""))
  
  for (k in 1:dim(p.cumulfreq)[1]) {
    
    cumulfreqmatvect = c(p.cumulfreq[k,1:dim(p.cumulfreq)[2]], u.cumulfreq[k,1:dim(u.cumulfreq)[2]])
    
    print(cumulfreqmatvect)
    colnumber = length(cumulfreqmatvect)/2
    
    trialmat = matrix(cumulfreqmatvect, nrow=2, ncol=colnumber, byrow = TRUE)
    
    trialmatnormalized = trialmat
    for (index in 1:dim(trialmat)[2]) {
      trialmatnormalized[,index] = trialmat[,index]/sum(trialmat[,index])
    }
    
    pdf(paste("./plots/method", method, "/retained/", anames[k], ".pdf", sep=""))
    barplot(trialmatnormalized, main=anames[k], xlab="value (L to R, 1 to 2 (or 3))", col=c("darkblue","red"), legend = c("positive", "unlabelled"))
    dev.off()
    tracker <- tracker + 1
    
  }
  
  
  
  # #ALL ON THE SAME FILE
  
  pdf(paste("./plots/method", method, "/retained/", "0ALLPLOT_after", ".pdf", sep=""))
  par(mfrow=c(2,2))
  for (k in 1:dim(p.cumulfreq)[1]) {
    
    cumulfreqmatvect = c(p.cumulfreq[k,1:dim(p.cumulfreq)[2]], u.cumulfreq[k,1:dim(u.cumulfreq)[2]])
    colnumber = length(cumulfreqmatvect)/2
    trialmat = matrix(cumulfreqmatvect, nrow=2, ncol=colnumber, byrow = TRUE)
    
    trialmatnormalized = trialmat
    for (index in 1:dim(trialmat)[2]) {
      trialmatnormalized[,index] = trialmat[,index]/sum(trialmat[,index])
    }
    
    barplot(trialmatnormalized, main=anames[k], xlab="value", col=c("darkblue","red"), legend = c("positive", "unlabelled"))
    tracker <- tracker + 1
    
  }
  dev.off()
  
  
}










##################################################
########### COMPUTE WITH ALL THE uORFs ###########
##################################################

# store the performance measures at the end of each iteration and model
performance_measures_df_in_complete_df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(performance_measures_df_in_complete_df) = c("method","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                      "recall", "f1_score_like_estimate")


stats_internal_iterations_df_complete_df <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(stats_internal_iterations_df_complete_df) = c("method", "internal_iteration", "number_of_positives", "threshold_positive",
                                           "number_of_unlabelled",
                                           "number_of_negatives", "threshold_negative",
                                           "difference_between_matrices",
                                           "mean_absolute_probability")


internal_internal_stats_df <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(internal_internal_stats_df) = c("internal_iteration", "number_of_positives", "threshold_positive",
                                           "number_of_unlabelled",
                                           "number_of_negatives", "threshold_negative",
                                           "difference_between_matrices",
                                           "mean_absolute_probability")




prior = 0.5

analysis_home_out_files_pred = paste(analysis_home_out_files, "/", Sys.Date(), sep = "")
system(paste("mkdir -p ", analysis_home_out_files_pred, sep = ""))


for (method in 0:4){
  analysis_home_out_files_method = paste(analysis_home_out_files_pred, "/method", method, sep = "")
  system(paste("mkdir -p ", analysis_home_out_files_method, sep = ""))
  
  cat.unl = list()
  cat.pos = list()
  cat.ret = list()
  
  if (method == 0){
    prob_matrices_list_final = compute.prob.matrix(unl.looping.disc, pos.looping.disc, max.len)
    print("Matrices 0 computed!")
    
  } else if (method == 1){
    quantile_threshold1 = 0.97
    update_positive_score_threshold1 = 10
    prob_matrices_list_final = train.matrices.updating.unlabelled(no_label_data, unl.looping.disc, pos.looping.disc,
                                                             max.len, prior,
                                                             #it_min = 1,
                                                             it_max = 50,
                                                             quantile_threshold = quantile_threshold1,
                                                             update_positive_score_threshold = update_positive_score_threshold1,
                                                             filename=paste(paste(analysis_home_prob_matrices, "/uORF_probability_matrices_updating_unlabelled", sep = ""), i,
                                                                            prior, method, quantile_threshold1, update_positive_score_threshold1, sep='_'))
    print("Matrices 1 computed!")
    
    
  } else if (method == 2){
    quantile_threshold2 = 0.97
    update_positive_score_threshold2 = 10
    prob_matrices_list_final = train.matrices.updating.unlabelled.positive(no_label_data,unl.looping.disc, pos.looping.disc,
                                                                      max.len, prior,
                                                                      #it_min = 5,
                                                                      it_max = 50,
                                                                      quantile_threshold = quantile_threshold2,
                                                                      update_positive_score_threshold = update_positive_score_threshold2,
                                                                      #skip_differences = T,
                                                                      filename=paste(paste(analysis_home_prob_matrices, "/uORF_probability_matrices_updating_unlabelled_and_positive", sep = ""), i, prior,  method, quantile_threshold2, update_positive_score_threshold2, sep='_'))
    print("Matrices 2 computed!")
    
  } else if (method == 3){
    
    quantile_threshold3 = 0.97
    update_positive_score_threshold3 = 10
    update_negative_score_threshold3 = -update_positive_score_threshold3
    prob_matrices_list_final = train.matrices.updating.negative(no_label_data, unl.looping.disc, pos.looping.disc,
                                                           max.len, prior,
                                                           #it_min = 4,
                                                           it_max = 50,
                                                           quantile_threshold_positive = quantile_threshold3,
                                                           update_positive_score_threshold = update_positive_score_threshold3,
                                                           quantile_threshold_negative = quantile_threshold3,
                                                           update_negative_score_threshold = update_negative_score_threshold3,
                                                           filename=paste(paste(analysis_home_prob_matrices,
                                                                                "/uORF_probability_matrices_updating_negative", sep = ""),
                                                                          i, prior, method, quantile_threshold3,
                                                                          update_positive_score_threshold3, sep='_'))
    print("Matrices 3 computed!")
    
    
  } else if (method == 4){
    quantile_threshold4 = 0.99
    update_positive_score_threshold4 = 10
    prob_matrices_list_final = train.matrices.updating.negative.positive(no_label_data, unl.looping.disc, pos.looping.disc,
                                                                    max.len, prior,
                                                                    #it_min = 4,
                                                                    it_max = 50,
                                                                    quantile_threshold_positive = 0.99,
                                                                    update_positive_score_threshold = 5,
                                                                    quantile_threshold_negative = 0.95,
                                                                    update_negative_score_threshold = -10,
                                                                    filename=paste(paste(analysis_home_prob_matrices,
                                                                                         "/uORF_probability_matrices_updating_negative_and_positive", sep = ""),
                                                                                   i, prior, method, quantile_threshold4, update_positive_score_threshold4, sep='_'))
    print("Matrices 4 computed!")
    
  }
  
  if (length(prob_matrices_list_final) == 2){
    
    u.prb_final = prob_matrices_list_final[[1]]
    p.prb_final = prob_matrices_list_final[[2]]
    
  } else {
    
    u.prb_final = prob_matrices_list_final[[1]]
    p.prb_final = prob_matrices_list_final[[2]]
    
    iter_stats_fragment = cbind(method, prob_matrices_list_final[[3]])
    colnames(iter_stats_fragment) <- c("method", "internal_iteration", "number_of_positives", "threshold_positive",
                                       "number_of_unlabelled",
                                       "number_of_negatives", "threshold_negative",
                                       "difference_between_matrices",
                                       "mean_absolute_probability")
    
    
    colnames(stats_internal_iterations_df_complete_df) <- c("method","internal_iteration", "number_of_positives", "threshold_positive",
                                                            "number_of_unlabelled",
                                                            "number_of_negatives", "threshold_negative",
                                                            "difference_between_matrices",
                                                            "mean_absolute_probability")
    stats_internal_iterations_df_complete_df = rbind(stats_internal_iterations_df_complete_df, iter_stats_fragment)
    colnames(stats_internal_iterations_df_complete_df) <- c("method","internal_iteration", "number_of_positives", "threshold_positive",
                                                            "number_of_unlabelled",
                                                            "number_of_negatives", "threshold_negative",
                                                            "difference_between_matrices",
                                                            "mean_absolute_probability")
    print("stats added")
    
  }
  
  
  #write.predictions
  unl.pred.file = write.predictions(all_unl.disc, u.prb_final, p.prb_final, prior,
                                     filename=paste(paste(analysis_home_out_files_method, "/uORF_function_unlabelled_predictions_sorted_", sep = ""), prior, method, sep='_'))
  pos.pred.file = write.predictions(all_pos.disc, u.prb_final, p.prb_final, prior,
                                     filename=paste(paste(analysis_home_out_files_method, "/uORF_function_positive_predictions_sorted_", sep = ""), prior, method, sep='_'))
  
  catfile.unl = paste(paste(analysis_home_out_files_method, "/uORF_function_unlabelled_predictions_sorted_", sep = ""), prior, method, sep='_')
  cat.unl.add = read.table(catfile.unl, header=FALSE,sep="\t",comment.char="", quote="")
  catfile.pos = paste(paste(analysis_home_out_files_method, "/uORF_function_positive_predictions_sorted_", sep = ""), prior, method, sep='_')
  cat.pos.add = read.table(catfile.pos, header=FALSE,sep="\t",comment.char="", quote="")
  
  
  df_for_perf_measures = rbind(cat.unl.add, cat.pos.add)
  df_for_perf_measures = cbind(df_for_perf_measures, mod.disc.data$class.label[df_for_perf_measures[,1]])
  df_for_perf_measures[,2] = as.integer(df_for_perf_measures[,2] > 0)
  
  
  perf_measures_fragment <- c(method, get.perf.measures(df_for_perf_measures[,2], df_for_perf_measures[,3]))
  names(perf_measures_fragment) = c("method","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                    "recall", "f1_score_like_estimate")
  
  colnames(performance_measures_df_in_complete_df) = c("method","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                                       "recall", "f1_score_like_estimate")
  performance_measures_df_in_complete_df <- rbind(performance_measures_df_in_complete_df,
                                                  perf_measures_fragment)
  
  colnames(performance_measures_df_in_complete_df) = c("method","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                                       "recall", "f1_score_like_estimate")
  print("performance measures added")
  
  
  print(paste("Method ", method, " DONE!", sep = ""))
  print(cat("\n\n\n\n"))

}


# Summarize the stats in plots
svg(paste(analysis_home_out_files_pred, '/stats_internal_iterations_by_method.svg', sep = ""), width = 18, height = 16)
par(mfrow = c(1,1))
plot(stats_internal_iterations_df_complete_df, col = as.factor(stats_internal_iterations_df_complete_df$method))
dev.off()

svg(paste(analysis_home_out_files_pred, '/stats_performance_by_method.svg', sep = ""), width = 18, height = 16)
par(mfrow = c(1,1))
plot(performance_measures_df_in_complete_df, col = as.factor(performance_measures_df_in_complete_df$method), pch = 19)
dev.off()








# method = 1
# 
# quantile_threshold1 = 0.97
# update_positive_score_threshold1 = 10
# prob_matrices_list1 = train.matrices.updating.unlabelled(no_label_data, unl.looping.disc, pos.looping.disc,
#                                                          max.len, prior,
#                                                          #it_min = 1,
#                                                          it_max = 50,
#                                                          quantile_threshold = quantile_threshold1,
#                                                          update_positive_score_threshold = update_positive_score_threshold1,
#                                                          filename=paste(paste(analysis_home_prob_matrices, "/uORF_probability_matrices_updating_unlabelled", sep = ""), i, prior, method, quantile_threshold1, update_positive_score_threshold1, sep='_'))
# print("Matrices 1 computed!")
# 
# 
# method = 2
# 
# quantile_threshold2 = 0.97
# update_positive_score_threshold2 = 10
# prob_matrices_list2 = train.matrices.updating.unlabelled.positive(no_label_data,unl.looping.disc, pos.looping.disc,
#                                                                   max.len, prior,
#                                                                   #it_min = 5,
#                                                                   it_max = 50,
#                                                                   quantile_threshold = quantile_threshold2,
#                                                                   update_positive_score_threshold = update_positive_score_threshold2,
#                                                                   #skip_differences = T,
#                                                                   filename=paste(paste(analysis_home_prob_matrices, "/uORF_probability_matrices_updating_unlabelled_and_positive", sep = ""), i, prior,  method, quantile_threshold2, update_positive_score_threshold2, sep='_'))
# print("Matrices 2 computed!")
# 
# 
# 
# method = 3
# 
# quantile_threshold3 = 0.97
# update_positive_score_threshold3 = 10
# update_negative_score_threshold3 = -update_positive_score_threshold3
# 
# prob_matrices_list3 = train.matrices.updating.negative(no_label_data, unl.looping.disc, pos.looping.disc,
#                                                        max.len, prior,
#                                                        #it_min = 4,
#                                                        it_max = 50,
#                                                        quantile_threshold_positive = quantile_threshold3,
#                                                        update_positive_score_threshold = update_positive_score_threshold3,
#                                                        quantile_threshold_negative = quantile_threshold3,
#                                                        update_negative_score_threshold = update_negative_score_threshold3,
#                                                        filename=paste(paste(analysis_home_prob_matrices,
#                                                                             "/uORF_probability_matrices_updating_negative", sep = ""),
#                                                                       i, prior, method, quantile_threshold3,
#                                                                       update_positive_score_threshold3, sep='_'))
# print("Matrices 3 computed!")
# 
# 
# 
# method = 4
# 
# quantile_threshold4 = 0.99
# update_positive_score_threshold4 = 10
# prob_matrices_list4 = train.matrices.updating.negative.positive(no_label_data, unl.looping.disc, pos.looping.disc,
#                                                                 max.len, prior,
#                                                                 it_min = 4,
#                                                                 it_max = 50,
#                                                                 quantile_threshold_positive = 0.99,
#                                                                 update_positive_score_threshold = 5,
#                                                                 quantile_threshold_negative = 0.95,
#                                                                 update_negative_score_threshold = -10,
#                                                                 filename=paste(paste(analysis_home_prob_matrices,
#                                                                                      "/uORF_probability_matrices_updating_negative_and_positive", sep = ""),
#                                                                                i, prior, method, quantile_threshold4, update_positive_score_threshold4, sep='_'))
# print("Matrices 4 computed!")
# 
# 

 

paste(paste(analysis_home_out_files_method, "/uORF_function_positive_predictions_sorted_", sep = ""), prior, method, sep='_')
paste(paste(analysis_home_out_files_method, "/uORF_function_unlabelled_predictions_sorted_", sep = ""), prior, method, sep='_')



# 0.53 was a good prior in the cross validation
# 0.61 is the prior I chose to compare normal and corrected discretized predictions
# Ahh, it turns out that in my corrected discretized predictions the CDS.start.to.uORF.end variable
# has all values for positive as 1, thus giving 0 probability for any unlabelled uORF that doesn't meet
# that criteria. While this may be cause for concern and motivate using one of the Naive Bayes
# modifications to keep any probability from being 0, I think ti is just fine to leave out these
# examples. However, I will get an infinity error for dividing by this 0 in outputing the result

setwd("./output/")

prior = 0.5

#write.predictions



#-------------------------------------

# Function Definitions 

hist.ratio = function(pred, xmin, xmax, breaks=100, main="No main provided", xlab="Positive/unlabelled probability ratio\n prior positive distribution is 0.61", col="lightgreen") {
  hist(pred[pred > xmin & pred < xmax], breaks=breaks, main=main, xlab=xlab, col=col)
  
}


#------------------------------------------

# Distribution plots and Box and Whisker plots to evaluate positive predictions

hist.ratio = function(pred, xmin, xmax, breaks=100, label="unlabelled", col='lightblue', prior=0.61) {
  hist(pred[pred > xmin & pred < xmax], xlim=c(xmin, xmax), col=col, breaks=breaks, main=paste("Class probability ratio distribution\n for", label, "uORFs\n predicted as positives"), xlab=paste("Positive/unlabelled probability ratio\n with", prior, "prior distribution"))
  
}

get.pos = function(data) {
  return(data[data[,2] > 1 & !is.na(data[,2]),])
}



#-------------------------------------
# setwd(paste(home, "Sandbox/", sep = ""))
library("ROCR")
require(ggplot2)

for (method in 0:4){
  analysis_home_out_files_method = paste(analysis_home_out_files_pred, "/method", method, sep = "")
  setwd(analysis_home_out_files_method)
  
  
  unl.pred.file = paste(paste(analysis_home_out_files_method, "/uORF_function_unlabelled_predictions_sorted_", sep = ""), prior, method, sep='_')
  pos.pred.file = paste(paste(analysis_home_out_files_method, "/uORF_function_positive_predictions_sorted_", sep = ""), prior, method, sep='_')
  
  
  
  unl.pred = read.table(unl.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the unlabelled uORFs
  pos.pred = read.table(pos.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the positive uORFs
  
  unl.pred$V2 = check_inf(unl.pred$V2)
  pos.pred$V2 = check_inf(pos.pred$V2)
  
  
  
  unl.pred.samp = unl.pred[sample(1:length(as.matrix(unl.pred[,1])),10000),]

  
  # bind the set of predictions for the positive uORFs with the sample of unlabelled ones selected
  target_pred = rbind(pos.pred,unl.pred.samp)
  
  # store the number of columns in the predictions of positive (2)
  ncols = ncol(pos.pred)
  
  # create a matrix of 1s for the positive set of uORFs and of 0s for the sample of unlabelled uORFs
  class.pos <- matrix(sample(1, (ncol(pos.pred)*nrow(pos.pred)), replace=T), ncol=ncols)
  class.unl <- matrix(sample(0, (ncol(unl.pred.samp)*nrow(unl.pred.samp)), replace=T), ncol=ncols)
  
  # we bind the two matrices, and this are the expected labels of the 
  target_class <-rbind(class.pos,class.unl)
  
  
  # this function prepares the data for the performance evaluation using ROCR library
  pred <- prediction(target_pred, target_class)
  
  
  # compute tpr and fpr performance measures
  perf <- performance(pred,"tpr","fpr")
  
  
  # pred_comb <- prediction(target_pred_comb, target_class_comb)
  # perf_comb <- performance(pred,"tpr","fpr")
  
  
  #This ROC curve only counts the positive data (excludes neutral from the analysis).
  roc_name = paste("./ROC_prime_method", method, ".svg", sep = "")
  svg(roc_name)
  par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
  plot(perf,col="black",lty=3, lwd=3)
  auc <- performance(pred,"auc")
  auc <- unlist(slot(auc, "y.values"))
  auc<-round(auc, digits = 2)
  auct <- paste(c("(AUC)  = "),auc,sep="")
  legend(0.3,0.6,c(auct,"\n"),border="white",cex=1.7,box.col = "white")
  dev.off()
  
  print(pos.pred[1,2])
  print(unl.pred[1,2])
  
  
  bench_hist = hist(unl.pred[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for unlabelled uORFs", col="lightblue")
  
  svg('rplot.svg')
  par(mfrow = c(2,1))
  hist(unl.pred[,2], breaks=bench_hist$breaks, main="Class probability ratio distribution for unlabelled uORFs", col="lightblue")
  hist(pos.pred[,2], breaks=bench_hist$breaks, main="Class probability ratio distribution for positive uORFs")
  #hist(neu.pred[,2], xlim[1], xlim[2], breaks=bench_hist$breaks, main="Class probability ratio distribution for neutral uORFs", col="lightblue")
  dev.off()
  
  # svg('rplot_neu.svg')
  # par(mfrow = c(1,1))
  # # xlim = c(1.00, 1000)
  # xlim = c(-100, 100)
  # dev.off()
  
  svg('rplot2.svg')
  par(mfrow = c(2,1))
  # xlim = c(1.00, 1000)
  xlim = c(0, max(c(unl.pred[,2], pos.pred[,2])))
  hist(unl.pred[,2], xlim = xlim, breaks=30, main="Class probability ratio distribution for unlabelled uORFs\n predicted as positives", col="lightblue")
  hist(pos.pred[,2], xlim = xlim, breaks=30, main="Class probability ratio distribution for positive uORFs\n predicted as positives")
  dev.off()
  
  
  # svg('rplot2_neu.svg')
  # par(mfrow = c(1,1))
  # # xlim = c(1.00, 1000)
  # xlim = c(0, 100)
  # hist(neu.pred[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for neutral uORFs\n predicted as positives", col="lightblue")
  # dev.off()
  
  high.thresh = quantile(unl.pred[,2], probs = c(0.93))
  high.preds = unl.pred[unl.pred[,2] >= high.thresh,]
  
  write.table(high.preds, file="top7%_predicted_functional_uORFs.txt", row.names=FALSE, quote=FALSE, sep="\t", col.names=FALSE)
  
  print(quantile(pos.pred[pos.pred[,2] > 1.00, 2], probs = c(0.45)))
  print(quantile(unl.pred[unl.pred[,2] > 1.00, 2], probs = c(0.45)))
  summary(pos.pred[pos.pred[,2] > 1.00, 2])
  
  summary(unl.pred[unl.pred[,2] > 1.00, 2])
  

  # Discretized results with a 0.61 prior.
  
  # setwd("~/Sandbox/")
  
  # disc.file = "uORF_function_unlabelled_predictions_sorted_discretized_2014-09-23"
  
  # disc.pred = read.table(disc.file, header=FALSE,sep="\t",comment.char="", quote="")
  
  # disc.pos.file = "uORF_function_positive_predictions_sorted_discretized_2014-09-23"
  
  # disc.pos = read.table(disc.pos.file, header=FALSE,sep="\t",comment.char="", quote="")
  
  d.unl = get.pos(unl.pred)
  d.pos = get.pos(pos.pred)
  
  svg('rplot3.svg')
  boxplot(list(d.unl[,2], d.pos[,2]) , outline=FALSE, boxwex=0.3, names=c("unlabelled", "Positive"), main="Class probability ratio quartiles\n for uORFs\n predicted as positives")
  # boxplot(d.pos[,2], outline=FALSE, add=TRUE)
  dev.off()
  
  svg('rplot4.svg')
  boxplot(list(unl.pred[,2], pos.pred[,2]) , outline=FALSE, boxwex=0.3, names=c("unlabelled", "Positive"), main="Class probability ratio quartiles\n for uORFs")
  # boxplot(d.pos[,2], outline=FALSE, add=TRUE)
  dev.off()
  
  #comb.pred = rbind(unl.pred,pos.pred,neu.pred)
  comb.pred = rbind(unl.pred,pos.pred)
  comb.pred.finite = (comb.pred[is.finite(comb.pred[,2]),])
  
  # sampling_number = min(c(length(neu.pred[,2]), length(pos.pred[,2])))
  sampling_number = min(c(length(pos.pred[,2])))
  
  pos.dist.funct = ecdf(sample(pos.pred[,2], sampling_number, replace = FALSE))
  unl.dist.funct = ecdf(sample(unl.pred[,2], sampling_number, replace = FALSE))
  # neu.dist.funct = ecdf(sample(neu.pred[,2], sampling_number, replace = FALSE))
  cumul.dist.funct = ecdf(sample(comb.pred[,2], sampling_number, replace = FALSE))
  
  svg('pos_unl_combined_CDF.svg')
  # par(mfrow = c(3,1))
  par(mfrow = c(2,1))
  plot(unl.dist.funct, xlim=c(-30, 30), ylim=c(0, 1), col='green', pch="", lwd=10)
  plot(pos.dist.funct, xlim=c(-30, 30), ylim=c(0, 1), col='red', pch="", lwd=10)
  # plot(neu.dist.funct, xlim=c(-30, 30), ylim=c(0, 1), col='blue', pch="", lwd=10)
  
  dev.off()
  
  xandytest = pos.dist.funct
  
  # plotting the cumulative distribution function, and first derivative of the CDF (based on spline fit) for the positive CDF
  
  xandytest = as.list(environment(pos.dist.funct))
  
  x = xandytest[[2]]
  y = xandytest[[3]]
  
  print(x)
  print(y)
  print(length(x))
  print(length(y))
  
  
  
  svg('pos_CDF.svg')
  spl = smooth.spline(x, y, df = 10)
  
  pred = predict(spl, deriv = 1)
  normpredy = pred$y / max(pred$y)
  
  pred1 = predict(spl, deriv = 2)
  normpred1y = pred1$y / abs(min(pred1$y))
  
  plot (x, y, ylim=c(-1, 1))
  lines(spl, col=1)
  lines(pred$x, normpredy, col=2)
  lines(pred1$x, normpred1y, col=3)
  
  legend(-20,-0.5, c("CDF","deriv","deriv2"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5),col=c(1,2,3))
  
  dev.off()
  
  
  # plotting the cumulative distribution function, and first derivative of the CDF (based on spline fit) for the positive CDF
  
  xandytest = as.list(environment(cumul.dist.funct))
  
  x = xandytest[[2]]
  y = xandytest[[3]]
  
  print(x)
  print(y)
  print(length(x))
  print(length(y))
  
  svg('unl_CDF.svg')
  spl = smooth.spline(x, y, df = 10)
  
  pred = predict(spl, deriv = 1)
  normpredy = pred$y / max(pred$y)
  
  pred1 = predict(spl, deriv = 2)
  normpred1y = pred1$y / max(pred1$y)
  
  plot (x, y, ylim=c(-1, 1))
  lines(spl, col=1)
  lines(pred$x, normpredy, col=2)
  lines(pred1$x, normpred1y, col=3)
  
  legend(10,-0.5, c("CDF","deriv","deriv2"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5),col=c(1,2,3))
  
  dev.off()
  
  
  print(xandytest[3])
  print(length(comb.pred[,2]))
  x = xandytest[[2]]
  y = xandytest[[3]]
  xandytest[3]
  print(length(xandytest))
  
}

unl.pred = read.table(unl.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the unlabelled uORFs
pos.pred = read.table(pos.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the positive uORFs

unl.pred$V2 = check_inf(unl.pred$V2)
pos.pred$V2 = check_inf(pos.pred$V2)

# Save this plots in the plots folder
setwd(analysis_home)
setwd("./plots/")

# sample 10000 uORFs from the whole set of unlabelled uORFs predictions
unl.pred.samp = unl.pred[sample(1:length(as.matrix(unl.pred[,1])),10000),]
# nunl.pred.samp = unl.pred.samp[unl.pred.samp$V2 < 0 ,] 
# punl.pred.samp = unl.pred.samp[unl.pred.samp$V2 >= 0 ,] 
# unl.pred.samp = rbind(nunl.pred.samp, head(punl.pred.samp))

library("ROCR")
# bind the set of predictions for the positive uORFs with the sample of unlabelled ones selected
target_pred = rbind(pos.pred,unl.pred.samp)

# store the number of columns in the predictions of positive (2)
ncols = ncol(pos.pred)

# create a matrix of 1s for the positive set of uORFs and of 0s for the sample of unlabelled uORFs
class.pos <- matrix(sample(1, (ncol(pos.pred)*nrow(pos.pred)), replace=T), ncol=ncols)
class.unl <- matrix(sample(0, (ncol(unl.pred.samp)*nrow(unl.pred.samp)), replace=T), ncol=ncols)

# we bind the two matrices, and this are the expected labels of the 
target_class <-rbind(class.pos,class.unl)


# this function prepares the data for the performance evaluation using ROCR library
pred <- prediction(target_pred, target_class)


# compute tpr and fpr performance measures
perf <- performance(pred,"tpr","fpr")


# pred_comb <- prediction(target_pred_comb, target_class_comb)
# perf_comb <- performance(pred,"tpr","fpr")

require(ggplot2)

#This ROC curve only counts the positive data (excludes neutral from the analysis).
svg("./ROC_prime.svg")
par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
plot(perf,col="black",lty=3, lwd=3)
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
auc<-round(auc, digits = 2)
auct <- paste(c("(AUC)  = "),auc,sep="")
legend(0.3,0.6,c(auct,"\n"),border="white",cex=1.7,box.col = "white")
dev.off()

print(pos.pred[1,2])
print(unl.pred[1,2])


bench_hist = hist(unl.pred[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for unlabelled uORFs", col="lightblue")

svg('rplot.svg')
par(mfrow = c(2,1))
hist(unl.pred[,2], breaks=bench_hist$breaks, main="Class probability ratio distribution for unlabelled uORFs", col="lightblue")
hist(pos.pred[,2], breaks=bench_hist$breaks, main="Class probability ratio distribution for positive uORFs")
#hist(neu.pred[,2], xlim[1], xlim[2], breaks=bench_hist$breaks, main="Class probability ratio distribution for neutral uORFs", col="lightblue")
dev.off()

# svg('rplot_neu.svg')
# par(mfrow = c(1,1))
# # xlim = c(1.00, 1000)
# xlim = c(-100, 100)
# dev.off()

svg('rplot2.svg')
par(mfrow = c(2,1))
# xlim = c(1.00, 1000)
xlim = c(0, max(c(unl.pred[,2], pos.pred[,2])))
hist(unl.pred[,2], xlim = xlim, breaks=30, main="Class probability ratio distribution for unlabelled uORFs\n predicted as positives", col="lightblue")
hist(pos.pred[,2], xlim = xlim, breaks=30, main="Class probability ratio distribution for positive uORFs\n predicted as positives")
dev.off()


# svg('rplot2_neu.svg')
# par(mfrow = c(1,1))
# # xlim = c(1.00, 1000)
# xlim = c(0, 100)
# hist(neu.pred[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for neutral uORFs\n predicted as positives", col="lightblue")
# dev.off()

high.thresh = quantile(unl.pred[,2], probs = c(0.93))
high.preds = unl.pred[unl.pred[,2] >= high.thresh,]

write.table(high.preds, file="top7%_predicted_functional_uORFs.txt", row.names=FALSE, quote=FALSE, sep="\t", col.names=FALSE)

print(quantile(pos.pred[pos.pred[,2] > 1.00, 2], probs = c(0.45)))
print(quantile(unl.pred[unl.pred[,2] > 1.00, 2], probs = c(0.45)))
summary(pos.pred[pos.pred[,2] > 1.00, 2])

summary(unl.pred[unl.pred[,2] > 1.00, 2])

#------------------------------------------

# Distribution plots and Box and Whisker plots to evaluate positive predictions

hist.ratio = function(pred, xmin, xmax, breaks=100, label="unlabelled", col='lightblue', prior=0.61) {
  hist(pred[pred > xmin & pred < xmax], xlim=c(xmin, xmax), col=col, breaks=breaks, main=paste("Class probability ratio distribution\n for", label, "uORFs\n predicted as positives"), xlab=paste("Positive/unlabelled probability ratio\n with", prior, "prior distribution"))
  
}

get.pos = function(data) {
  return(data[data[,2] > 1 & !is.na(data[,2]),])
}

# Discretized results with a 0.61 prior.

# setwd("~/Sandbox/")

# disc.file = "uORF_function_unlabelled_predictions_sorted_discretized_2014-09-23"

# disc.pred = read.table(disc.file, header=FALSE,sep="\t",comment.char="", quote="")

# disc.pos.file = "uORF_function_positive_predictions_sorted_discretized_2014-09-23"

# disc.pos = read.table(disc.pos.file, header=FALSE,sep="\t",comment.char="", quote="")

d.unl = get.pos(unl.pred)
d.pos = get.pos(pos.pred)

svg('rplot3.svg')
boxplot(list(d.unl[,2], d.pos[,2]) , outline=FALSE, boxwex=0.3, names=c("unlabelled", "Positive"), main="Class probability ratio quartiles\n for uORFs\n predicted as positives")
# boxplot(d.pos[,2], outline=FALSE, add=TRUE)
dev.off()

svg('rplot4.svg')
boxplot(list(unl.pred[,2], pos.pred[,2]) , outline=FALSE, boxwex=0.3, names=c("unlabelled", "Positive"), main="Class probability ratio quartiles\n for uORFs")
# boxplot(d.pos[,2], outline=FALSE, add=TRUE)
dev.off()

#comb.pred = rbind(unl.pred,pos.pred,neu.pred)
comb.pred = rbind(unl.pred,pos.pred)
comb.pred.finite = (comb.pred[is.finite(comb.pred[,2]),])

# sampling_number = min(c(length(neu.pred[,2]), length(pos.pred[,2])))
sampling_number = min(c(length(pos.pred[,2])))

pos.dist.funct = ecdf(sample(pos.pred[,2], sampling_number, replace = FALSE))
unl.dist.funct = ecdf(sample(unl.pred[,2], sampling_number, replace = FALSE))
# neu.dist.funct = ecdf(sample(neu.pred[,2], sampling_number, replace = FALSE))
cumul.dist.funct = ecdf(sample(comb.pred[,2], sampling_number, replace = FALSE))

svg('pos_unl_combined_CDF.svg')
# par(mfrow = c(3,1))
par(mfrow = c(2,1))
plot(unl.dist.funct, xlim=c(-30, 30), ylim=c(0, 1), col='green', pch="", lwd=10)
plot(pos.dist.funct, xlim=c(-30, 30), ylim=c(0, 1), col='red', pch="", lwd=10)
# plot(neu.dist.funct, xlim=c(-30, 30), ylim=c(0, 1), col='blue', pch="", lwd=10)

dev.off()

xandytest = pos.dist.funct

# plotting the cumulative distribution function, and first derivative of the CDF (based on spline fit) for the positive CDF

xandytest = as.list(environment(pos.dist.funct))

x = xandytest[[2]]
y = xandytest[[3]]

print(x)
print(y)
print(length(x))
print(length(y))



svg('pos_CDF.svg')
spl = smooth.spline(x, y, df = 10)

pred = predict(spl, deriv = 1)
normpredy = pred$y / max(pred$y)

pred1 = predict(spl, deriv = 2)
normpred1y = pred1$y / abs(min(pred1$y))

plot (x, y, ylim=c(-1, 1))
lines(spl, col=1)
lines(pred$x, normpredy, col=2)
lines(pred1$x, normpred1y, col=3)

legend(-20,-0.5, c("CDF","deriv","deriv2"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5),col=c(1,2,3))

dev.off()


# plotting the cumulative distribution function, and first derivative of the CDF (based on spline fit) for the positive CDF

xandytest = as.list(environment(cumul.dist.funct))

x = xandytest[[2]]
y = xandytest[[3]]

print(x)
print(y)
print(length(x))
print(length(y))

svg('unl_CDF.svg')
spl = smooth.spline(x, y, df = 10)

pred = predict(spl, deriv = 1)
normpredy = pred$y / max(pred$y)

pred1 = predict(spl, deriv = 2)
normpred1y = pred1$y / max(pred1$y)

plot (x, y, ylim=c(-1, 1))
lines(spl, col=1)
lines(pred$x, normpredy, col=2)
lines(pred1$x, normpred1y, col=3)

legend(10,-0.5, c("CDF","deriv","deriv2"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5),col=c(1,2,3))

dev.off()


print(xandytest[3])
print(length(comb.pred[,2]))
x = xandytest[[2]]
y = xandytest[[3]]
xandytest[3]
print(length(xandytest))


# This does not work
# svg('rplot7.svg')
# ycs.prime = diff(xandytest)/diff(comb.pred[,2])
# pred.prime = predict(spl, deriv=1)
# 
# plot(ycs.prime)
# lines(pred.prime$y, col=2)
# 
# 
# dev.off()

# svg('rplot6.svg')
# cumul.dist.funct.deriv = predict(c(comb.pred[,2],xandytest), cumul.dist.funct[,1], 1)
# plot(c(comb.pred[,2],xandytest), xlim=c(-30, 30), main="Cumulative Distribution Derivative Function for all uORFs")
# dev.off()









# CHECK OVERLAPPING WITH REFERENCE SETS #
check_reference <- function(x){
  # this function checks if a value is in the set of reference uORFs
  return (x %in% positive_location)
}


uORFs_candidates <- paste(input_home, "allENSG_15-5-2020_at_17-12_complete_uORFs_with_exp_homATG.tsv", sep = "")

uorfs <- read.delim(uORFs_candidates)



# store the performance measures at the end of each iteration and model
performance_measures_df_reference_sets <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(performance_measures_df_reference_sets) = c("method","reference_set","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                                     "recall", "f1_score_like_estimate")




for (method in 0:4){
  analysis_home_out_files_method = paste(analysis_home_out_files_pred, "/method", method, sep = "")
  setwd(analysis_home_out_files_method)
  
  
  unl.pred.file = paste(paste(analysis_home_out_files_method, "/uORF_function_unlabelled_predictions_sorted_", sep = ""), prior, method, sep='_')
  pos.pred.file = paste(paste(analysis_home_out_files_method, "/uORF_function_positive_predictions_sorted_", sep = ""), prior, method, sep='_')
  
  unl.pred = read.table(unl.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the unlabelled uORFs
  pos.pred = read.table(pos.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the positive uORFs
  
  unl.pred$V2 = check_inf(unl.pred$V2)
  pos.pred$V2 = check_inf(pos.pred$V2)
  
  

  a <- c("heart_science","tools_heart","tools_science")
  b <- c("heart_science_intersection","tools_heart_intersection","tools_science_intersection")
  c <- c("science_non_canonicalORFs_matched_with_candidates.tsv", "uORF_tools_matched_with_candidates.tsv",
         "human_heart_uORFs_matched_with_candidates_proper_format.tsv")
  d <- c(a, b, c)
  
  for (name_ref in d){
    #for (name_ref in c("tools_science")){
    dat <- uorfs
    print(name_ref)
    

    # name_ref <-"heart_science"
    # Read the data from the known positive uORFs
    positive_cases <- paste(ref_home, name_ref, sep = "")
    positive <- read.delim(positive_cases, header = F)
    # Paste the columns to make the intersection with location column
    positive_location <- paste(paste(paste(positive$V1, positive$V2, sep = ":"), positive$V3,
                                     sep = "-"), positive$V4, sep = ":")
    
    dat$label <- sapply(dat$location, check_reference)
    #dat$label <- as.logical(as.integer(dat$label))
    
    
    df_for_perf_measures = rbind(unl.pred, pos.pred)
    df_for_perf_measures = cbind(df_for_perf_measures, dat$label[df_for_perf_measures[,1]])
    df_for_perf_measures[,2] = as.integer(df_for_perf_measures[,2] > 0)
    
    
    perf_measures_fragment <- c(method, name_ref, get.perf.measures(df_for_perf_measures[,2], df_for_perf_measures[,3]))
    names(perf_measures_fragment) = c("method","reference_set","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                      "recall", "f1_score_like_estimate")
    
    colnames(performance_measures_df_reference_sets) = c("method","reference_set","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                                         "recall", "f1_score_like_estimate")
    performance_measures_df_reference_sets <- rbind(performance_measures_df_reference_sets,
                                                    perf_measures_fragment)
    
    colnames(performance_measures_df_reference_sets) = c("method","reference_set","positives_as_positives","unlabelled_as_positives","percentage_positives",
                                                        "recall", "f1_score_like_estimate")
    print("performance measures added")
    
  }
  
  
  print(paste("Method ", method, " DONE!", sep = ""))
  print(cat("\n\n\n\n"))
  
    
}

# as the reference set column is not being created the right way, I will create it manually
performance_measures_df_reference_sets$reference_set <- as.factor(rep(d, 5))
#performance_measures_df_reference_sets$reference_set <- as.factor(performance_measures_df_reference_sets$reference_set)
performance_measures_df_reference_sets$method <- as.factor(performance_measures_df_reference_sets$method)

svg(paste(analysis_home_out_files_pred, '/stats_performance_by_method_by_reference_set.svg', sep = ""), width = 18, height = 16)
par(mfrow = c(1,1))
plot(performance_measures_df_reference_sets, col = performance_measures_df_reference_sets$method, pch = 19)
dev.off()

svg(paste(analysis_home_out_files_pred, '/stats_performance_by_method_by_reference_set_by_reference_set.svg', sep = ""), width = 18, height = 16)
par(mfrow = c(1,1))
plot(performance_measures_df_reference_sets, col = performance_measures_df_reference_sets$reference_set, pch = 19)
dev.off()



# together in the same file
pdf(paste(analysis_home_out_files_pred, '/stats_performance_by_method_by_reference_set.pdf', sep = ""), width = 18, height = 16)
par(mfrow = c(1,1))
plot(performance_measures_df_reference_sets, col = performance_measures_df_reference_sets$method, pch = 19)
plot(performance_measures_df_reference_sets, col = performance_measures_df_reference_sets$reference_set, pch = 19)
dev.off()



# Show recall and percentage of positives mean value of each method
aggregate(performance_measures_df_reference_sets$recall, list(performance_measures_df_reference_sets$method), mean)
aggregate(performance_measures_df_reference_sets$percentage_positives, list(performance_measures_df_reference_sets$method), mean)



##################################################
### SHOW OVERLAP BETWEEN METHODS IN ALL uORFS ####
##################################################

df_uORFs <- data.frame(matrix(ncol = 5, nrow = nrow(mod.disc.data)))
colnames(df_uORFs) = c("method_0","method_1","method_2","method_3","method_4")


for (method in 0:4){
  analysis_home_out_files_method = paste(analysis_home_out_files_pred, "/method", method, sep = "")
  setwd(analysis_home_out_files_method)
  
  
  unl.pred.file = paste(paste(analysis_home_out_files_method, "/uORF_function_unlabelled_predictions_sorted_", sep = ""), prior, method, sep='_')
  pos.pred.file = paste(paste(analysis_home_out_files_method, "/uORF_function_positive_predictions_sorted_", sep = ""), prior, method, sep='_')
  
  unl.pred = read.table(unl.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the unlabelled uORFs
  pos.pred = read.table(pos.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the positive uORFs
  
  unl.pred$V2 = check_inf(unl.pred$V2)
  pos.pred$V2 = check_inf(pos.pred$V2)
  
  
  df_method = rbind(unl.pred, pos.pred)
  df_method[,2] = as.integer(df_method[,2] > 0)
  df_method = df_method[order(df_method[,1]),]
  df_uORFs[,paste("method_", method, sep = "")] = df_method[,2]
  
}

library(limma)
a <- vennCounts(df_uORFs)
a

svg(paste(analysis_home_out_files_pred, '/overlap_between_predictions_diff_methods.svg', sep = ""), width = 18, height = 16)
par(mfrow = c(1,1), cex = 1.5)
vennDiagram(a, main = "Positive uORFs", cex=0.7, show.include = T)
dev.off()









##################################################
############# SAVE uORFs WITH SCORES #############
##################################################

initial_uorfs_data_home <- "/home/fcalvet/Documents/replicate_yale_analysis/data/input"
uORFs_candidates <- read.delim(paste(initial_uorfs_data_home, "/allENSG_15-5-2020_at_17-12_complete_uORFs_with_exp_homATG.tsv", sep = ""))

for (method in 0:4){
  analysis_home_out_files_method = paste(analysis_home_out_files_pred, "/method", method, sep = "")
  setwd(analysis_home_out_files_method)
  
  
  unl.pred.file = paste(paste(analysis_home_out_files_method, "/uORF_function_unlabelled_predictions_sorted_", sep = ""), prior, method, sep='_')
  pos.pred.file = paste(paste(analysis_home_out_files_method, "/uORF_function_positive_predictions_sorted_", sep = ""), prior, method, sep='_')
  
  unl.pred = read.table(unl.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the unlabelled uORFs
  pos.pred = read.table(pos.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the positive uORFs

  # remove inf values
  unl.pred$V2 = check_inf(unl.pred$V2)
  pos.pred$V2 = check_inf(pos.pred$V2)
  
  
  joined.pred <- rbind(unl.pred, pos.pred)
  joined.pred <- joined.pred[order(joined.pred$V2, decreasing = T),]
  
  uORFs_candidates_sorted_with_prob <- uORFs_candidates[c(joined.pred$V1),]
  uORFs_candidates_sorted_with_prob$scores <- joined.pred$V2
  write.table(x = uORFs_candidates_sorted_with_prob,
              file = paste(paste(analysis_home_out_files_method, "/uORFs_scored_", sep = ""), prior, "method", method,".tsv", sep='_'),
              quote = F, sep = "\t", row.names = F, col.names = T)
  
  
}


