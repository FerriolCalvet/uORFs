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


# plotting histogram ...
hist.ratio = function(pred, xmin, xmax, breaks=100) {
  hist(pred[pred > xmin & pred < xmax], breaks=breaks)
}

# plotting histogram ...
hist.trial = function(pred, xmin, xmax, breaks=30, main="No main provided", xlab="value", col="lightgreen") {
  hist(pred[pred > xmin & pred < xmax], breaks=breaks, main=main, xlab=xlab, col=col)
}





#-------------------------------------
# Function Definitions

percent_pos = function(training, pred, class, freq=F) {
  # If both are 1 in the same place, then we have percent of positives
  n_class = sum(training == class)
  n_pos = sum(((training == class) + (pred == 1)) == 2)
  if (freq) return(n_pos / n_class)
  else return(n_pos)
}


# get.perf.measures.yao = function(pred, labels) {
#   
#   p_pos_pos = percent_pos(labels, pred, T)
#   p_unl_pos = percent_pos(labels, pred, F)
#   
#   
#   # Precision calculation
#   pos_pos = percent_pos(labels, pred, T, freq=T)
#   unl_pos = percent_pos(labels, pred, F, freq=T)
#   prob_positive = (pos_pos + unl_pos) / total_length
#   recall = p_pos_pos			# These are the same thing
#   
#   F.stat_estimate = recall * recall
#     
#   return(c(p_pos_pos, p_unl_pos, prec, recall, F.stat))
#   
#   # Precision = true-positives / (true-positives + false-positives)
#   # recall = true-positives / (true-positives + false-negatives)
#   
#   # F statistic is the geometric mean of Precision and Recall:
#   # F = 2 * (precision * recall) / (precision + recall)
# }


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
  print(paste("Finished one prior: ", prior))
  return(get.perf.measures(pred, labels))
}



evaluate.diff.priors <- function(data_with_labels, unlabelled_prob, positive_prob){
  scores_diff_prior = list()
  for (pri in c(seq(0.05, 0.95, 0.1))){
    scores_diff_prior = rbind(scores_diff_prior, c(pri,n.Bayes(data_with_labels, u.prb = unlabelled_prob, p.prb = positive_prob, prior = pri)))
  }
  return(scores_diff_prior)
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



sample.mat = function(M, p) {
  n = dim(M)[1]
  m = dim(M)[2]
  f = 100000
  p.data = as.integer(p * f)
  log = c(rep(T, p.data), rep(F, f - p.data))
  log.s = sample(log, n, replace=TRUE)
  return(M[log.s,])
}

write.predictions = function(data, u.prb, p.prb, prior, filename="uORF_function_unlabeled_predictions_sorted_") {
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

#The maximum number of values for a variable is measured.
max.len = 1
for (k in 1:d.N) {
  len = length(table(disc.data[,k]))
  if (len > max.len) max.len <- len
}
print(max.len)


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



# compare.prob.matrix <- function(u.prb_compare, p.prb_compare){
#   u.N = dim(u.prb_compare)[1] # number of variables
#   u.M = dim(u.prb_compare)[2] # number of possible values
#   
#   total_diff <- 0
#   
#   # compare the values for each cell in both matrices
#   for (k in 1:u.N) {
#     for (j in 1:max.len) {
#       total_diff <- total_diff + abs(u.prb_compare[k, j] - p.prb_compare[k, j])
#     }
#   }
#   return(total_diff)
# }


check_inf <- function(x){
  max_not_inf <- max(x[x != Inf])
  x[x == Inf] = max_not_inf

  min_not_inf <- min(x[x != -Inf])
  x[x == -Inf] = min_not_inf

  return(x)
}

write.predictions_unlabelled = function(data, unlabelled, positive, max.len, prior, it_max,
                                        quantile_threshold = 0.97,
                                        filename="uORF_function_unlabeled_predictions_sorted_") {
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
  
  # data = unl.looping.disc
  # unlabelled = unl.looping.disc
  # positive = pos.looping.disc
  # max.len = 4
  
  
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
    
    print("uORFs sorted")
    
    positive_scored <- b[b[,1] > 0,]
    print(paste("Number of positives", nrow(positive_scored)))
    # print(summary(positive_scored))
    
    threshold_positive <- quantile(positive_scored[,1], quantile_threshold) # I should adjust this threshold
    print(paste("Threshold quantile", threshold_positive))
    
    
    selected_from_unlabelled <- rownames(sorted.scores)[sorted.scores[,1] > threshold_positive] # get the positions of the uORFs above the threshold
    
    names_unlabelled_chosen_positives <- c(names_unlabelled_chosen_positives,
                                           selected_from_unlabelled)

    print(paste("chosen_positive", length(names_unlabelled_chosen_positives)))
    names_unlabelled <- setdiff(names_unlabelled, names_unlabelled_chosen_positives)
    
    print(paste("Unlabelled", length(names_unlabelled)))
    
    it_unlabelled <- it_unlabelled[names_unlabelled,]
    print("Unlabelled filtered")
    
    rows_from_unlabelled_to_positive <- data[as.integer(selected_from_unlabelled),]
    
    it_positive <- rbind(it_positive, rows_from_unlabelled_to_positive)
    print("Positive updated")
    
    
    new_difference <- sum(abs(u.prb - p.prb))
    print(paste("The new difference is", new_difference))
    
    
    mean_prob_sum_int <- mean(abs(check_inf(b[,1])))
    print(paste("Mean probability of all remaining unlabelled (absolute value)", mean_prob_sum_int))
    
    
    if (i > 1){
      if (new_difference < diff_between_matrices + 0.1){
        break
      } else {
        diff_between_matrices <- new_difference
      }
      
      if (threshold_positive < 1){
        break
      }
    }
    
    
    
  }
  
  
  print("Unlabelled 'purification' finished")
  
  # Doing the final iteration running with the improved matrices
  # and all the data.
  prob.ratios = n.Bayes.score(data, u.prb, p.prb, prior)
  
  print("Scores computed")
  
  # store the scores into a dataframe and sort them in the proper order
  uORFs.ratio = data.frame(prob.ratios)
  rownames(uORFs.ratio) <- rownames(data) # these are the rownames from the initial data matrix
  #                                         where labelled and unlabelled were together
  uORFs.ratio[,2] <- prob.ratios				# To keep rownames when I sort the data frame
  
  # sort the uORFs according to the scores, and put them into a dataframe
  b = uORFs.ratio[order(-uORFs.ratio[,1]),]
  sorted.scores = data.frame(b[,1])
  rownames(sorted.scores) <- rownames(b)
  
  
  mean_prob_sum <- mean(abs(check_inf(b[,1])))
  print(paste("Final mean probability (absolute value)", mean_prob_sum))
  
  
  
  # write the predictions into the desired file
  output = paste(filename, sep="")
  write.table(sorted.scores, file=output, quote=FALSE, sep="\t", col.names=FALSE)
  
  
  cat("\n\n\n\n")
  
  # return the filename where we saved the data
  return(list(output, u.prb, p.prb))
}


























train.matrices.updating.unlabelled = function(unlabelled, positive, max.len, prior, it_max,
                                        quantile_threshold = 0.97,
                                        filename="probability_matrix_updating_unlabelled_") {
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
  
  # data = unl.looping.disc
  # unlabelled = unl.looping.disc
  # positive = pos.looping.disc
  # max.len = 4
  
  
  names_unlabelled <- rownames(unlabelled)
  names_unlabelled_chosen_positives <- c()
  it_unlabelled <- unlabelled
  
  
  diff_between_matrices <- 0
  
  
  p.prb = compute.single.prob.matrix(positive, max.len)
  
  #for (i in 1:10){
  for (i in 1:it_max){
    cat(paste("\n\nStarting iteration", i))
    print("")
    prob_matrices = compute.single.prob.matrix(it_unlabelled, max.len = max.len)
    #prob_matrices = compute.prob.matrix(it_unlabelled, positive, max.len = max.len)
    
    u.prb = prob_matrices[[1]]
    
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
    
    it_unlabelled <- it_unlabelled[names_unlabelled,]
    print("Unlabelled filtered")
    

    new_difference <- sum(abs(u.prb - p.prb))
    print(paste("The new difference is", new_difference))
    
    
    mean_prob_sum_int <- mean(abs(check_inf(b[,1])))
    print(paste("Mean probability of all remaining unlabelled (absolute value)", mean_prob_sum_int))
    
    
    if (i > 1){
      if (new_difference < diff_between_matrices + 0.1){
        break
      } else {
        diff_between_matrices <- new_difference
      }
      
      if (threshold_positive < 1){
        break
      }
    }
  
    
  }
  
  # return the two probability matrices
  return(list(u.prb, p.prb))
}













train.matrices.updating.unlabelled.positive = function(unlabelled, positive, max.len, prior, it_max,
                                                       quantile_threshold = 0.97,
                                                       filename="probability_matrix_updating_unlabelled_positive_") {
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
  
  # data = unl.looping.disc
  # unlabelled = unl.looping.disc
  # positive = pos.looping.disc
  # max.len = 4
  
  
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
    
    it_unlabelled <- it_unlabelled[names_unlabelled,]
    print("Unlabelled filtered")
    
    rows_from_unlabelled_to_positive <- unlabelled[as.integer(selected_from_unlabelled),]
    
    it_positive <- rbind(it_positive, rows_from_unlabelled_to_positive)
    print("Positive updated")
    
    
    new_difference <- sum(abs(u.prb - p.prb))
    print(paste("The new difference is", new_difference))
    
    
    mean_prob_sum_int <- mean(abs(check_inf(b[,1])))
    print(paste("Mean probability of all remaining unlabelled (absolute value)", mean_prob_sum_int))
    
    
    if (i > 1){
      if (new_difference < diff_between_matrices + 0.1){
        break
      } else {
        diff_between_matrices <- new_difference
      }
      
      if (threshold_positive < 1){
        break
      }
    }

  }
  
  # return the two probability matrices
  return(list(u.prb, p.prb))
}





















#the positive and unlabelled data sets are separated based on the class column, the class column is removed.

all_pos.disc = (disc.data[disc.data[,d.N] == 1,])[, -d.N]	# shaves off the class column
all_unl.disc = (disc.data[disc.data[,d.N] == 0,])[, -d.N]	# shaves off the class column



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
  len = length(table(disc.data[,k]))
  if (len > max.len) max.len <- len
}
print(max.len)

#for (prior in c(seq(0.05, 0.95, 0.05))){
#for (prior in c(seq(0.25, 0.95, 0.1))){
prior = 0.5
for (i in 1:8) {
  
  unl.looping.disc = bound[[i]]
  pos.looping.disc = positives[[i]]
  ret.looping.disc = retained[[i]]
  
  #The number of variables is measured.
  d.N = dim(disc.data)[2]
  max.len = 1
  for (k in 1:d.N) {
    len = length(table(disc.data[,k]))
    if (len > max.len) max.len <- len
  }
  # print(max.len)
  
  
  prob_matrices_list0 = compute.prob.matrix(unl.looping.disc, pos.looping.disc, max.len)
  print("Matrices 0 computed!")
  
  quantile_threshold1 = 0.95
  update_positive_score_threshold1 = 10
  prob_matrices_list1 = train.matrices.updating.unlabelled(unl.looping.disc, pos.looping.disc,
                                                           max.len, prior,
                                                           it_min = 4,
                                                           it_max = 50,
                                                           quantile_threshold = quantile_threshold1,
                                                           update_positive_score_threshold = update_positive_score_threshold1,
                                                           filename=paste(paste(analysis_home_prob_matrices, "/uORF_probability_matrices_updating_unlabelled", sep = ""), i, prior, "1", quantile_threshold1, update_positive_score_threshold1, sep='_'))
  print("Matrices 1 computed!")
  
  
  quantile_threshold2 = 0.99
  update_positive_score_threshold2 = 10
  
  prob_matrices_list2 = train.matrices.updating.unlabelled.positive(unl.looping.disc, pos.looping.disc,
                                                                    max.len, prior,
                                                                    it_min = 5,
                                                                    it_max = 50,
                                                                    quantile_threshold = quantile_threshold2,
                                                                    update_positive_score_threshold = update_positive_score_threshold2,
                                                                    skip_differences = T,
                                                                    filename=paste(paste(analysis_home_prob_matrices, "/uORF_probability_matrices_updating_unlabelled_and_positive", sep = ""), i, prior, "2", quantile_threshold2, update_positive_score_threshold2, sep='_'))
  print("Matrices 2 computed!")
  
  
  u.prb0 = prob_matrices_list0[[1]]
  p.prb0 = prob_matrices_list0[[2]]
  
  u.prb1 = prob_matrices_list1[[1]]
  p.prb1 = prob_matrices_list1[[2]]
  
  u.prb2 = prob_matrices_list2[[1]]
  p.prb2 = prob_matrices_list2[[2]]
  
  
  #write.predictions
  unl.pred.file0 = write.predictions(unl.looping.disc, u.prb0, p.prb0, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, "0", sep='_'))
  pos.pred.file0 = write.predictions(pos.looping.disc, u.prb0, p.prb0, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior,"0", sep='_'))
  ret.pred.file0 = write.predictions(ret.looping.disc, u.prb0, p.prb0, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior,"0", sep='_'))
  
  
  unl.pred.file1 = write.predictions(unl.looping.disc, u.prb1, p.prb1, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, "1", sep='_'))
  pos.pred.file1 = write.predictions(pos.looping.disc, u.prb1, p.prb1, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior,"1", sep='_'))
  ret.pred.file1 = write.predictions(ret.looping.disc, u.prb1, p.prb1, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior,"1", sep='_'))
  
  unl.pred.file2 = write.predictions(unl.looping.disc, u.prb2, p.prb2, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_unlabelled_predictions_sorted_discretized_looping", sep = ""), i, prior, "2", sep='_'))
  pos.pred.file2 = write.predictions(pos.looping.disc, u.prb2, p.prb2, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_positive_predictions_sorted_discretized_looping", sep = ""), i, prior,"2", sep='_'))
  ret.pred.file2 = write.predictions(ret.looping.disc, u.prb2, p.prb2, prior,
                                     filename=paste(paste(analysis_home_out_files_test, "/uORF_function_retained_predictions_sorted_discretized_looping", sep = ""), i, prior,"2", sep='_'))    
  
  print(paste(i, "DONE!"))
  print(cat("\n\n\n\n"))
}
#}

# #Concatenate the files together:

# setwd("~/Sandbox/")


hist.looping = function(pred, xmin, xmax, breaks=100, main="No main provided", xlab="Positive/unlabelled probability ratio\n prior positive distribution is 0.5", col="lightgreen") {
  hist(pred[pred > xmin & pred < xmax], breaks=breaks, main=main, xlab=xlab, col=col)
}

for (method in c(0,1,2)){
  
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
  
  
  svg(paste(analysis_home_plots, '/ret_pos_unl_distro_method_second_round_', method,'.svg', sep = ""))
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
  
  
}
