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


setwd("/home/fcalvet/Documents/replicate_yale_analysis/scripts")
setwd("..")
home <- getwd()

# day_date = Sys.Date()
day_date = "2020-05-18"

analysis_home <- paste(home,day_date, sep = "/analysis_")
system(paste("mkdir -p", analysis_home))


ref_home <- paste(home, "/data/reference_sets/", sep = "")
input_home <- paste(home, "/data/input/", sep = "")

analysis_home_data <- paste(analysis_home,"data", sep = "/")
system(paste("mkdir -p", analysis_home_data))

# name_ref <- "heart_science"
name_ref <- "heart_science_intersection"
disc.data <- read.delim(paste(analysis_home_data, "/", name_ref, "_discretized_uORF_table_no_row_names.tsv", sep = ""))

analysis_home_plots <- paste(analysis_home,"plots", sep = "/")
analysis_home_plots_histograms <- paste(analysis_home_plots,"var_histogram", sep = "/")
system(paste("mkdir -p", analysis_home_plots))
system(paste("mkdir -p", analysis_home_plots_histograms))


analysis_home_out_files <- paste(analysis_home,"output", sep = "/")
system(paste("mkdir -p", analysis_home_out_files))

setwd(analysis_home)








#############################################

#-------------------------------------
# Function Definitions

percent_pos = function(training, pred, class, freq=F) {
  # If both are 1 in the same place, then we have percent of positives
  n_class = sum(training == class)
  n_pos = sum(((training == class) + (pred == 1)) == 2)
  if (freq) return(n_pos)
  else return(n_pos / n_class)
}


get.perf.measures.yao = function(pred, labels) {
  p_pos_pos = percent_pos(labels, pred, T)
  p_unl_pos = percent_pos(labels, pred, F)
  
  # Precision calculation
  pos_pos = percent_pos(labels, pred, T, freq=T)
  unl_pos = percent_pos(labels, pred, F, freq=T)
  prec = pos_pos / (pos_pos + unl_pos)
  recall = p_pos_pos			# These are the same thing
  F.stat = 2 * (prec*recall) / (prec + recall)
  return(c(p_pos_pos, p_unl_pos, prec, recall, F.stat))
  
  # Precision = true-positives / (true-positives + false-positives)
  # recall = true-positives / (true-positives + false-negatives)
  
  # F statistic is the geometric mean of Precision and Recall:
  # F = 2 * (precision * recall) / (precision + recall)
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
  return(get.perf.measures.yao(pred, labels))
}

n.Bayes.score = function(data, u.prb, p.prb, prior) {
  scores = data[,1]			# To retain the row names
  M = dim(data)[1]
  for (i in 1:M) {
    scores[i] = pred.Bayes(as.numeric(data[i,]), u.prb, p.prb, prior, ratio=T)
  }
  return(scores)
}

# p is the probability of positive
pred.Bayes = function(x, u.prb, p.prb, p, ratio=FALSE) {
  p.unl = log10((1-p))
  n = dim(u.prb)[1]
  for (j in 1:n) {
    p.unl = p.unl + log10(u.prb[j, x[j]])
    # deb = p.prb[j, x[j]]
    # print(paste("Deb:", deb, "j:", j, "x[j]", x[j]))
    # if (x[j] == 1) {
    # p.unl = p.unl * u.prb[j]
    # }
    # else if (x[j] == 2) {
    # p.unl = p.unl * (1 - u.prb[j])
    # }
    # else print("ERROR, feature value not discretized as expected")
    
  }
  # print(p.prb)
  # print(u.prb)
  p.pos = log10(p)
  for (j in 1:n) {
    p.pos = p.pos + log10(p.prb[j, x[j]])
    # deb = p.prb[j, x[j]]
    # print(paste("Deb:", deb, "j:", j, "x[j]", x[j]))
    # if (x[j] == 1) {
    # p.pos = p.pos * p.prb[j]
    # }
    # else if (x[j] == 2) {
    # p.pos = p.pos * (1 - p.prb[j])
    # }
    # else print("ERROR, feature value not discretized as expected")
    
  }
  # print(p.pos)
  # print(p.unl)
  if (ratio == TRUE) {
    return(p.pos - p.unl)
  }
  if (p.pos > p.unl) return(TRUE)
  else return(FALSE)
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
  prob.ratios = n.Bayes.score(data, u.prb, p.prb, prior)
  uORFs.ratio = data.frame(prob.ratios)
  rownames(uORFs.ratio) <- rownames(data)
  uORFs.ratio[,2] <- prob.ratios				# To keep rownames when I sort the data frame
  
  b = uORFs.ratio[order(-uORFs.ratio[,1]),]
  sorted.scores = data.frame(b[,1])
  rownames(sorted.scores) <- rownames(b)
  
  output = paste(filename, sep="")
  write.table(sorted.scores, file=output, quote=FALSE, sep="\t", col.names=FALSE)
  return(output)
}

hist.ratio = function(pred, xmin, xmax, breaks=100) {
  hist(pred[pred > xmin & pred < xmax], breaks=breaks)
  
}


hist.trial = function(pred, xmin, xmax, breaks=30, main="No main provided", xlab="value", col="lightgreen") {
  hist(pred[pred > xmin & pred < xmax], breaks=breaks, main=main, xlab=xlab, col=col)
  
}






# Imaging the table ------------------------------------

#The names of the features before removal of homogenous features are written in ./featuresbefore.txt

allfeaturenamesbefore = colnames(disc.data)
write.table(allfeaturenamesbefore, file=paste("./output/featuresbefore.txt", sep=""), quote=FALSE, sep="\t",
            col.names=FALSE, row.names=FALSE)

#The data from each of the discretized features (both unl and pos combined) is plotted next, on an individual basis.

xlim = c(0, 10000)
xlab = paste("trial", sep="")
tracker = 1

for (aname in allfeaturenamesbefore) {
  
  #samplename = paste("./plots/", allfeaturenames[k], ".pdf", sep="")
  pdf(paste("./plots/var_histogram/", aname, ".pdf", sep=""))
  hist.trial(disc.data[,tracker], xlim[1], xlim[2], breaks=30, main=aname, xlab, col="lightblue")
  dev.off()
  tracker <- tracker + 1
  
}

# Imaging the table-------------------------------------

print("Removing homogenous features...")

# this is not working, but we can skip it because there is variation in all variables
disc.data = screen(disc.data)			# Throw out homogenous features (all 1's or all 2's)

write.table(disc.data, file=paste("./output/discretized_data_all_1_screened.txt", sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)

#The names of each feature, retained following discretization, are written to the file ./featuresafter.txt

allfeaturenamesafter = colnames(disc.data)
write.table(allfeaturenamesafter, file=paste("./output/featuresafter.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

#the names of the rejected features (homogenous features) are written to the file ./rejected.txt

rejected = setdiff(allfeaturenamesbefore, allfeaturenamesafter)
write.table(rejected, file=paste("./output/rejected.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

print("Homogenous features rejected.")

#the dimensions of the discretized data are measured

d.M = dim(disc.data)[1]
d.N = dim(disc.data)[2]

#the positive and unlabeled data sets are separated based on the class column, the class column is removed.

all_pos.disc = (disc.data[disc.data[,d.N] == 1,])[, -d.N]	# shaves off the class column
all_unl.disc = (disc.data[disc.data[,d.N] == 0,])[, -d.N]	# shaves off the class column







#A matrix is constructed, allowing for measurement of the t-test value for each variable.
allfeaturenamesafter <- names(disc.data)

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

#condense the tissue types variables:
kstestmatordered_GTEXcondensed = kstestmatordered
removal_vector = c("Bone.Marrow","Liver","Pituitary","Spleen","Bladder","Skin","Stomach","Lung","Nerve","Small.Intestine","Blood.Vessel","Muscle","Adipose.Tissue","Pancreas","Salivary.Gland","Esophagus","Blood","Brain","Thyroid","Fallopian.Tube","Vagina","Kidney","Prostate","Uterus","Cervix.Uteri","Colon","Breast","Heart","Testis","Ovary","Adrenal.Gland")
forcombined_average = kstestmatordered_GTEXcondensed[rownames(kstestmatordered_GTEXcondensed) %in% removal_vector, ]

combined_average = mean(as.numeric(forcombined_average[,1]))
combined_pvalue = mean(as.numeric(forcombined_average[,2]))
combined_type = forcombined_average[1,3]
combined_test = forcombined_average[1,4]
combined_data = forcombined_average[1,5]

kstestmatordered_GTEXcondensed = kstestmatordered_GTEXcondensed[!rownames(kstestmatordered_GTEXcondensed) %in% removal_vector, ]
kstestmatordered_GTEXcondensed = rbind(kstestmatordered_GTEXcondensed, "GTEX_combined" = c(combined_average, combined_pvalue, combined_type,combined_test, combined_data))
kstestmatordered_GTEXcondensed=kstestmatordered_GTEXcondensed[order(abs(as.numeric(kstestmatordered_GTEXcondensed[,1]))),]

write.table(kstestmatordered_GTEXcondensed, file=paste("./output/kstestmatordered_GTEXcondensed.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)

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
barplot(abs(as.numeric(kstestmatordered_GTEXcondensed[1:10,1])), main="ks-test", horiz=TRUE, names.arg=tail(dimnames(kstestmatordered_GTEXcondensed)[[1]], 10), las=1)
dev.off()




#Legend for the t-test table is constructed, so that the numbers on the table can be associated with variables.
ttestlegend = ttestmatordered[,1:2]
ttestlegend[,2] = ttestmatordered[,1]
ttestlegend[,1] = 1:dim(ttestmatordered)[1]
write.table(ttestlegend, file=paste("./output/ttestlegend.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
write.table(ttestlegend, file=paste("./plots/ttestlegend.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=TRUE)
#---------------------------------------------------




#The dimensions of both the unlabeled, and positive, matrices of discretized values are measured.

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
  len = length(table(disc.data[,k]))
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


#A table with the number of instances of each value of each variable is constructed (both positive and unlabeled)
u.count = matrix(0, nrow=u.N, ncol=max.len)
p.count = matrix(0, nrow=u.N, ncol=max.len)
dimnames(u.count)[1] = dimnames(all_unl.disc)[2]
dimnames(p.count)[1] = dimnames(all_pos.disc)[2]


#A table with the frequency (percentwise) of each value of for each variable is constructed (both positive and unlabeled)
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
  barplot(trialmatnormalized, main=anames[k], xlab="value (L to R, 1 to 2 (or 3))", col=c("darkblue","red"), legend = c("positive", "unlabeled"))
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
  
  barplot(trialmatnormalized, main=anames[k], xlab="value", col=c("darkblue","red"), legend = c("positive", "unlabeled"))
  tracker <- tracker + 1
  
}
dev.off()

# 0.53 was a good prior in the cross validation
# 0.61 is the prior I chose to compare normal and corrected discretized predictions
# Ahh, it turns out that in my corrected discretized predictions the CDS.start.to.uORF.end variable
# has all values for positive as 1, thus giving 0 probability for any unlabeled uORF that doesn't meet
# that criteria. While this may be cause for concern and motivate using one of the Naive Bayes
# modifications to keep any probability from being 0, I think ti is just fine to leave out these
# examples. However, I will get an infinity error for dividing by this 0 in outputing the result

setwd("./output/")

prior = 0.5

#write.predictions

unl.pred.file = write.predictions(all_unl.disc, u.prb, p.prb, prior, filename=paste("uORF_function_unlabeled_predictions_sorted_discretized_", date(), sep=''))
pos.pred.file = write.predictions(all_pos.disc, u.prb, p.prb, prior, filename=paste("uORF_function_positive_predictions_sorted_discretized_", date(), sep=''))

# I used these ones to check the effect of using an improvement of the matrices
unl1.pred.file = write.predictions(all_unl.disc, u.prb, p.prb, prior, filename=paste("1uORF_function_unlabeled_predictions_sorted_discretized_", date(), sep=''))
pos1.pred.file = write.predictions(all_pos.disc, u.prb, p.prb, prior, filename=paste("1uORF_function_positive_predictions_sorted_discretized_", date(), sep=''))


# ret.pred.file = write.predictions(retpos.disc, u.prb, p.prb, prior, filename=paste("uORF_function_retained_predictions_sorted_discretized_", date, sep=''))

# unl.pred.file = paste("uORF_function_unlabeled_predictions_sorted_discretized_", date(), sep='')
# pos.pred.file = paste("uORF_function_positive_predictions_sorted_discretized_", date(), sep='')

# ret.pred.file = paste("uORF_function_retained_predictions_sorted_discretized_", date, sep='')













#-------------------------------------

# Function Definitions 

hist.ratio = function(pred, xmin, xmax, breaks=100, main="No main provided", xlab="Positive/unlabeled probability ratio\n prior positive distribution is 0.61", col="lightgreen") {
  hist(pred[pred > xmin & pred < xmax], breaks=breaks, main=main, xlab=xlab, col=col)
  
}


check_inf <- function(x){
  if (sum(x == Inf)){
    max_not_inf <- max(x[x != Inf])
    x[x == Inf] = max_not_inf
  }
  return(x)
}

#-------------------------------------
# setwd(paste(home, "Sandbox/", sep = ""))

unl.pred = read.table(unl.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the unlabelled uORFs
pos.pred = read.table(pos.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the positive uORFs

# Again, same as before
# unl.pred = read.table(unl1.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the unlabelled uORFs
# pos.pred = read.table(pos1.pred.file, header=FALSE,sep="\t",comment.char="", quote="") # predictions for the positive uORFs


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


pred_comb <- prediction(target_pred_comb, target_class_comb)
perf_comb <- performance(pred,"tpr","fpr")

require(ggplot2)

#This ROC curve only counts the positive data (excludes neutral from the analysis).
svg("./ROC_prime111.svg")
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


bench_hist = hist(unl.pred[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for unlabeled uORFs", col="lightblue")

svg('rplot.svg')
par(mfrow = c(2,1))
hist(unl.pred[,2], breaks=bench_hist$breaks, main="Class probability ratio distribution for unlabeled uORFs", col="lightblue")
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
hist(unl.pred[,2], xlim = xlim, breaks=30, main="Class probability ratio distribution for unlabeled uORFs\n predicted as positives", col="lightblue")
hist(pos.pred[,2], xlim = xlim, breaks=30, main="Class probability ratio distribution for positive uORFs\n predicted as positives")
dev.off()


# svg('rplot2_neu.svg')
# par(mfrow = c(1,1))
# # xlim = c(1.00, 1000)
# xlim = c(0, 100)
# hist(neu.pred[,2], xlim[1], xlim[2], breaks=30, main="Class probability ratio distribution for neutral uORFs\n predicted as positives", col="lightblue")
# dev.off()

setwd(analysis_home_out_files)
high.thresh = quantile(unl.pred[,2], probs = c(0.93))
high.preds = unl.pred[unl.pred[,2] >= high.thresh,]

write.table(high.preds, file="top7%_predicted_functional_uORFs.txt", row.names=FALSE, quote=FALSE, sep="\t", col.names=FALSE)

print(quantile(pos.pred[pos.pred[,2] > 0, 2], probs = c(0.45)))
print(quantile(unl.pred[unl.pred[,2] > 0, 2], probs = c(0.45)))
summary(pos.pred[pos.pred[,2] > 1.00, 2])

summary(unl.pred[unl.pred[,2] > 1.00, 2])

#------------------------------------------

# Distribution plots and Box and Whisker plots to evaluate positive predictions

hist.ratio = function(pred, xmin, xmax, breaks=100, label="unlabeled", col='lightblue', prior=0.61) {
  hist(pred[pred > xmin & pred < xmax], xlim=c(xmin, xmax), col=col, breaks=breaks, main=paste("Class probability ratio distribution\n for", label, "uORFs\n predicted as positives"), xlab=paste("Positive/unlabeled probability ratio\n with", prior, "prior distribution"))
}

get.pos = function(data) {
  positive_threshold <- 0
  return(data[data[,2] > positive_threshold & !is.na(data[,2]),])
}

# Discretized results with a 0.61 prior.

# setwd("~/Sandbox/")

# disc.file = "uORF_function_unlabeled_predictions_sorted_discretized_2014-09-23"

# disc.pred = read.table(disc.file, header=FALSE,sep="\t",comment.char="", quote="")

# disc.pos.file = "uORF_function_positive_predictions_sorted_discretized_2014-09-23"

# disc.pos = read.table(disc.pos.file, header=FALSE,sep="\t",comment.char="", quote="")

d.unl = get.pos(unl.pred)
d.pos = get.pos(pos.pred)

svg('rplot3.svg')
boxplot(list(d.unl[,2], d.pos[,2]) , outline=FALSE, boxwex=0.3, names=c("Unlabeled", "Positive"), main="Class probability ratio quartiles\n for uORFs\n predicted as positives")
# boxplot(d.pos[,2], outline=FALSE, add=TRUE)
dev.off()

svg('rplot4.svg')
boxplot(list(unl.pred[,2], pos.pred[,2]) , outline=FALSE, boxwex=0.3, names=c("Unlabeled", "Positive"), main="Class probability ratio quartiles\n for uORFs")
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

setwd(analysis_home_plots)
svg('pos_unl_combined_CDF.svg')
# par(mfrow = c(3,1))
par(mfrow = c(2,1))
plot(unl.dist.funct, xlim=c(-30, 30), ylim=c(0, 1), col='green', pch="", lwd=10, main = "Unlabelled")
plot(pos.dist.funct, xlim=c(-30, 30), ylim=c(0, 1), col='red', pch="", lwd=10, main = "Positive")
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


##
#
# Get the genes for the top predictions
home <- "/home/ferriol/Desktop/replicate_yale_analysis/"
uORFs_candidates <- read.delim(paste(home, "complete_uORFs_with_exp_homATG.tsv", sep = ""))

top_genes <- uORFs_candidates$gene_id[c(high.preds$V1)]
write.table(x = top_genes, file = paste(home, "genes_of_top_predicted_uORFs.tsv", sep = ""),
            quote = F, sep = "\t", row.names = F, col.names = F)

joined.pred <- rbind(unl.pred, pos.pred)
joined.pred <- joined.pred[order(joined.pred$V2, decreasing = T),]

uORFs_candidates_sorted_with_prob <- uORFs_candidates[c(joined.pred$V1),]
uORFs_candidates_sorted_with_prob$scores <- joined.pred$V2
write.table(x = uORFs_candidates_sorted_with_prob, file = paste(home, "scored_uORFs.tsv", sep = ""),
            quote = F, sep = "\t", row.names = F, col.names = F)
