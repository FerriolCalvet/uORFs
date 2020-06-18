home <- "/home/fcalvet/Documents/replicate_yale_analysis/data/reference_sets/"
science_non_canonical <- read.table(paste(home, "science_non_canonicalORFs_matched_with_candidates.tsv", sep = ""))
human_heart <- read.table(paste(home, "human_heart_uORFs_matched_with_candidates_proper_format.tsv", sep = ""))
uorfs_tools <- read.table(paste(home, "uORF_tools_matched_with_candidates.tsv", sep = ""))

science_non_canonical$location <- paste(paste(paste(science_non_canonical$V1, science_non_canonical$V2, sep = ":"), science_non_canonical$V3,
                                 sep = "-"), science_non_canonical$V4, sep = ":")

human_heart$location <- paste(paste(paste(human_heart$V1, human_heart$V2, sep = ":"), human_heart$V3,
                                              sep = "-"), human_heart$V4, sep = ":")

uorfs_tools$location <- paste(paste(paste(uorfs_tools$V1, uorfs_tools$V2, sep = ":"), uorfs_tools$V3,
                                              sep = "-"), uorfs_tools$V4, sep = ":")

all_ids <- union(union(science_non_canonical$location, human_heart$location),uorfs_tools$location)


# install.packages("limma")
library(limma)

check_reference <- function(x, y){
  # this function checks if a value (x) is in a vector (y)
  return (x %in% y)
}


data <- data.frame(matrix(nrow=length(all_ids), ncol=0))
data$science <- as.integer(sapply(all_ids, check_reference, y=science_non_canonical$location))
data$heart <- as.integer(sapply(all_ids, check_reference, y=human_heart$location))
data$tools <- as.integer(sapply(all_ids, check_reference, y=uorfs_tools$location))
a <- vennCounts(data)
a
vennDiagram(a, main = "Reference uORFs", cex=0.7, show.include = T)



