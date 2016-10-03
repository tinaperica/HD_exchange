# install.packages("getopt")
library(getopt)  # for setting up the options system
# install.packages("plyr")
library(plyr)
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
library(Biostrings)  ## for matching peptides to the sequence
na.mean <- function (x) {
  return (max(x, na.rm = T))
}
spec = matrix(c(
  'data_dir', 'd', 2, 'character',
  'protein_sequence', 's', 1, 'character',
  'protein_sequence2', 't', 2, 'character',
  'ms_data1', 'm', 1, 'character',
  'ms_data2', 'n', 2, 'character',
  'out_prefix', 'o', 2, 'character'
), byrow = T, ncol = 4)
opt <- getopt(spec)
### run code with these options
## e.g.
## RScript plot_HD_exchange.R --ms_data1 "SOS hdpc peptides.csv" --out_prefix "SOS_hdcp"
if ( is.null(opt$data_dir) ) { opt$data_dir <- "."}  ## if no dir is defined, data needs to be in the same dir as the program
if ( is.null(opt$out_prefix) ) { opt$out_prefix <- "out" }  ## if no out_prefix, default to out
outfilename <- paste(opt$out_prefix, "_HD_exchange_stacked_barplot.pdf", sep = "")
#### temp development code
#opt$data_dir = "~/Documents/HD_exchange_MS/"
#opt$ms_data = 'SOS hdpc peptides.csv'
#######################################
setwd(opt$data_dir)
ms_data1 <- read.csv(opt$ms_data, stringsAsFactors = F)
#### this function will add the D% per sequence position (residue number) and count the position occurence (to calculate average)
parse_ms_data <- function (ms_df_row) {
  sequence <- ms_df_row$Sequence
  sequence_vector <- unlist(strsplit(sequence, split = ""))
  percentD <- ms_df_row[,8]
  return(data.frame("peptide" = sequence ,"time" = ms_df_row[,5], "seq_position" = ms_df_row[,2]:ms_df_row[,3], "amino_acid" = sequence_vector, percentD))
}
parsed_ms_df <- adply(ms_data1, .margins = 1, .fun = parse_ms_data, .expand = F) 
per_position_coverage <- with (parsed_ms_df, aggregate(amino_acid, by=list(seq_position = seq_position, time = time), length))
per_position_mean_Dpercent <- with(parsed_ms_df, aggregate(percentD, by=list(amino_acid = amino_acid, seq_position = seq_position, time = time), na.mean))
bin.labels <- c("<10%","10-20%","20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", ">90%")
per_position_mean_Dpercent$bin.label <- cut(per_position_mean_Dpercent$x, breaks = seq(0,100, 10), labels = bin.labels)
per_position_mean_Dpercent$bin <- cut(per_position_mean_Dpercent$x, breaks = seq(0,100,10), labels = as.numeric(1:10))
seq.breaks <- floor(seq(0, max(per_position_mean_Dpercent$seq_position), length = ceiling(max(per_position_mean_Dpercent$seq_position/50))))
seq_bins <- 1:(length(seq.breaks)-1)
per_position_mean_Dpercent$seq_bin <- cut(per_position_mean_Dpercent$seq_position, breaks= seq.breaks , labels = seq_bins)
names(per_position_mean_Dpercent)[4] <- "mean_Dpercent"
colors <- rev(heat.colors(length(unique(per_position_mean_Dpercent$time)), alpha = 0.5))
opar <- par
pdf(file = outfilename, width = 10)
op<-par(mfrow = c(4,1))
par(mar = c(3,5,2,2))
for (i in seq_bins) {
  temp_pos_bin_subset <- subset(per_position_mean_Dpercent, seq_bin == i)
  for.barplot <- matrix(as.numeric(temp_pos_bin_subset$bin), ncol = length(unique(temp_pos_bin_subset$seq_position)), byrow = T)
  #rownames(for.barplot) <- unique(temp_pos_bin_subset$time)
  #colnames(for.barplot) <- unique(temp_pos_bin_subset$seq_position)
  diff.for.barplot <- diff(for.barplot)
  first.row <- matrix(for.barplot[1,], byrow = T, ncol = length(unique(temp_pos_bin_subset$seq_position)))
  #rownames(first.row) <- min(temp_pos_bin_subset$time)
  final.for.barplot <- rbind(first.row, diff.for.barplot)
  my_barplot <- barplot(final.for.barplot, col = colors, ylim = c(0,10), axes = F)
  axis(1, line = 1, labels = as.character(unique(temp_pos_bin_subset$seq_position)), at = my_barplot, cex.axis = 0.8)
  axis(2, labels = bin.labels, at = seq(0,9), las = 2, cex.axis = 0.7)
}
dev.off()
par(opar)
