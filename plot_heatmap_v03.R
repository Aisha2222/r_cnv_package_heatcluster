library(tidyverse)
library(copynumber)
final_bed_25_samples <- read_tsv("intersected_allsamp_v05.txt", col_types = "ciicdc", na =".") #read in dataframe from intersection and the specify that chrom is character using col_types to prevent X and Y from turning to NA
final_bed_25_samples$n.probes=10 #add in the column for probes
final_bed_25_samples_v02 <- final_bed_25_samples [, c(6,1,4,2,3,7,5)] #rearrange the columns to fit the copy number data frame format 
table(final_bed_25_samples_v02[,1]) #to check if i have NA in the sample id column
dist_chrom <- distinct(final_bed_25_samples_v02, sampleID)
#####at this point I have seen that there are some rows without sample name, so I am going to use bedtools to cut thos rows of using the sed command, I will export the data first with write_tsv####
write_tsv(final_bed_25_samples_v02, "preprocessed_final_final_bed_25_samples_v02.txt")
#i used the command grep -v '^NA' preprocessed_final_final_bed_25_samples_v02.txt > linux_processed_final_final_bed_25_samples_v02.txt to remove the columns with NA sample ID because I know that wont be tolerated by the plot_heatmap function
#now i am reading in the modified final_bed_25_samples_v02 which i named linux_processed_final_final_bed_25_samples_v02.txt 

final_bed_25_samples_v03 <- read_tsv("newest folder for processed data for final analysis/linux_processed_final_final_bed_25_samples_v02.txt", col_types = "ccciidd") #it is important that you specify coltypes otherwise R will return X and Y as NA

####data quality check before moving on####
dist_sampleID<- distinct(final_bed_25_samples_v03, sampleID) #just 25 samples yay!
dist_chrom<- distinct(final_bed_25_samples_v03, chrom) 
dist_mean<-  distinct(final_bed_25_samples_v03, mean) 
dist_n.probes <- distinct(final_bed_25_samples_v03, n.probes)
dist_arm <- distinct(final_bed_25_samples_v03, arm)

#great, data looks good so I can now moove on to correct for ploidy and log tranformation with mean column (CN)

final_bed_25_samples_v04 <- final_bed_25_samples_v03 
final_bed_25_samples_v04$mean <- log2(final_bed_25_samples_v04$mean/2) #replace the mean column with log 2 of the value correct for ploidy 
table(final_bed_25_samples_v04[,7]) #see how the copy numbers spread for each of the sample to help with determining what my limit will be for the plotheatmap function

##check for NA's comparing before and after log transformation##
dist_mean_v04 <- distinct(final_bed_25_samples_v04, mean)
dist_mean_v03 <- distinct(final_bed_25_samples_v03, mean)

#first I have to solve the -inf problem because the plotheatmap function cant handle this, have to convert the mean at -inf to  -2
final_bed_25_samples_v05<- final_bed_25_samples_v04 %>% 
  mutate(mean = ifelse(mean == -Inf, -2, mean)) 
table(final_bed_25_samples_v05[,7]) #to confirm if my conversion to -2 truly worked 


#then now I can make my plotheatmap function 
plotHeatmap(as.data.frame(final_bed_25_samples_v05), upper.lim=c(1,5), lower.lim = c(-2,-1)) #the limit that you choose for the plotheatmap really depends on the on the range of values for copy numbers where your sample lie
plotHeatmap(as.data.frame(final_bed_25_samples_v05), upper.lim= -0, lower.lim = 5)
