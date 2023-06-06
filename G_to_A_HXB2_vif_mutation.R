#library importation
library(readr)
library(data.table)
library(reshape2)
library(tidyverse)

#Parameters
#######ION TORRENT, 500 for min, 3 for min freq alt
min_coverage_major = 1000
min_coverage_minor = 1000
min_freq_alt = 5
max_freq_alt = 49.99999
min_freq_indel = 30
#start and end of the gene
start_pol = 5041
end_pol = 5619

#import reference sequence
ref_seq <- read.table("vif_reference.txt", header = TRUE, sep = "\t")
#create dataframe for hypermut analysis
hypermut = data.frame(Id_sample=character(), nb_motif=integer(), nb_control=integer())

#import all tsv files
filenames = dir(pattern="*.tsv")

#function
'%!in%'<- function(x,y)!('%in%'(x,y))

#for each sample
for (my_file in filenames){
  print(my_file)
  #get name file
  file_name = sub("\\..*", "", my_file)
  #c for character, i for integer
  my_sample <- read_tsv(my_file, col_names = TRUE,
                        col_types = c("c", "i", "i", "i", "i", "i", "i",
                                      "c", "c", "c"))
  
  #add lacking nucleotide
  for(i in start_pol:end_pol){
    if(i %!in% my_sample$Position){
      my_sample[nrow(my_sample) + 1,] = list(my_file, i, 0, 0, 0, 0, 0, "N", "N", NA)
    }
  }
  my_sample = my_sample[order(my_sample$Position),]
  
  #fixing nucleotide with zero depth
  for(i in 1:nrow(my_sample)){
    if(my_sample$n_A[i] == 0 && my_sample$n_T[i] == 0 && my_sample$n_G[i] == 0 && my_sample$n_C[i] == 0){
      my_sample$n_A[i] = 1
    }
  }
  
  #initialize variables
  my_sample$complement <- NULL
  my_sample$qual <- NULL
  
  #table combination
  full_data <- cbind(my_sample, ref_seq)
  full_data$depth_total <- rowSums(full_data[3:7])
  full_data$major_depth <- apply(full_data[3:7], 1, FUN=max)
  
  #for each position of the gene
  for (i in 1:nrow(full_data)){
    #calculate freq of each base and indel
    full_data$freq_A[i] <- full_data$n_A[i]/full_data$depth_total[i]*100
    full_data$freq_T[i] <- full_data$n_T[i]/full_data$depth_total[i]*100
    full_data$freq_G[i] <- full_data$n_G[i]/full_data$depth_total[i]*100
    full_data$freq_C[i] <- full_data$n_C[i]/full_data$depth_total[i]*100
    full_data$freq_indel[i] <- full_data$n_indel[i]/full_data$depth_total[i]*100
    #calculate freq of major and minor alleles
    full_data$freq_alt[i] <- max(full_data[i,14:18][full_data[i,14:18] != max(full_data[i,14:18])])
    full_data$freq_maj[i] <- max(full_data[i,14:18])
    
    #define the consensus base
    if(full_data$depth_total[i] >= min_coverage_major){
      if(full_data$freq_A[i] == full_data$freq_maj[i]){
        full_data$major_seq[i] = "A"
      } else if(full_data$freq_T[i] == full_data$freq_maj[i]){
        full_data$major_seq[i] = "T"
      } else if(full_data$freq_G[i] == full_data$freq_maj[i]){
        full_data$major_seq[i] = "G"
      } else if(full_data$freq_C[i] == full_data$freq_maj[i]){
        full_data$major_seq[i] = "C"
      } else {
        full_data$major_seq[i] = "del"
      }
    } else {
      full_data$major_seq[i] = "N"
    }
    
    #define the minor base if possible
    if(full_data$freq_alt[i] == full_data$freq_indel[i]){
      if(full_data$freq_alt[i] < min_freq_indel && full_data$depth_total[i] >= min_coverage_major){
        full_data$minor_seq[i] = full_data$major_seq[i]
      } else if (full_data$freq_indel[i] >= min_freq_indel && full_data$freq_indel[i] < max_freq_alt){
        full_data$minor_seq[i] = "del"
      } else {
        full_data$minor_seq[i] = "N"
      }
    } else {
      if(full_data$freq_alt[i] < min_freq_alt && full_data$depth_total[i] >= min_coverage_major){
        full_data$minor_seq[i] = full_data$major_seq[i]
      } else if (full_data$depth_total[i] >= min_coverage_minor && full_data$freq_alt[i] >= min_freq_alt){
        if(full_data$freq_A[i] >= min_freq_alt && full_data$freq_A[i] < max_freq_alt){
          full_data$minor_seq[i] = "A"
        } else if (full_data$freq_T[i] >= min_freq_alt && full_data$freq_T[i] < max_freq_alt){
          full_data$minor_seq[i] = "T"
        } else if (full_data$freq_G[i] >= min_freq_alt && full_data$freq_G[i] < max_freq_alt){
          full_data$minor_seq[i] = "G"
        } else if (full_data$freq_C[i] >= min_freq_alt && full_data$freq_C[i] < max_freq_alt){
          full_data$minor_seq[i] = "C"
        } 
      } else {
        full_data$minor_seq[i] = "N"
      }
    }
    
    #for each codon, check if sufficient data for minor mutation
    if(i%%3 == 0){
      if(full_data$minor_seq[i-2] == "N" || full_data$minor_seq[i-1] == "N" || full_data$minor_seq[i] == "N"){
        full_data$N_minor[i] = "too_low"
      } else if (full_data$minor_seq[i-2] == "del" || full_data$minor_seq[i-1] == "del" || full_data$minor_seq[i] == "del"){
        full_data$N_minor[i] = "del"
      } else {
        full_data$N_minor[i] = "pass"      
      }
    } else{
      full_data$N_minor[i] = NA
    }
    
    #for each codon, check if sufficient data for major mutation
    if(i%%3 == 0){
      if(full_data$major_seq[i-2] == "N" || full_data$major_seq[i-1] == "N" || full_data$major_seq[i] == "N"){
        full_data$N_major[i] = "too_low"
      } else if (full_data$major_seq[i-2] == "del" || full_data$major_seq[i-1] == "del" || full_data$major_seq[i] == "del"){
        full_data$N_major[i] = "del"
      } else {
        full_data$N_major[i] = "pass"      
      }
    } else {
      full_data$N_major[i] = NA
    }
    
    #perform translation for each major and minor codons
    #also check for GRD codon motif
    if(i%%3 == 0){
      if(full_data$minor_seq[i-2] == "T"){
        if(full_data$minor_seq[i-1] == "T"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "F"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "F"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "L"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "L"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "C"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "S"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "S"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "S"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "S"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "A"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "Y"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "Y"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "*"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "*"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "G"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "C"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "C"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "*"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "W"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-2] == "C"){
        if(full_data$minor_seq[i-1] == "T"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "L"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "L"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "L"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "L"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "C"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "P"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "P"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "P"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "P"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "A"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "H"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "H"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "Q"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "Q"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "G"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "R"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "R"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "R"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "R"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-2] == "A"){
        if(full_data$minor_seq[i-1] == "T"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "I"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "I"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "I"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "M"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "C"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "T"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "T"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "T"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "T"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "A"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "N"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "N"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "K"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "K"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "G"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "S"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "S"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "R"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "R"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-2] == "G"){
        if(full_data$minor_seq[i-1] == "T"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "V"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "V"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "V"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "V"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "C"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "A"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "A"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "A"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "A"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "A"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "D"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "D"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "E"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "E"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "G"){
          if(full_data$minor_seq[i] == "T"){
            full_data$minor[i] = "G"
          } else if(full_data$minor_seq[i] == "C"){
            full_data$minor[i] = "G"
          } else if(full_data$minor_seq[i] == "A"){
            full_data$minor[i] = "G"
          } else if(full_data$minor_seq[i] == "G"){
            full_data$minor[i] = "G"
          } else if(full_data$minor_seq[i] == "del") {
            full_data$minor[i] = "del"
          } else {
            full_data$minor[i] = "too_low"
          }
        } else if(full_data$minor_seq[i-1] == "del") {
          full_data$minor[i] = "del"
        } else {
          full_data$minor[i] = "too_low"
        }
      } else if(full_data$minor_seq[i-2] == "del") {
        full_data$minor[i] = "del"
      } else {
        full_data$minor[i] = "too_low"
      }
    } else {
      full_data$minor[i] = NA
    }
    
    if(i%%3 == 0){
      if(full_data$major_seq[i-2] == "T"){
        if(full_data$major_seq[i-1] == "T"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "F"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "F"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "L"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "L"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "C"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "S"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "S"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "S"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "S"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "A"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "Y"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "Y"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "*"
            full_data$GDR[i] = "STOP"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "*"
            full_data$GDR[i] = "STOP"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "G"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "C"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "C"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "*"
            full_data$GDR[i] = "STOP"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "W"
            full_data$GDR[i] = "STOP possible"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "del") {
          full_data$major[i] = "del"
          full_data$GDR[i] = "No"
        } else {
          full_data$major[i] = "too_low"
          full_data$GDR[i] = "No"
        }
      } else if(full_data$major_seq[i-2] == "C"){
        if(full_data$major_seq[i-1] == "T"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "L"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "L"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "L"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "L"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "C"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "P"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "P"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "P"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "P"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "A"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "H"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "H"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "Q"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "Q"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "G"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "R"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "R"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "R"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "R"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "del") {
          full_data$major[i] = "del"
          full_data$GDR[i] = "No"
        } else {
          full_data$major[i] = "too_low"
          full_data$GDR[i] = "No"
        }
      } else if(full_data$major_seq[i-2] == "A"){
        if(full_data$major_seq[i-1] == "T"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "I"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "I"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "I"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "M"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "C"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "T"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "T"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "T"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "T"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "A"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "N"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "N"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "K"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "K"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "G"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "S"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "S"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "R"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "R"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "del") {
          full_data$major[i] = "del"
          full_data$GDR[i] = "No"
        } else {
          full_data$major[i] = "too_low"
          full_data$GDR[i] = "No"
        }
      } else if(full_data$major_seq[i-2] == "G"){
        if(full_data$major_seq[i-1] == "T"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "V"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "V"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "V"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "V"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "C"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "A"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "A"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "A"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "A"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "A"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "D"
            full_data$GDR[i] = "G>A possible"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "D"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "E"
            full_data$GDR[i] = "G>A possible"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "E"
            full_data$GDR[i] = "G>A possible"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "G"){
          if(full_data$major_seq[i] == "T"){
            full_data$major[i] = "G"
            full_data$GDR[i] = "G>A possible"
          } else if(full_data$major_seq[i] == "C"){
            full_data$major[i] = "G"
            full_data$GDR[i] = "No"
          } else if(full_data$major_seq[i] == "A"){
            full_data$major[i] = "G"
            full_data$GDR[i] = "G>A possible"
          } else if(full_data$major_seq[i] == "G"){
            full_data$major[i] = "G"
            full_data$GDR[i] = "G>A possible"
          } else if(full_data$major_seq[i] == "del") {
            full_data$major[i] = "del"
            full_data$GDR[i] = "No"
          } else {
            full_data$major[i] = "too_low"
            full_data$GDR[i] = "No"
          }
        } else if(full_data$major_seq[i-1] == "del") {
          full_data$major[i] = "del"
          full_data$GDR[i] = "No"
        } else {
          full_data$major[i] = "too_low"
          full_data$GDR[i] = "No"
        }
      } else if(full_data$major_seq[i-2] == "del") {
        full_data$major[i] = "del"
        full_data$GDR[i] = "No"
      } else {
        full_data$major[i] = "too_low"
        full_data$GDR[i] = "No"
      }
    } else {
      full_data$major[i] = NA
      full_data$GDR[i] = "No"
    }
    
    #definitive observation for minor codon
    if(i%%3 == 0){
      if(full_data$N_minor[i] == "pass"){
        full_data$final_minor[i] = full_data$minor[i]
      } else if (full_data$N_minor[i] == "del") {
        if(full_data$minor_seq[i-1] == "del" && full_data$minor_seq[i-2] == "del"){
          full_data$final_minor[i] = "del"
        }else{
          full_data$final_minor[i] = full_data$AA_ref[i]
        }
      } else {
        full_data$final_minor[i] = "too_low"
      }
      
      #definitive observation for major codon
      if(full_data$N_major[i] == "pass"){
        full_data$final_major[i] = full_data$major[i]
      } else if (full_data$N_major[i] == "del") {
        if(full_data$major_seq[i-1] == "del" && full_data$major_seq[i-2] == "del"){
          full_data$final_major[i] = "del"
        }else{
          full_data$final_major[i] = full_data$AA_ref[i]
        }
      } else {
        full_data$final_major[i] = "too_low"
      }
    }
    
    #Sample APOBEC detection
    #Context GRD, at least 5% A, and G > A not observed on A at ref
    #we discard the two first bases
    #looking for minor APOBEc mutation
    if(i>2){
      if(full_data$major_seq[i-2] == "G"){
        if(full_data$major_seq[i-1] == "A" || full_data$major_seq[i-1] == "G"){
          if(full_data$major_seq[i] == "A" || full_data$major_seq[i] == "G" || full_data$major_seq[i] == "T") {
            if(full_data$major_seq[i-2] == "G" && full_data$freq_A[i-2] >= min_freq_alt && full_data$nucl_ref[i-2] != "A"){
              full_data$sample_Apobec[i-2] = "Apobec minor"
              full_data$sample_Apobec_freq[i-2] = full_data$freq_A[i-2]
            } else {
              full_data$sample_Apobec[i-2] = "No"
              full_data$sample_Apobec_freq[i-2] = 0
            }
          } else {
            full_data$sample_Apobec[i-2] = "No"
            full_data$sample_Apobec_freq[i-2] = 0
          }
        } else {
          full_data$sample_Apobec[i-2] = "No"
          full_data$sample_Apobec_freq[i-2] = 0
        } 
      } else {
        full_data$sample_Apobec[i-2] = "No"
        full_data$sample_Apobec_freq[i-2] = 0
      }
    }
    else {
      full_data$sample_Apobec[i] = "No"
      full_data$sample_Apobec_freq[i] = 0
    }
    
    #Check for APOBEC mutations compared to the reference sequence
    if(i>2){
      if(full_data$nucl_ref[i-2] == full_data$major_seq[i-2] && full_data$nucl_ref[i-1] == full_data$major_seq[i-1] && full_data$nucl_ref[i] == full_data$major_seq[i]){
        full_data$fixed_GRD[i-2] = "No"
        full_data$fixed_GRD_freq[i-2] = 0
      }
      else {
        if(full_data$nucl_ref[i-2] == "G"){
          if(full_data$nucl_ref[i-1] == "A" || full_data$nucl_ref[i-1] == "G"){
            if(full_data$nucl_ref[i] == "A" || full_data$nucl_ref[i] == "G" || full_data$nucl_ref[i] == "T") {
              if(full_data$nucl_ref[i-2] == "G" && full_data$major_seq[i-2] == "A"){
                full_data$fixed_GRD[i-2] = "Apobec"
                full_data$fixed_GRD_freq[i-2] = full_data$freq_A[i-2]
              } else {
                full_data$fixed_GRD[i-2] = "No"
                full_data$fixed_GRD_freq[i-2] = 0
              }
            }else {
              full_data$fixed_GRD[i-2] = "No"
              full_data$fixed_GRD_freq[i-2] = 0
            }
          }else {
            full_data$fixed_GRD[i-2] = "No"
            full_data$fixed_GRD_freq[i-2] = 0
          }
        }else {
          full_data$fixed_GRD[i-2] = "No"
          full_data$fixed_GRD_freq[i-2] = 0
        }
      }
    }
    else {
      full_data$fixed_GRD[i] = "No"
      full_data$fixed_GRD_freq[i] = 0
    }
    
    #STOP major codon detection
    if(full_data$depth_total[i] >= min_coverage_major){
      if(i%%3 == 0){
        if(full_data$major_seq[i-2] == "T"){
          if(full_data$major_seq[i-1] == "G"){
            if(full_data$major_seq[i] == "A"){
              full_data$STOP_sample[i-2] = "STOP"
              full_data$STOP_sample_freq[i-2] = min(full_data$freq_T[i-2], full_data$freq_G[i-1],full_data$freq_A[i])
            } else {
              full_data$STOP_sample[i-2] = "No"
              full_data$STOP_sample_freq[i-2] = 0
            }
          } else if(full_data$major_seq[i-1] == "A"){
            if(full_data$major_seq[i] == "A" || full_data$major_seq[i] == "G"){
              full_data$STOP_sample[i-2] = "STOP"
              if(full_data$major_seq[i] == "A"){
                full_data$STOP_sample_freq[i-2] = min(full_data$freq_T[i-2], full_data$freq_A[i-1],full_data$freq_A[i])
              } else {
                full_data$STOP_sample_freq[i-2] = min(full_data$freq_T[i-2], full_data$freq_A[i-1],full_data$freq_G[i])
              }
            } else {
              full_data$STOP_sample[i-2] = "No"
              full_data$STOP_sample_freq[i-2] = 0
            }
          } else {
            full_data$STOP_sample[i-2] = "No"
            full_data$STOP_sample_freq[i-2] = 0
          }
        } else {
          full_data$STOP_sample[i-2] = "No"
          full_data$STOP_sample_freq[i-2] = 0
        }
      } else {
        full_data$STOP_sample[i] = "No"
        full_data$STOP_sample_freq[i] = 0
      }
    } else {
      full_data$STOP_sample[i] = "No"
      full_data$STOP_sample_freq[i] = 0
    }

    #STOP minor codon detection
    if(full_data$depth_total[i] >= min_coverage_major){
      if(i%%3 == 0){
        if(full_data$STOP_sample[i-2] == "STOP"){
          full_data$STOP_sample2[i-2] = "No"
          full_data$STOP_sample2_freq[i-2] = 0
        } else {
          if(full_data$freq_T[i-2] >= min_freq_alt){
            if(full_data$freq_G[i-1] >= min_freq_alt){
              if(full_data$freq_A[i] >= min_freq_alt){
                full_data$STOP_sample2[i-2] = "STOP minor"
                full_data$STOP_sample2_freq[i-2] = min(full_data$freq_T[i-2], full_data$freq_G[i-1],full_data$freq_A[i])
              } else {
                full_data$STOP_sample2[i-2] = "No"
                full_data$STOP_sample2_freq[i-2] = 0
              }
            } else if(full_data$freq_A[i-1] >= min_freq_alt){
              if(full_data$freq_A[i] >= min_freq_alt || full_data$freq_G[i] >= min_freq_alt){
                full_data$STOP_sample2[i-2] = "STOP minor"
                if(full_data$freq_A[i] >= min_freq_alt){
                  full_data$STOP_sample2_freq[i-2] = min(full_data$freq_T[i-2], full_data$freq_A[i-1],full_data$freq_A[i])
                } else {
                  full_data$STOP_sample2_freq[i-2] = min(full_data$freq_T[i-2], full_data$freq_A[i-1],full_data$freq_G[i])
                }
              } else {
                full_data$STOP_sample2[i-2] = "No"
                full_data$STOP_sample2_freq[i-2] = 0
              }
            } else {
              full_data$STOP_sample2[i-2] = "No"
              full_data$STOP_sample2_freq[i-2] = 0
            }
          } else {
            full_data$STOP_sample2[i-2] = "No"
            full_data$STOP_sample2_freq[i-2] = 0
          }
        }
      } else {
        full_data$STOP_sample2[i] = "No"
        full_data$STOP_sample2_freq[i] = 0
      }
    }else {
      full_data$STOP_sample2[i] = "No"
      full_data$STOP_sample2_freq[i] = 0
    }
    
    #we discard insufficiently covered codons for the search of minor mutations
    if(i%%3 == 0){
      if(full_data$final_major[i] == "too_low"){
        full_data$mutation[i] = "discarded"
      } else if(full_data$final_major[i] != full_data$AA_ref[i]) {
        full_data$mutation[i] = gsub(" ", "", paste(full_data$AA_ref[i],full_data$pos_pol_AA[i],full_data$final_major[i]))
      } else if(full_data$final_minor[i] == "too_low"){
        full_data$mutation[i] = "discarded_minor"
      } else if(full_data$final_minor[i] != full_data$AA_ref[i]) {
        full_data$mutation[i] = gsub(" ", "", paste(full_data$AA_ref[i],full_data$pos_pol_AA[i],full_data$final_minor[i]))
      } else {
        full_data$mutation[i] = "no_mutation"
      }
    }
    
    #we check if the mutation is non-synonymous
    #We do not consider synonymous mutations
    if(i%%3 == 0){
      if(full_data$mutation[i] != "no_mutation"){
        if(full_data$mutation[i] == "discarded"){
          full_data$obs[i] = "too_low"
        } else if(full_data$mutation[i] == "discarded_minor") {
          full_data$obs[i] = "too_low_minor"
        } else if (full_data$major_seq[i] != full_data$nucl_ref[i] || full_data$major_seq[i-1] != full_data$nucl_ref[i-1] || full_data$major_seq[i-2] != full_data$nucl_ref[i-2]){
          if(full_data$final_major[i] != full_data$AA_ref[i]){
            full_data$obs[i] = "major"
          }
          else {
            full_data$obs[i] = "minor_investigation"
          }
        }
        else {
          full_data$obs[i] = "minor_investigation"
        }
        if(full_data$obs[i] == "minor_investigation"){
          if(full_data$minor_seq[i] != full_data$nucl_ref[i] || full_data$minor_seq[i-1] != full_data$nucl_ref[i-1] || full_data$minor_seq[i-2] != full_data$nucl_ref[i-2]){
            full_data$obs[i] = "minor"
          } else {
            full_data$obs[i] = "-"
          }
        }
      } else {
        full_data$obs[i] = "-"
      }
    }
    
    #we attribute the frequency of the observation
    if(i%%3 == 0){
      if(full_data$obs[i] == "minor"){
        full_data$freq[i] = max(full_data$freq_alt[i-2], full_data$freq_alt[i-1], full_data$freq_alt[i])
      } else if(full_data$obs[i] == "major") {
        full_data$freq[i] = max(full_data$freq_maj[i-2], full_data$freq_maj[i-1], full_data$freq_maj[i])
      } else if(full_data$obs[i] == "too_low_minor"){
        full_data$freq[i] = "too_low_minor"
      } else if(full_data$obs[i] == "too_low"){
        full_data$freq[i] = "too_low"
      } else {
        full_data$freq[i] = NA
      }
    }
    
    #we attribute the depth related to the mutation
    if(i%%3 == 0){
      if(full_data$obs[i] == "minor" || full_data$obs[i] == "major"){
        if(full_data$major_seq[i-2] != full_data$nucl_ref[i-2] || full_data$minor_seq[i-2] != full_data$nucl_ref[i-2]){
          full_data$depth[i] = full_data$depth_total[i-2]
        } else if (full_data$major_seq[i-1] != full_data$nucl_ref[i-1] || full_data$minor_seq[i-1] != full_data$nucl_ref[i-1]){
          full_data$depth[i] = full_data$depth_total[i-1]
        } else if (full_data$major_seq[i] != full_data$nucl_ref[i] || full_data$minor_seq[i] != full_data$nucl_ref[i]){
          full_data$depth[i] = full_data$depth_total[i]
        }
      } else {
        full_data$depth[i] = NA
      }
    } else {
      full_data$depth[i] = NA
    }
    
    #we check for major G>A mutations in the APOBEC context
    if(i%%3 == 0){
      if(full_data$nucl_ref[i-2] == "G"){
        if(full_data$nucl_ref[i-1] == "A" || full_data$nucl_ref[i-1] == "G"){
          if(full_data$nucl_ref[i] == "A" || full_data$nucl_ref[i] == "G" || full_data$nucl_ref[i] == "T") {
            if(full_data$major_seq[i-2] == "A"){
              full_data$GDR[i-2] = "G>A"
              full_data$GDR[i-1] = "G>A"
              full_data$GDR[i] = "G>A"
            } 
          }
        }
      }
    }
    
    #we update the data for each position of the codon
    if(i%%3 == 0){
      if(full_data$GDR[i] == "STOP"){
        full_data$GDR[i-1] = "STOP"
        full_data$GDR[i-2] = "STOP"
      } else if (full_data$GDR[i] == "G>A possible"){
        full_data$GDR[i-1] = "G>A possible"
        full_data$GDR[i-2] = "G>A possible"
      } else if(full_data$GDR[i] == "STOP possible"){
        full_data$GDR[i-1] = "STOP possible"
        full_data$GDR[i-2] = "STOP possible"
      } else if(full_data$GDR[i] == "G>A"){
        full_data$GDR[i-1] = "G>A"
        full_data$GDR[i-2] = "G>A"
      }
    }
  }
  
  #for each position of the gene
  for (i in 1:nrow(full_data)){
    
    #we check for G>A mutations NOT IN the APOBEC context
    if(full_data$fixed_GRD[i] != "Apobec"){
      if(full_data$nucl_ref[i] == "G" && full_data$major_seq[i] == "A"){
        full_data$fixed_NON_GRD[i] = "Non Apobec"
        full_data$fixed_NON_GRD_freq[i] = full_data$freq_A[i]
      } else {
        full_data$fixed_NON_GRD[i] = "No"
        full_data$fixed_NON_GRD_freq[i] = 0
      }
    } else {
      full_data$fixed_NON_GRD[i] = "No"
      full_data$fixed_NON_GRD_freq[i] = 0
    }
    
    #Non APOBEC minor mutation detection
    #at least 5% A, and G > A not observed on A at ref, no GRD
    if(full_data$sample_Apobec[i] != "Apobec"){
      if(full_data$major_seq[i] == "G" && full_data$freq_A[i] >= min_freq_alt && full_data$nucl_ref[i] != "A") {
        full_data$sample_NON_Apobec[i] = "Non Apobec minor"
        full_data$sample_NON_Apobec_freq[i] = full_data$freq_A[i]
      } else {
        full_data$sample_NON_Apobec[i] = "No"
        full_data$sample_NON_Apobec_freq[i] = 0
      }
    }
    else {
      full_data$sample_NON_Apobec[i] = "No"
      full_data$sample_NON_Apobec_freq[i] = 0
    }
  }
  
  #for each position, we attribute the final observation
  for (i in 1:nrow(full_data)){
    if(full_data$fixed_GRD[i] == "Apobec"){
      full_data$FINAL[i] = "Apobec"
      full_data$FINAL_freq[i] = full_data$fixed_GRD_freq[i]
    } else if(full_data$sample_Apobec[i] == "Apobec minor"){
      full_data$FINAL[i] = "Apobec minor"
      full_data$FINAL_freq[i] = full_data$sample_Apobec_freq[i]
    } else if(full_data$fixed_NON_GRD[i] == "Non Apobec"){
      full_data$FINAL[i] = "Non Apobec"
      full_data$FINAL_freq[i] = full_data$fixed_NON_GRD_freq[i]
    } else if(full_data$sample_NON_Apobec[i] == "Non Apobec minor") {
      full_data$FINAL[i] = "Non Apobec minor"
      full_data$FINAL_freq[i] = full_data$sample_NON_Apobec_freq[i]
    } else if(full_data$STOP_sample[i] == "STOP"){
      full_data$FINAL[i] = "STOP"
      full_data$FINAL_freq[i] = full_data$STOP_sample_freq[i]
    } else if(full_data$STOP_sample2[i] == "STOP minor"){
      full_data$FINAL[i] = "STOP minor"
      full_data$FINAL_freq[i] = full_data$STOP_sample2_freq[i]
    } else {
      full_data$FINAL[i] = "No"
      full_data$FINAL_freq[i] = 0
    }
  }
  
  #save all the data for each sample
  write.table(x = full_data, file = paste(my_file, "_done.txt", sep=""), sep="\t")
  
  #perform the hypermut analysis
  nb_motif = 0
  nb_control = 0
  possible_GRD=c("GAA","GAG","GAT","GGA","GGG","GGT")
  for(my_position in 1:nrow(full_data)){
    if(my_position>2){
      pos1_ref = full_data$nucl_ref[my_position-2]
      pos2_ref = full_data$nucl_ref[my_position-1]
      pos3_ref = full_data$nucl_ref[my_position]
      pos1_seq = full_data$major_seq[my_position-2]
      pos2_seq = full_data$major_seq[my_position-1]
      pos3_seq = full_data$major_seq[my_position]
      motif_ref = paste(pos1_ref, pos2_ref, pos3_ref, "")
      motif_ref = gsub(" ", "", motif_ref, fixed = T)
      if(motif_ref %in% possible_GRD){
        if(pos1_seq %in% c("A","T","G","C") || pos2_seq %in% c("A","T","G","C") && pos3_seq %in% c("A","T","G","C")){
          nb_motif = nb_motif + 1
        }
      }else{
        if(pos1_ref=="G")
          if(pos1_seq %in% c("A","T","G","C")){
            nb_control = nb_control + 1
          }
      }
    }
  }
  
  #complete the table
  hypermut[nrow(hypermut) + 1,] = list(file_name, nb_motif, nb_control)
}

#write the hypermut result in a file
write.table(x = hypermut, file = "hypermut.txt", sep="\t", row.names = F)

####produce fasta consensus sequence for each sample
temp = list.files(pattern="*_done.txt")
for (my_file in temp){
  file_name = sub("\\_.*", "", my_file)
  print(my_file)
  myfiles = lapply(my_file, read.delim)
  all_data = rbindlist(myfiles, fill=FALSE, idcol=NULL)
  for(i in 1:nrow(all_data)){
    if(all_data$major_seq[i]=="del"){
      all_data$major_seq[i]="-"
    }
  }
  my_seq = all_data$major_seq
  my_seq = paste(my_seq,collapse="")
  final_df = data.frame(1)
  final_df[nrow(final_df) + 1,] = gsub(" ", "", paste(">", my_file))
  final_df[nrow(final_df) + 1,] = my_seq
  final_df = final_df[-1,]
  #save the fasta in a file
  write.table(x = final_df, file = paste(file_name, ".fasta", sep=""),
              sep="\t", col.names = F, row.names = F, quote = F)
}

#concatenate all results
temp = list.files(pattern="*_done.txt")
myfiles = lapply(temp, read.delim)
all_data = rbindlist(myfiles, fill=FALSE, idcol=NULL)

#retain only specific columns
my_sub_data = all_data[,c(1,2,9,10,11,39,40,41,46,47)]
my_sub_data = my_sub_data[complete.cases(my_sub_data), ]
my_sub_data = my_sub_data[my_sub_data$freq!="too_low",]

#we consider a mutation as major when the frequency is higher than 20%
for(i in 1:nrow(my_sub_data)){
  if(as.numeric(my_sub_data$freq[i])>=20){
    if(my_sub_data$obs[i]=="major" || my_sub_data$obs[i]=="minor"){
      my_sub_data$obs2[i]="major"
    }
  }else{
    my_sub_data$obs2[i]=my_sub_data$obs[i]
  }
}

my_sub_data = my_sub_data[my_sub_data$freq>1,]
my_sub_data = my_sub_data[my_sub_data$freq!="too_low",]

#generate a matrix showing the frequency of mutations according to patients
matrix_mutation = dcast(my_sub_data, Id_sample ~ mutation, value.var="freq")
matrix_mutation = matrix_mutation[,order(parse_number(colnames(matrix_mutation)))]
write.table(matrix_mutation, file="list_mutations_VIH_HIV1.txt", col.names = TRUE, sep="\t", row.names = FALSE)
my_sub_data_save = my_sub_data[,c(1,2,4,5,6,8,11)]
write.table(my_sub_data_save, file="results_mutations_VIH_HIV1.txt", col.names = TRUE, sep="\t", row.names = FALSE)

#retain specific columns for APOBEC analysis
GDR_data = all_data[,c(1,2,9,10,11,39,40,41,46,47)]

#generate a count of apobec mutations and stop codons
Apobec_GDR = GDR_data[GDR_data$FINAL == "Apobec",]
Apobec_minor_GDR = GDR_data[GDR_data$FINAL == "Apobec minor",]
Non_Apobec_GDR = GDR_data[GDR_data$FINAL == "Non Apobec",]
Non_Apobec_minor_GDR = GDR_data[GDR_data$FINAL == "Non Apobec minor",]
STOP_codon = GDR_data[GDR_data$FINAL == "STOP",]
STOP_codon_minor = GDR_data[GDR_data$FINAL == "STOP minor",]
GDR_data = rbind(Apobec_GDR, Apobec_minor_GDR, Non_Apobec_GDR, Non_Apobec_minor_GDR, STOP_codon, STOP_codon_minor)
GDR_data2 = GDR_data[!duplicated(GDR_data)]
GDR_data2 = GDR_data2[,c(1,5,9,10)]
GDR_data2_clean = GDR_data2[!duplicated(GDR_data2)]

#we consider an APOBEC mutation or a STOP codon as major when the frequency is higher than 20%
for(i in 1:nrow(GDR_data2_clean)){
  if(GDR_data2_clean$FINAL_freq[i]>=20){
    if(GDR_data2_clean$FINAL[i]=="Apobec minor" || GDR_data2_clean$FINAL[i]=="Apobec"){
      GDR_data2_clean$FINAL2[i]="Apobec"
    } else if(GDR_data2_clean$FINAL[i]=="Non Apobec" || GDR_data2_clean$FINAL[i]=="Non Apobec minor"){
      GDR_data2_clean$FINAL2[i]="Non Apobec"
    } else if(GDR_data2_clean$FINAL[i]=="STOP" || GDR_data2_clean$FINAL[i]=="STOP minor"){
      GDR_data2_clean$FINAL2[i]="STOP"
    }
    
  }else{
    GDR_data2_clean$FINAL2[i]=GDR_data2_clean$FINAL[i]
  }
}

#generate a matrix of APOBEC, non-APOBEC and stop codons
matrix_GDR = dcast(GDR_data2_clean, Id_sample ~ FINAL2)
write.table(matrix_GDR, file="list_apobec.txt", col.names = TRUE, sep="\t", row.names = FALSE)
write.table(GDR_data2_clean, file="details_mutations.txt", col.names = TRUE, sep="\t", row.names = FALSE)

#Hypermut analysis and statistical test
hypermut = read.table("hypermut.txt", sep="\t", header=T)
hypermut$Id_sample = sub("\\_.*", "", hypermut$Id_sample)

matrix_GDR$Id_sample = sub("\\..*", "", matrix_GDR$Id_sample)
matrix_GDR$Id_sample = sub("\\_.*", "", matrix_GDR$Id_sample)

hypermut_analysis = merge(matrix_GDR, hypermut, by="Id_sample")
hypermut_analysis = hypermut_analysis[,c(1,2,3,4,5,8,9)]

hypermut_analysis$A_major = hypermut_analysis$Apobec
hypermut_analysis$A_minor = (hypermut_analysis$Apobec + hypermut_analysis$`Apobec minor`)
hypermut_analysis$B_major = (hypermut_analysis$nb_motif - hypermut_analysis$Apobec)
hypermut_analysis$B_minor = (hypermut_analysis$nb_motif - (hypermut_analysis$Apobec + hypermut_analysis$`Apobec minor`))
hypermut_analysis$C_major = hypermut_analysis$`Non Apobec`
hypermut_analysis$C_minor = (hypermut_analysis$`Non Apobec` + hypermut_analysis$`Non Apobec minor`)
hypermut_analysis$D_major = (hypermut_analysis$nb_control - hypermut_analysis$`Non Apobec`)
hypermut_analysis$D_minor = (hypermut_analysis$nb_control - (hypermut_analysis$`Non Apobec` + hypermut_analysis$`Non Apobec minor`))

#perform a fisher test for each patient to detect hypermutated viruses
for(i in 1:nrow(hypermut_analysis)){
  current_matrix_major = matrix(c(hypermut_analysis$A_major[i], hypermut_analysis$B_major[i], hypermut_analysis$C_major[i], hypermut_analysis$D_major[i]), nrow=2)
  current_matrix_minor = matrix(c(hypermut_analysis$A_minor[i], hypermut_analysis$B_minor[i], hypermut_analysis$C_minor[i], hypermut_analysis$D_minor[i]), nrow=2)
  my_test_major = fisher.test(current_matrix_major)
  my_test_minor = fisher.test(current_matrix_minor)
  hypermut_analysis$hypermut_major[i] = my_test_major$p.value
  hypermut_analysis$hypermut_minor[i] = my_test_minor$p.value
}

hypermut_analysis=hypermut_analysis[,c(1,16,17)]

#write the pvalue results
write.table(hypermut_analysis, file="hypermut_pval.txt", col.names = TRUE, sep="\t", row.names = FALSE)