# Code by Isabelle Du Plessis, 2024

# R version 4.3.2 (2023-10-31)

require(dplyr)
require(reshape2)

#### Load Data
path = getwd() # Set path

# load sample metadata
metadata = read.csv(paste0(path, "/intermediate_files/metadata.csv"))
#change sample id so that we can compare
metadata$Sample_ID = paste("sample_",metadata$Sample_ID,sep="")
#change order so that we can compare
metadata = metadata[order(metadata$Sample_ID),]

# load votu data
votus_cov75thres = read.delim(paste0(path,"/intermediate_files/votus_cov75thres.txt"), header=FALSE, comment.char="#")
colnames(votus_cov75thres) = c("sample","votus","coverage","meandepth","rpkm")


# Create new data frame with sample, votu, compartment, and phenotype
df = votus_cov75thres[,c(1,2)]
df["Compartment"] = rep("", nrow(df))
df["Phenotype"] = rep("", nrow(df))

# Match compartment and phenotype to each votu based on metadata
for (i in 1:nrow(df)){
  df[i,3] = metadata[metadata$Sample_ID == df[i,1], 4]
  df[i,4] = metadata[metadata$Sample_ID == df[i,1],5]
}

# Get phenotype counts


newdf = data.frame(matrix(ncol = 3, nrow = 0))
colnames(newdf) = c("votus", "Compartment", "Phenotype")


for (i in unique(df$votus)){ # For each vOTU
  # Convert rows with two phenotypes present to label both
  votu_subset = df[(df$votus == i),c(2,3,4)]
  compartments = unique(votu_subset$Compartment)
  for (j in compartments){ # For each compartment, find out if vOTU is in one phenotype or both
    votu_subset_comp = unique(votu_subset[(votu_subset$Compartment == j),])
    # Change label to both if in both, keep only one row if in one phenotype
    if (("Tall" %in% votu_subset_comp$Phenotype) && ("Short" %in% votu_subset_comp$Phenotype)){ # This works
      phen_label = "Both"
    } else {
      phen_label = unique(votu_subset_comp$Phenotype)
    }
    votu_subset_comp$Phenotype = phen_label
    votu_subset_comp = unique(votu_subset_comp)
    newdf = rbind(newdf, votu_subset_comp) # We get a new df that has replaced rows with both phenotypes to one row labeled "Both"
  }
}

# Now do the same for compartments
finaldf = data.frame(matrix(ncol = 3, nrow = 0))
colnames(finaldf) = c("votus", "Compartment", "Phenotype")


for (i in unique(df$votus)){ 
  votu_subset = newdf[(newdf$votus == i),]
  phenotypes = unique(votu_subset$Phenotype)
  for (j in phenotypes){
    votu_subset_phen = unique(votu_subset[(votu_subset$Phenotype == j),])
    
    if (("Bulk sediment" %in% votu_subset_phen$Compartment) && ("Rhizosphere" %in% votu_subset_phen$Compartment) && ("Endosphere" %in% votu_subset_phen$Compartment)){
      comp_label = "BRE"
    } else if (("Bulk sediment" %in% votu_subset_phen$Compartment) && ("Rhizosphere" %in% votu_subset_phen$Compartment)) {
      comp_label = "BR"
    } else if (("Rhizosphere" %in% votu_subset_phen$Compartment) && ("Endosphere" %in% votu_subset_phen$Compartment)){
      comp_label = "RE"
    } else if (("Bulk sediment" %in% votu_subset_phen$Compartment) && ("Endosphere" %in% votu_subset_phen$Compartment)){
      comp_label = "BE"
    } else {
      comp_label = unique(votu_subset_phen$Compartment)
    }
    votu_subset_phen$Compartment = comp_label
    votu_subset_phen = unique(votu_subset_phen)
    finaldf = rbind(finaldf, votu_subset_phen)
  }
}

# Matrix looks like this
#      Tall Short Both
# Bulk  0   0     0
# Rhizo 0   0     0
# Endo  0   0     0
# BR    0   0     0
# BE    0   0     0
# RE    0   0     0
# BRE   0   0     0

m = matrix(0, nrow = 7, ncol = 3)

# Get individual counts
for (i in 1:nrow(finaldf)){
  if (finaldf[i,2] == "Bulk sediment"){
    row_index = 1
  } else if (finaldf[i,2] == "Rhizosphere"){
    row_index = 2
  } else if (finaldf[i,2] == "Endosphere"){
    row_index = 3
  } else if (finaldf[i,2] == "BR") {
    row_index = 4
  } else if (finaldf[i,2] == "BE") {
    row_index = 5
  } else if (finaldf[i,2] == "RE") {
    row_index = 6
  } else {
    row_index = 7
  }
  
  if (finaldf[i,3] == "Tall") {
    col_index = 1
  } else if (finaldf[i,3] == "Short") {
    col_index = 2
  } else {
    col_index = 3
  }
  m[row_index,col_index] = m[row_index,col_index] + 1
}


#      Tall Short Both
#Bulk   15   59    5
#Endo   59  230   22
#Rhizo  20   58   14
#BR     13   41    5
#BE     0    1    0
#RE     23  205   14
#BRE    1   22    2

# Totals won't add up to 769 because some vOTUs are repeated

# Next calculate totals
observations = data.frame(m)
observations = rbind(observations, colSums(observations))
observations = cbind(observations, X4 = rowSums(observations))

# Keep only totals
expectations = data.frame(matrix(0,nrow = 8, ncol = 4)) # Initialize empty data frame
expectations$X4 = observations$X4 # Get row sums
expectations[8,] = observations[8,] # Get col sums

# Calculate expected values
rownum = nrow(expectations)
colnum = ncol(expectations)

for (i in 1:(rownum-1)) {
  for (j in 1:(colnum-1)){
    expectations[i,j] = round((expectations[i,colnum])*(expectations[rownum, j])/(expectations[rownum,colnum]))
  }
}

# Q for Marian - Should I round since they aren't whole values or keep for statistical tests

rows = c("Bulk sediment", "Endosphere","Rhizosphere","Bulk/Rhizo","Bulk/Endo","Rhizo/Endo","Bulk/Endo/Rhizo","Total")
cols = c("Tall", "Short", "Both", "Total")

rownames(observations) = rows
colnames(observations) = cols
rownames(expectations) = rows
colnames(expectations) = cols

### CHI SQUARED TEST
  # Q - why do we say none are present in root/bulk in figure when there is one?
  # Q - should I still keep compartment separated in the 7 categories for this, or the 3?

# Organize values into expected/observed table
melted_obs = melt(as.matrix(observations[1:7,1:3]))
melted_exp = melt(as.matrix(expectations[1:7,1:3]))
my_table = cbind(melted_obs, melted_exp$value)
colnames(my_table) = c("Compartment", "Phenotype", "Observed", "Expected")

# Things we need to test
# 1. Proportion of Short vs tall for those found in one phenotype (keep only rows where phenotype is not both)
# 2. Proportion of bulk, root, rhizo for those found in one compartment (keep only rows where compartment is singular)
# 3. Proportion of bulk/rhizo and root/rhizo for those found in two compartments (keep only rows where compartment is double)

# 1
sub_table = my_table[(my_table$Phenotype != "Both"),]
chisq.test(sub_table$Observed, sub_table$Expected)
# Result: X-squared = 154, df = 143, p-value = 0.2503

# 2
sub_table = rbind(my_table[(my_table$Compartment == "Bulk sediment"),], my_table[(my_table$Compartment == "Rhizosphere"),], my_table[(my_table$Compartment == "Endosphere"),])
chisq.test(sub_table$Observed, sub_table$Expected)
# Result: X-squared = 63, df = 56, p-value = 0.2425

# 3
sub_table = rbind(my_table[(my_table$Compartment == "Bulk/Rhizo"),], my_table[(my_table$Compartment == "Bulk/Endo"),], my_table[(my_table$Compartment == "Rhizo/Endo"),])
chisq.test(sub_table$Observed, sub_table$Expected)
# Result: X-squared = 63, df = 49, p-value = 0.08625


# Warning for all of them: Chi-squared approximation may be incorrect


