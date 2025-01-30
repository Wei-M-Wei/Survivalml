library(reshape2)
library(ggplot2)
s = 10
quarter = 4
degree = 2
#################################
data = read.csv2('period_final_new10.csv')
p  = dim(data)[2] - 5
name_v = c()
for (j in seq(1, p, s*quarter)){
  name_v = c(name_v, rep(names(data)[j],3))
}
######################################
t = 13
cof = read.csv2('MIDAS10_cof_AUC13.csv')
cof
cof_e = c()
for (i in seq(2,dim(cof)[1])){
  split_string <- strsplit(cof[i,1], ",")[[1]]
  number_after_comma <- as.numeric(split_string[2])
  cof_e = c(cof_e, number_after_comma)
}

names(cof_e) = name_v

cof_analysis = matrix(0, length(cof_e[seq(1,length(cof_e),degree+1)]), 3)
cof_analysis[which(cof_e[seq(1,length(cof_e),degree+1)] !=0 ),1] = as.numeric(1)

#######################################################
t = 13.5
cof = read.csv2('MIDAS10_cof_AUC13.5.csv')
cof
cof_e = c()
for (i in seq(2,dim(cof)[1])){
  split_string <- strsplit(cof[i,1], ",")[[1]]
  number_after_comma <- as.numeric(split_string[2])
  cof_e = c(cof_e, number_after_comma)
}

names(cof_e) = name_v

cof_analysis[which(cof_e[seq(1,length(cof_e),degree+1)] !=0 ),2] = as.numeric(1)

##################################################

t = 14
cof = read.csv2('MIDAS10_cof_AUC14.csv')
cof
cof_e = c()
for (i in seq(2,dim(cof)[1])){
  split_string <- strsplit(cof[i,1], ",")[[1]]
  number_after_comma <- as.numeric(split_string[2])
  cof_e = c(cof_e, number_after_comma)
}

names(cof_e) = name_v
cof_analysis[which(cof_e[seq(1,length(cof_e),degree+1)] !=0 ),3] = as.numeric(1)
############################################
names_pic <- sub("\\.\\.\\..*", "", names(cof_e[seq(1,length(cof_e),degree+1)]))
names_pic <- gsub("\\.", " ", names_pic)

#############

col_labels <- c('13 years','13.5 years','14 years')
names_pic = c( "Current assets turnover ratio"
               , "Fixed assets turnover rate"
               ,"Total assets turnover ratio"

               , "Current ratio"
               ,"Quick ratio"
               ,"Equity ratio"
               ,"Total tangible assetss / Total liabilities"
               ,"Tangible assetss / Net debt"
               ,"Earnings before interest, tax, depreciation, and amortization / Total liabilities"

               ,"Return on total assetss (ROA)"
               ,"Net profit on assetss"
               ,"Net profit / Total operating income"
               ,"Earnings before interest and taxes / Total operating income"
               , "Total operating income (year on year growth rate)"
               ,"Gross profit (year on year growth rate)"
               ,"Operating profit (year on year growth rate)"
               ,"Total profit (year on year growth rate)"

               ,"Total debt (year on year growth rate)"
               ,"Total assetss (year on year growth rate)"

               ,"X1 Working capital / Total assets"
               ,"X2 Retained earnings / Total assets"
               ,"X3 Earnings before interest and taxes / Total assets"
               ,"X4 Market value of equity / Book value of total liabilities"
               ,"X5 Sales / Total assets"

               ,"Total shareholders equity / Total liabilities"
               ,"Debt ratio"
               ,"Interest bearing debt ratio"
               ,"Equity multiplier"
               ,"Current assetss / Total assetss"
               ,"Current liabilities / Total liabilities"

               ,"Net cash flow from operating activities per share"
               ,"Operating income per share"
               ,"Profit before tax per share"
               ,"Net assetss per share BPS"

               ,"Net cash flow from operating activities / Operating income"
               ,"Net operating Cash flow / Operating income")
# Assign labels to the matrix


rownames(cof_analysis) <- names_pic
colnames(cof_analysis) <- col_labels
melted_data <- melt(cof_analysis)
colnames(melted_data) <- c("Variables", "Time", "Value")
# Create heatmap with ggplot2
ggplot(melted_data, aes(Time, Variables, fill = factor(Value))) +  # Use factor for discrete values
  geom_tile() +  # Add border to tiles
  scale_fill_manual(values = c("0" = "white", "1" = "red")) +  # Custom colors
  labs(
       x = "Prediction horizon",
       y = "Variable") +
  theme_minimal() +
  theme(legend.position = "none")


####################################################################################
names_pic = c('Operation ability related', 'Debt related', 'Profit related', 'Potential related',
              'Z-score related', 'Captial related', 'Stock related', 'Cash related')

#######################################
cof_final = matrix(0, 8, 3)
cof_final[1,] = colSums(cof_analysis[1:3,])/6
cof_final[2,] = colSums(cof_analysis[4:9,])/10
cof_final[3,] = colSums(cof_analysis[10:17,])/16
cof_final[4,] = colSums(cof_analysis[18:19,])/6
cof_final[5,] = colSums(cof_analysis[20:24,])/5
cof_final[6,] = colSums(cof_analysis[25:30,])/6
cof_final[7,] = colSums(cof_analysis[31:34,])/5
cof_final[8,] = colSums(cof_analysis[35:36,])/3
##############################################

rownames(cof_final) <- names_pic
colnames(cof_final) <- col_labels
melted_data <- melt(cof_final)
colnames(melted_data) <- c("Variables", "Time", "Value")
#melted_data$Time = as.numeric(melted_data$Time)
#melted_data$Variables = as.numeric(melted_data$Variables)
melted_data$Value = as.numeric(melted_data$Value)
# Create heatmap with ggplot2
ggplot(melted_data, aes(Time, Variables, fill = Value)) +  # Use factor for discrete values
  geom_tile() +  # Add border to tiles
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +  # Custom colors
  labs(
       x = "Prediction horizon",
       y = "Variable",
       fill = "Value") +
  theme_minimal()



