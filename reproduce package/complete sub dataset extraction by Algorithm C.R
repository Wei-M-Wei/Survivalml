#######
We need package 'Survivalml' including the empirical dataset to run this .r file, please jump to 'complete sub dataset extraction by Algorithm C (import raw dataset).R'.
#######


rm(list = ls())
library(readxl)
library(lubridate)
library(dplyr)
library(stringr)
library(openxlsx)
library(Survivalml)

#
data("company_data_raw")

dim(company_data_raw)

data_process_all = company_data_raw
dim(data_process_all) # 5363 8898

############change '--' to NA
data_process_all[data_process_all == "--"] <- NA

###############extract the listing date
start_day = as.Date(data_process_all$`First listing date`)
start_year = year(start_day)

################select firms whose IPO date is from 1985 to 2015
index = which(start_year >= 1985 & start_year <= 2015)
data_process = data_process_all[index,]

#################delete firms from financial sector
index_fin = which(data_process[,6]=='J')
data_process = data_process[-index_fin,]
dim(data_process)
################extract the ST day
bankrupt_day = ifelse(data_process$`risk time` == "--", 1, str_extract(data_process$`risk time`, "(?<=:)[^;]+"))
bankrupt_day = as.Date(bankrupt_day, format = '%Y%m%d')
data_process$bankrupt_day = bankrupt_day
data_process$start_day = data_process$`First listing date`
######################check if the firm is censored in our observation period
end_observation = as.Date('2020-12-31')
data_process$status = 0
for (i in seq_along(data_process$bankrupt_day)) {
  if ( is.na(data_process$bankrupt_day[i])  || data_process$bankrupt_day[i] >= end_observation ) {
    data_process$status[i] = 0
    data_process$bankrupt_day[i] = end_observation
  }
  else{
    data_process$status[i] = 1
  }
}

sum(data_process$status)
dim(data_process)

###################calculate the survival time T
days_difference <- as.numeric( data_process$bankrupt_day - as.Date(data_process$start_day) )

# change the unit of difference as year
years_difference <- days_difference / 365
data_process$survival_time = as.numeric(years_difference)
data_process$censoringtime = as.numeric(as.numeric( as.Date(end_observation) - as.Date(data_process$start_day) )/365)
ind = which(data_process$status == 1)
length(ind)
mean(data_process$censoringtime[ind])
mean(data_process$censoringtime[-ind])
sort(data_process$survival_time)
data_process[order(data_process$survival_time),]$start_day
data_process[order(data_process$survival_time),]$bankrupt_day
########################################################
data_before = data_process
id = which(data_process[,6]=='C')
data_before = data_process[id,]
dim(data_before)
sum(data_before$status)
#write.csv2(data_before, 'data_before.csv', row.names = TRUE,)
######################based on s, extract sub-dataset from the original dataset
# if let s =10, you should change line 127 as lag_use_year = 4, which is described in the empirical section of this paper
dim(data_before)
for (s in c(6)){
  print(s)
  data_before = data_process
  # manufacturing sector
  id = which(data_process[,6]=='C')
  data_before = data_process[id,]
  id_s = which(data_before$survival_time >= s)
  data_pre = data_before[id_s,]
  n = dim(data_pre)[1]
  end_in = 2023
  sat_in = 1985

  ################################
  operation_num = 6 # 456 / l # 6
  index_operation = seq( 7 , 7  + ((end_in - sat_in) + 1)*4*operation_num - 1, 1 )
  total_operation = length(index_operation)
  ###################################################
  debt_num = 10 # 912 / l # 12
  index_debt = seq( 7 , 7  + ((end_in - sat_in) + 1)*4*debt_num - 1, 1 ) + total_operation
  total_debt = length(index_debt)
  #################################################################
  profit_num = 16# 684 / l # 9
  index_profit = seq( 7, 7  + ((end_in - sat_in) + 1)*4*profit_num - 1, 1 ) + total_operation + total_debt
  total_profit = length(index_profit)
  ##############################
  potential_num = 6 # 988 / l # 13
  index_potential = seq( 7, 7  + ((end_in - sat_in) + 1)*4*potential_num - 1, 1 ) + total_operation + total_debt + total_profit
  total_potential = length(index_potential)
  #############################################################
  zscore_num = 5 # 456 / l # 6
  index_zscore = seq( 7, 7  + ((end_in - sat_in) + 1)*4*zscore_num - 1, 1 ) + total_operation + total_debt + total_profit + total_potential
  total_zscore = length(index_zscore)
  #############################################################
  captial_num = 6
  index_captial = seq( 7, 7  + ((end_in - sat_in) + 1)*4*captial_num - 1, 1 ) + total_operation + total_debt + total_profit + total_potential + total_zscore
  total_captial = length(index_captial)
  #############################################################
  stock_num = 5 # 6
  index_stock = seq( 7, 7  + ((end_in - sat_in) + 1)*4*stock_num - 1, 1 ) + total_operation + total_debt + total_profit + total_potential + total_zscore + total_captial
  total_stock = length(index_stock)
  #############################################################
  cash_num = 3 # 6
  index_cash = seq( 7, 7  + ((end_in - sat_in) + 1)*4*cash_num - 1, 1 ) + total_operation + total_debt + total_profit + total_potential + total_zscore + total_captial + total_stock
  total_cash = length(index_cash)
  ##############################################
  macro_num = 98
  index_macro = seq( 7, 7  + ((end_in - sat_in) + 1)*4*macro_num - 1, 1 ) + total_operation + total_debt + total_profit + total_potential + total_zscore + total_captial + total_stock + total_cash
  total_macro = length(index_macro)
  ######################################################
  lag_use_year = 6
  lags = lag_use_year * 4
  X_g = data.frame(matrix(ncol = 0, nrow = 0))

  for (i in seq(n)) {
    print(i)
    nex_year = next_quarter(year(data_pre$start_day[i]), month(data_pre$start_day[i]))$year
    nex_qu = next_quarter(year(data_pre$start_day[i]), month(data_pre$start_day[i]))$quarter
    start_index = (nex_year - sat_in ) * 4 + nex_qu + 4 * (s-lag_use_year)
    if (i == 1){
      X_g = cbind(
        data_extract(data_pre[i, index_operation], start_index, start_index + lags - 1, operation_num),
        data_extract(data_pre[i, index_debt], start_index, start_index + lags - 1, debt_num),
        data_extract(data_pre[i, index_profit], start_index, start_index + lags - 1, profit_num),
        data_extract(data_pre[i, index_potential], start_index, start_index + lags - 1, potential_num),
        data_extract(data_pre[i, index_zscore], start_index, start_index + lags - 1, zscore_num),
        data_extract(data_pre[i, index_captial], start_index, start_index + lags - 1, captial_num),
        data_extract(data_pre[i, index_stock], start_index, start_index + lags - 1, stock_num),
        data_extract(data_pre[i, index_cash], start_index, start_index + lags - 1, cash_num))
      #data_extract(data_pre[i, index_macro], start_index, start_index + lags - 1, macro_num))
    }
    else{
      X_g = rbind(as.matrix(X_g), as.matrix(cbind(
        data_extract(data_pre[i, index_operation], start_index, start_index + lags - 1, operation_num),
        data_extract(data_pre[i, index_debt], start_index, start_index + lags - 1, debt_num),
        data_extract(data_pre[i, index_profit], start_index, start_index + lags - 1, profit_num),
        data_extract(data_pre[i, index_potential], start_index, start_index + lags - 1, potential_num),
        data_extract(data_pre[i, index_zscore], start_index, start_index + lags - 1, zscore_num),
        data_extract(data_pre[i, index_captial], start_index, start_index + lags - 1, captial_num),
        data_extract(data_pre[i, index_stock], start_index, start_index + lags - 1, stock_num),
        data_extract(data_pre[i, index_cash], start_index, start_index + lags - 1, cash_num))
        #,data_extract(data_pre[i, index_macro], start_index, start_index + lags - 1, macro_num))
      ))
    }
  }

  X = data.frame(X_g)
  X$start_day = data_pre$start_day
  X$bankrupt_day = data_pre$bankrupt_day
  X$survival_time = data_pre$survival_time
  X$status = data_pre$status
  X$censoringtime = data_pre$censoringtime


  # Algorithm C
  pre_row = seq(25, dim(X)[1], 50)
  pre_col = seq(25, dim(X)[2], 50)
  judge1 = matrix(0, length(pre_row), length(pre_col))
  judge2 = matrix(0, length(pre_row), length(pre_col))

  for (i in seq(length(pre_row))){
    for ( j in seq(length(pre_col))){
      print(c(i,j))
      X2 <- X
      # delete firms whose missing value is larger than the threshold
      er = which(rowSums(is.na(X2)) > pre_col[j])
      if (length(er) == 0){
        X2 = X2
      }else{
        X2 = X2[-er,]
      }
      dim(X2)
      p_lags = lags
      index_uc = which(X[-er,]$status == 1)


      # delete variables whose missing value is larger than the threshold
      empty_columns <- unique(which(colSums(is.na(X2)) > pre_row[i]))
      ind_remove = NULL
      for (k in empty_columns) {
        i_index = k / (p_lags)
        if (i_index  %% 1 == 0) {
          ind_remove = c(ind_remove, i_index)
        } else {
          ind_remove = c(ind_remove, floor(i_index) + 1)
        }
      }
      ind_remove
      remove <- unique(ind_remove)
      remove_all = NULL
      for (h in remove) {
        id = seq((h-1)*p_lags + 1, h*p_lags)
        remove_all = c(remove_all, id)
      }
      if (length(er) == 0 || length(remove_all) == 0){
        X_design = X
      }else{
        X_design = X[-er , -remove_all]
      }

      ##########to get the final complete data, delete firms which still have missing values
      empty_rows = which(rowSums(is.na( X_design )) > 0)
      X_final = X_design[-empty_rows,]
      dim(X_final)

      # store the censoring rate of this subdataset.
      if (dim(X_final)[1] == 0){
        judge1[i, j] = 0
      }else{
        judge1[i, j] = sum(X_final$status)
      }

      # control the number of covariates
      if( (dim(X_final)[1] / dim(X)[1]) >= 0.5 && (dim(X_final)[2] / dim(X)[2]) >= 0.5){ #&& (dim(X_final)[2] / dim(X_final)[1]) >= 0.5
        judge2[i, j] = sum(X_final$status)
      }else{
        judge2[i, j] = 0
      }
    }
  }

  which(judge1 == max(judge1), arr.ind = TRUE)
  which(judge2 == max(judge2), arr.ind = TRUE)
  X2 <- X
  er = which(rowSums(is.na(X2)) > pre_col[ which(judge2 == max(judge2), arr.ind = TRUE)[1,][2]])
  if (length(er) == 0){
    X2 = X
  }else{
    X2 = X[-er,]
  }
  dim(X2)
  p_lags = lags

  empty_columns <- unique(which(colSums(is.na(X2)) > pre_row[ which(judge2 == max(judge2), arr.ind = TRUE)[1,][1]] ))

  ind_remove = NULL
  for (k in empty_columns) {
    i_index = k / (p_lags)

    if (i_index  %% 1 == 0) {
      ind_remove = c(ind_remove, i_index)
    } else {
      ind_remove = c(ind_remove, floor(i_index) + 1)
    }
  }
  ind_remove
  remove <- unique(ind_remove)
  remove_all = NULL
  for (h in remove) {
    id = seq((h-1)*p_lags + 1, h*p_lags)
    remove_all = c(remove_all, id)
  }

  if (length(er) == 0 || length(remove_all) == 0){
    X_design = X
  }else{
    X_design = X[-er , -remove_all]
  }
  empty_rows = which(rowSums(is.na( X_design )) > 0)
  X_final = X_design[-empty_rows,]
  dim(X_final)

  ########################
  ####################################################### remove the last lag
  quarter = 4
  lag_of_covariates = (lag_use_year-0.25)
  interval <- lag_of_covariates*quarter

  # Get the total number of columns in the data frame
  n_cols <- ncol(X_final)

  # Create a vector to store the selected columns
  selected_columns <- c()

  # Generate the column indices with the specified interval
  for (i in seq(1, n_cols, by = interval + 1)) {
    selected_columns <- c(selected_columns, i:(min(i + interval - 1, n_cols)))
  }
  X_final = X_final[, selected_columns]
  dim(X_final)
  sum(X_final$status)/dim(X_final)[1]
  min(X_final$censoringtime)
  table_name <- paste0("period_final_new", s, ".xlsx")
  write.xlsx(X_final, table_name)
}

