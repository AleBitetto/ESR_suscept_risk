
memory.limit(size=100000000000)
library(haven)
# library(maps)
library(ggplot2)
library(nipals)
library(softImpute)
library(tensorBF)
library(factoextra)
library(rpca)
library(grid)
library(gtable)
library(sparsepca)
library(gridtext) # devtools::install_github("r-lib/xml2")   devtools::install_github("clauswilke/gridtext")
library(ggrepel)
library(gridExtra)
library(plot3D)
library(magick)
library(ggridges)
# library(psychNET)
library(MARSS)
library(sparsevar)
library(Matrix)
library(FKF)
library(dse)
library(tseries)
# library(EnvStats)
# library(parallelMap)
# library(parallel)
# library(party)
# library(randomForest)
# library(gbm)
# library(mlr)
library(ppcor)
require(Hmisc)
library(data.table)
library(dplyr)
library(tidyverse)

source('./Help_ESR.R')
source('./Help_mlr.R')

# compile functions
{
  library(compiler)
  enableJIT(3)
  setCompilerOptions(optimize=3)
  setCompilerOptions(suppressAll = TRUE)
  funlist=lsf.str()
  for (i in c(1:length(funlist))){
    comfun=cmpfun(eval(parse(text=funlist[i])),options=list(suppressUndefined=T))
    assign(funlist[i],comfun)
  }
}


# load data and evaluate perimeter - skip if reloading
force_remove_variable = c('rain', 'temp', 'var1', 'bcrisis', 'ehr1', 'ehr2', 'ehr3')  # totally full of missing for 2017 to 2019
{
  df_orig = read_dta("Data/esr_basic_data_2010-2019.dta") %>%
    select(-id, -country1, -iso2, -iso3, -income1, -income2, -region, -regionname)
  
  # extract description for each variable
  lab_table = c()
  for (col in colnames(df_orig)){
    att = attributes(df_orig %>% select(all_of(col)) %>% pull(col))$label
    if (!is.null(att)){
      lab_table = lab_table %>%
        bind_rows(data.frame(var = col, label = att, stringsAsFactors = F))
    }
  }
  write.table(lab_table, 'Variable_description.csv', sep = ';', row.names = F)
  
  # check missing percentage
  df_orig = df_orig %>%
    mutate_all(function(x) { attributes(x) <- NULL; x })
  
  summary_missing_by_country = df_orig %>%
    gather('variable', 'val', -c(year, country)) %>%
    group_by(country) %>%
    summarize(Tot_non_NA = sum(!is.na(val)),
              Tot_NA = sum(is.na(val)), .groups = 'drop') %>%
    mutate(NA_perc = round(Tot_NA / (Tot_NA + Tot_non_NA) * 100, 1)) %>%
    arrange(NA_perc)
  
  summary_missing_by_year = df_orig %>%
    gather('variable', 'val', -c(year, country)) %>%
    group_by(year) %>%
    summarize(Tot_non_NA = sum(!is.na(val)),
              Tot_NA = sum(is.na(val)), .groups = 'drop') %>%
    mutate(NA_perc = round(Tot_NA / (Tot_NA + Tot_non_NA) * 100, 1)) %>%
    arrange(year)
  
  summary_missing_by_variable = df_orig %>%
    gather('variable', 'val', -c(year, country)) %>%
    group_by(variable) %>%
    summarize(Tot_non_NA = sum(!is.na(val)),
              Tot_NA = sum(is.na(val)), .groups = 'drop') %>%
    mutate(NA_perc = round(Tot_NA / (Tot_NA + Tot_non_NA) * 100, 1)) %>%
    arrange(NA_perc)
  
  write.table(summary_missing_by_country, "./Stats/00_summary_missing_by_country.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  write.table(summary_missing_by_year, "./Stats/00_summary_missing_by_year.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  write.table(summary_missing_by_variable, "./Stats/00_summary_missing_by_variable.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  # select perimeter - first select years, then countries by perc and finally variables by perc
  years = c(2010:2019)
  country_max_perc = 40
  variable_max_perc = 45
  
  df_perim = df_orig %>%
    gather('variable', 'val', -c(year, country)) %>%
    filter(year %in% years)
  
  country_perim = df_perim %>%
    group_by(country) %>%
    summarize(Tot_non_NA = sum(!is.na(val)),
              Tot_NA = sum(is.na(val)), .groups = 'drop') %>%
    mutate(NA_perc = round(Tot_NA / (Tot_NA + Tot_non_NA) * 100, 1)) %>%
    filter(NA_perc <= country_max_perc) %>%
    pull(country)
  
  variable_perim = df_perim %>%
    filter(country %in% country_perim) %>%
    group_by(variable) %>%
    summarize(Tot_non_NA = sum(!is.na(val)),
              Tot_NA = sum(is.na(val)), .groups = 'drop') %>%
    mutate(NA_perc = round(Tot_NA / (Tot_NA + Tot_non_NA) * 100, 1)) %>%
    filter(NA_perc <= variable_max_perc) %>%
    pull(variable)
  
  df = df_orig %>%
    filter(year %in% years) %>%
    filter(country %in% country_perim) %>%
    select(country, year, all_of(variable_perim), -all_of(force_remove_variable)) %>%
    as.data.frame()
  
  # basic statistics
  write.table(basicStatistics(df), "./Stats/01_Statistics_pre_missing_imputation.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  saveRDS(df, './Checkpoints/df.rds')
  rm(lab_table, summary_missing_by_country, summary_missing_by_year, summary_missing_by_variable, df_perim)
}


df = readRDS('./Checkpoints/df.rds')

# select variable to be removed because of multicollinearity
corr_thresh = 0.4
{
  correlation_list = evaluate_correlation(df %>% select(-country, -year))
  var_to_remove = correlation_list %>%
    filter(abs > corr_thresh) %>%
    select(Var1, Var2) %>%
    unlist() %>%
    table() %>%
    as.data.frame(stringsAsFactors = F) %>%
    setNames(c('Variable', 'Occurrence')) %>%
    arrange(desc(Occurrence)) %>%
    filter(Occurrence >= 3) %>%
    pull(Variable) %>%
    c(correlation_list %>%
        filter(abs(Corr) > 0.7) %>%
        select(Var1, Var2) %>%
        unlist()) %>%
    unique()
  
  cat('\n--- full set variables:', ncol(df) - 2)
  cat('\n--- full set variables:', ncol(df) - 2 - length(var_to_remove), '\n\n')
}
restricted_df_var = c('growth', 'interest1', 'labor1', 'penn19', 'penn32', 'penn33', 'penn42',
                           'religdiv1', 'religdiv2', 'chrs_pct', 'musl_pct',
                           df %>% select(starts_with("var")) %>% colnames(),
                           df %>% select(starts_with("health")) %>% colnames(),
                           df %>% select(starts_with("hef")) %>% colnames(),
                           df %>% select(starts_with("pop")) %>% colnames())
restricted_df_var = setdiff(restricted_df_var, c('health3', 'pop1', 'pop2', 'pop3', 'var8'))


### recover missing data
dummy_var = c()#c('bcrisis')    # apply median
recov_method_set = c('SOFT_IMPUTE', 'TENSOR_BF', 'FULL_AVERAGE')#, 'COUNTRY_AVERAGE')#,'NIPALS')
df_set = c('Original', 'Restricted')  # "Restricted" is without multicollinearity
{
  # save variable description for Restricted dataset
  lab_table = read.csv("Variable_description.csv", sep=";", stringsAsFactors=FALSE) %>%
    filter(var %in% (df %>% select(country, year, all_of(restricted_df_var)) %>% colnames()))
  write.table(lab_table, 'Variable_description_Restricted_Dataset.csv', sep = ';', row.names = F)
  
  
  res_recov = c()
  year_match = data.frame(year = unique(df$year)) %>% arrange(year) %>% mutate(CODE = c(1:n()))
  sink(paste0('./Log/Recov_missing_log_', format(Sys.time(), "%Y-%m-%dT%H%M"), '.txt'), append = F, split = T, type = 'output')
  cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
  for (recov_met in recov_method_set){
    cat('\n\n\n           ------------########################    evaluating', recov_met, '  ########################------------\n\n')
    for (df_type in df_set){
      
      df_t = c()
      if (df_type == 'Original'){df_t = df}
      if (df_type == 'Restricted'){df_t = df %>% select(country, year, all_of(restricted_df_var))}#{df_t = df %>% select(-all_of(var_to_remove))}
      
      # recover missing year by year
      if (recov_met %in% c('SOFT_IMPUTE', 'NIPALS')){
        
        for (yr in year_match$year){
          
          slice = df_t %>% filter(year == yr) %>% arrange(country) %>% select(-year, -country) %>%
            `rownames<-`(sort(unique(df_t$country)))
          na_ind = is.na(slice)
          slice_recov = c()
          
          if (sum(na_ind) == 0){
            cat('\n --- No missing for year ', yr, '- data', df_type)
            slice_recov = slice
          } else {
            
            # NIPALS
            if (recov_met == 'NIPALS'){
              nip = try_nipals(slice)
              if (!is.null(nip$no_error) | !is.null(nip$warn_text)){
                slice_recov = slice
                slice_recov[na_ind] = nip$value$fitted[na_ind]
                if (!is.null(nip$warn_text)){cat('\n ---- Warning in NIPALS - year', yr, '- data', df_type, ':\n         ', nip$warn_text)}
              } else if(!is.null(nip$error_text)){
                slice_recov = slice
                cat('\n #### Error in NIPALS - year', yr, '- data', df_type, ':\n         ', nip$error_text)
              }
              err = sum(slice_recov == slice, na.rm = T) + sum(na_ind) - nrow(slice) * ncol(slice)
              if (err != 0){cat('\n #### Error on non missing elements in NIPALS - year', yr, '- data', df_type, ':', err)}
            }
            
            # SOFT IMPUTE - uses https://arxiv.org/pdf/1410.2596.pdf - https://cran.r-project.org/web/packages/softImpute/vignettes/softImpute.html
            if (recov_met == 'SOFT_IMPUTE'){
              sf = try_soft_imp(slice)
              if (!is.null(sf$no_error) | !is.null(sf$warn_text)){
                slice_recov = sf$value
                if (!is.null(sf$warn_text)){cat('\n ---- Warning in SOFT IMPUTE - year', yr, '- data', df_type, ':\n         ', sf$warn_text)}
              } else if(!is.null(sf$error_text)){
                slice_recov = slice
                cat('\n #### Error in SOFT IMPUTE - year', yr, '- data', df_type, ':\n         ', sf$error_text)
              }
              err = sum(slice_recov == slice, na.rm = T) + sum(na_ind) - nrow(slice) * ncol(slice)
              if (err != 0){cat('\n #### Error on non missing elements in SOFT IMPUTE - year', yr, '- data', df_type, ':', err)}
            }
            
            # save results
            res_recov = res_recov %>% rbind(
              as.data.frame(as.table(na_ind)) %>% setNames(c("Var1", "Var2", "NA_ind")) %>%
                left_join(as.data.frame(as.table(as.matrix(slice_recov))), by = c("Var1", "Var2")) %>%
                mutate(year = yr,
                       method = recov_met,
                       data = df_type)
            )
            
          } # else no missing
          
        } # yr
        
        # recover missing for all years together
        # TENSOR BF - https://www.biorxiv.org/content/biorxiv/early/2016/12/29/097048.full.pdf
      } else if (recov_met == 'TENSOR_BF'){
        
        # create tensor
        tens = array(numeric(),c(uniqueN(df_t$country), ncol(df_t) - 2, max(year_match$CODE)))
        tens_na = array(logical(),c(uniqueN(df_t$country), ncol(df_t) - 2, max(year_match$CODE)))
        for (i in year_match$CODE){
          t = df_t %>% filter(year == (year_match %>% filter(CODE == i))$year) %>% arrange(country) %>% select(-year, -country) %>%
            `rownames<-`(sort(unique(df_t$country))) %>% as.matrix()
          names_col = colnames(t)
          names_row = rownames(t)
          tens[, , i] = t
          tens_na[, , i] = is.na(tens[, , i])
        }
        
        # recover mising
        # tens_recov = c()
        # ten = try_tensorBF(tens, K = 10)
        # if (!is.null(ten$no_error) | !is.null(ten$warn_text)){
        #   tens_recov = ten$value[na_ind]
        #   if (!is.null(ten$warn_text)){cat('\n ---- Warning in TENSORBF - year', yr, '- data', df_type, ':\n         ', ten$warn_text)}
        # } else if(!is.null(ten$error_text)){
        #   tens_recov = tens
        #   cat('\n #### Error in TENSORBF - data', df_type, ':\n         ', ten$error_text)
        # }
        cat('\n\n      ---- data:', df_type, '\n\n')
        opts <- getDefaultOpts()
        opts$iter.burnin <- 5000
        set.seed((10))
        tbf = tensorBF(tens,
                       K = 15,
                       fiberCentering = 1,
                       slabScaling = 2,
                       noiseProp = c(0.5, 0.5),
                       opts = opts)
        tens_recov = predictTensorBF(tens, tbf)
        if (sum(is.na(tens_recov)) > 0){cat('\n\n\n #### Error TENSOR BF: remaining missing data:', sum(is.na(tens_recov)))}
        
        # reshape results
        for (i in year_match$CODE){
          res_recov = res_recov %>% rbind(
            as.data.frame(as.table(tens_na[, , i] %>% `colnames<-`(names_col) %>% `rownames<-`(names_row))) %>% setNames(c("Var1", "Var2", "NA_ind")) %>%
              left_join(as.data.frame(as.table(as.matrix(tens_recov[, , i]  %>% `colnames<-`(names_col) %>% `rownames<-`(names_row)))), by = c("Var1", "Var2")) %>%
              mutate(year = (year_match %>% filter(CODE == i))$year,
                     method = recov_met,
                     data = df_type)
          )
        }
        # average over all rows
      } else if (recov_met == 'FULL_AVERAGE'){
        recov_df = df_t %>%
          mutate_at(all_of(dummy_var), ~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
          mutate_at(all_of(setdiff(colnames(df_t), c('country', 'year', dummy_var))), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)) %>%
          ungroup() %>%
          gather('Var2', 'Freq', -c(year, country)) %>%
          rename(Var1 = country)
        
        na_ind = df_t %>%
          gather('variable', 'val', -c(year, country)) %>%
          mutate(NA_ind = is.na(val)) %>%
          rename(Var1 = country,
                 Var2 = variable) %>%
          select(Var1, Var2, year, NA_ind)

        if (sum(is.na(recov_df)) > 0){cat('\n\n\n #### Error FULL_AVERAGE: remaining missing data:', sum(is.na(recov_df)))}
        
        res_recov = res_recov %>% rbind(
          na_ind %>%
            left_join(recov_df, by = c("Var1", "Var2", "year")) %>%
            mutate(method = recov_met,
                   data = df_type)
        )
        # average by country
      } else if (recov_met == 'COUNTRY_AVERAGE'){
        recov_df = df_t %>%
          group_by(country) %>%
          mutate_at(all_of(dummy_var), ~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
          ungroup() %>%
          group_by(country) %>%
          mutate_at(all_of(setdiff(colnames(df_t), c('country', 'year', dummy_var))), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)) %>%
          ungroup() %>%
          gather('Var2', 'Freq', -c(year, country)) %>%
          rename(Var1 = country)
        
        na_ind = df_t %>%
          gather('variable', 'val', -c(year, country)) %>%
          mutate(NA_ind = is.na(val)) %>%
          rename(Var1 = country,
                 Var2 = variable) %>%
          select(Var1, Var2, NA_ind)
        
        if (sum(is.na(recov_df)) > 0){cat('\n\n\n #### Error COUNTRY_AVERAGE: remaining missing data:', sum(is.na(recov_df)))}
        
        res_recov = res_recov %>% rbind(
          na_ind %>%
            left_join(recov_df, by = c("Var1", "Var2")) %>%
            mutate(method = recov_met,
                   data = df_type)
        )
      }
      
    } # df_type
  } # recov_met
  cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
  sink()
  
  colnames(res_recov)[c(1:2, 4)] = c('country', 'variable', 'val')
  res_recov$country = as.character(res_recov$country)
  res_recov$variable = as.character(res_recov$variable)
  
  # check residual missing and impute year mean for country+variable
  residual_missing = sum(is.na(res_recov))
  if (residual_missing > 0){
    variable_missing = res_recov %>%
      filter(is.na(val)) %>%
      pull(variable) %>%
      unique()
    missing_report = res_recov %>%
      filter(variable %in% variable_missing) %>%
      group_by(variable, data, method, year) %>%
      summarize(Total_missing = sum(is.na(val)),
                Total_not_missing = sum(!is.na(val)),
                Total_country = uniqueN(country), .groups = 'drop') %>%
      filter(Total_missing > 0) %>%
      as.data.frame()
    cat('\n\n--- Total missing still present:', residual_missing, '\n')
    print(missing_report)
    if (sum(missing_report$Total_missing) != residual_missing){cat('\n\n############## residual missing calculation mismatch\n')}
    
    cat('\n\n--- Replacing with annual average (by variable+country)\n')
    res_recov = res_recov %>%
      group_by(variable, country, data, method) %>%
      mutate_at('val', ~ifelse(is.na(.), mean(., na.rm = TRUE), .))
    if (sum(is.na(res_recov)) > 0){cat('\n\n############## missing in res_recov still present\n')}
  }
  saveRDS(res_recov, './Checkpoints/res_recov.rds')
  
} # skip if reloading
res_recov = readRDS('./Checkpoints/res_recov.rds')



### differentiate time series for stationarity
dummy_var = c()#c('bcrisis')    # don't differentiate
evaluate_stationarity_test = F
{
  res_recov_added_extended = res_recov
  
  total_new_rows = res_recov_added_extended %>%
    group_by(data, method, country, variable) %>%
    summarise(TOT_YEAR = n(), .groups = 'drop') %>%
    mutate(TOT_NEW_YEAR = TOT_YEAR - 1) %>%
    pull(TOT_NEW_YEAR) %>%
    sum()
  res_stationarity = c()
  column_order = res_recov_added_extended %>% select(-NA_ind) %>% colnames()
  res_recov_added_diff = matrix(NA, ncol = length(column_order), nrow = total_new_rows)
  start = 1
  end = 0
  for (df_type in unique(res_recov_added_extended$data)){
    t1 = res_recov_added_extended %>%
      filter(data == df_type)
    # for (recov_met in unique(res_recov_added_extended$method)){
    for (recov_met in unique(t1$method)){
      t2 = t1 %>%
        filter(method == recov_met)
      # for (cc in unique(res_recov_added_extended$country)){
      for (cc in unique(t2$country)){
        t3 = t2 %>%
          filter(country == cc)
        # for (var in unique(res_recov_added_extended$variable)){
        for (var in unique(t3$variable)){
          
          new_lab = ifelse(df_type == 'Original', 'Difference', 'RestrictedDifference')
          
          # tt = res_recov_added_extended %>%
          #   filter(data == df_type) %>%
          #   filter(method == recov_met) %>%
          #   filter(country == cc) %>%
          #   filter(variable == var) %>%
          #   mutate(data = new_lab) %>%
          #   arrange(year)
          tt = t3 %>%
            filter(variable == var) %>%
            mutate(data = new_lab) %>%
            arrange(year)
          
          
          if (nrow(tt) > 0){
            
            if (startsWith(var, 'GEO') | (var %in% dummy_var)){
              tt_to_add = tt %>%
                select(-NA_ind) %>%
                filter(row_number() != 1)
            } else {
              # test Augmented Dickey-Fuller (small p-value means stationarity)
              # test Ljung-Box (high p-value means stationarity -> 1-p-val is saved)
              ts = tt$val
              ts_diff = diff(tt$val)
              if (evaluate_stationarity_test){
                if (range(ts)[1] == range(ts)[2]){
                  ADF_before = ADF_after = LB_before = LB_after = 0
                } else {
                  ADF_before = suppressWarnings(adf.test(ts, k = 1)$p.value)
                  ADF_after = suppressWarnings(adf.test(ts_diff, k = 1)$p.value)
                  LB_before = 1 - Box.test(ts, lag = 1, type="Ljung-Box")$p.value
                  LB_after = 1 - Box.test(ts_diff, lag = 1, type="Ljung-Box")$p.value
                }
                res_stationarity = res_stationarity %>%
                  bind_rows(tt %>%
                              select(-NA_ind, -val, -year) %>%
                              filter(row_number() == 1) %>%
                              mutate(ADF_before = ADF_before,
                                     ADF_after = ADF_after,
                                     LB_before = LB_before,
                                     LB_after = LB_after) %>%
                              mutate(Legend = 'Low means stationarity')
                  )
              }
              
              tt_to_add = tt %>%
                select(-NA_ind, -val) %>%
                filter(row_number() != 1) %>%
                mutate(val = ts_diff)
            }
            
            tt_to_add = tt_to_add %>%
              select(all_of(column_order)) %>%
              as.matrix()
            
            end = start + nrow(tt_to_add) - 1
            # cat('\n', start, ' - ', end, '               ', nrow(tt_to_add))
            res_recov_added_diff[start:end, ] = tt_to_add
            start = end + 1
            
          } # nrow(tt) > 0
        } # var
      } # cc
    } # df_type
  } # recov_met
  res_recov_added_diff = res_recov_added_diff %>%
    as.data.frame(stringsAsFactors = F) %>%
    setNames(column_order) %>%
    mutate(val = as.numeric(val),
           year = as.numeric(year))
  if (sum(is.na(res_recov_added_diff)) > 0){cat('\n\n############### missing in res_recov_added_diff')}
  
  tot_year_recov = res_recov %>% group_by(data, method, country, variable) %>% summarize(COUNT = n(), .groups = 'drop') %>% pull(COUNT) %>% unique()
  tot_year_recov_diff = res_recov_added_diff %>% group_by(data, method, country, variable) %>% summarize(COUNT = n(), .groups = 'drop') %>% pull(COUNT) %>% unique()
  
  if (length(tot_year_recov) != 1){cat('\n\n############### different number of years in res_recov:', tot_year_recov)}
  if (length(tot_year_recov_diff) != 1){cat('\n\n############### different number of years in res_recov_added_diff:', tot_year_recov)}
  if (tot_year_recov_diff != tot_year_recov - 1){cat('\n\n############### dimensions after differentiation do not match')}
  
  write.table(res_stationarity, "./Stats/1_Stats_stationarity_test.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  saveRDS(res_recov_added_diff, './Checkpoints/res_recov_added_diff.rds')
  rm(t, t1, t2, t3)
}
res_recov_added_diff = readRDS('./Checkpoints/res_recov_added_diff.rds')



### merge dataset
dummy_var = c()#c('bcrisis')    # check
{
  df_final = res_recov_added_diff %>% 
    bind_rows(res_recov_added_extended %>% select(-NA_ind))
  
  dummy_check = df_final %>%
    filter(variable %in% dummy_var) %>%
    group_by(variable, method, data) %>%
    summarize(unique_vals = uniqueN(val),
              values = paste0(unique(val), collapse = ","), .groups = 'drop')
  dummy_error = dummy_check %>%
    filter(unique_vals != 2)
  
  if (nrow(dummy_error) > 0){
    cat('\n\n#### dummy variable with more than 2 values:\n\n')
    print(dummy_error %>% as.data.frame())
  }
  
  saveRDS(df_final, './Checkpoints/df_final.rds')
}
df_final = readRDS('./Checkpoints/df_final.rds')



### PCA and robust PCA
cv_eval = F # to Cross-Validate PCA for suggested number of PC
k_fold = 3
pca_met_set = c('PCA', 'RobPCA')#, 'RobSparPCA')
dummy_var = c()#c('bcrisis')    # don't scale
{
  year_set = c(sort(unique(df_final$year)))#, 'Concat')
  res_PCA_list = list()
  res_PCA_importance = res_PCA_loadings = max_variance_PCA_list = sparseness_count = c()
  res_PCA_concat = list(max_variance_list = c(), max_variance_list = c())
  sink(paste0('./Log/PCA_fitting_log_', format(Sys.time(), "%Y-%m-%dT%H%M"), '.txt'), append = F, split = T, type = 'output')
  cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
  for (pc_met in pca_met_set){
    cat('\n\n\n           ------------########################    evaluating', pc_met, '  ########################------------\n\n')
    for (df_type in setdiff(unique(df_final$data), c('ExtendedOriginal', 'ExtendedDifference'))){
      for (recov_met in unique(df_final$method)){
        
        if (!df_type %in% c('ExtendedOriginal', 'ExtendedDifference')){
          year_set_t = year_set[!grepl('\\.', year_set)]   # exclude "fractional" years from interpolated series
        } else {
          year_set_t = year_set
        }
        if (grepl('Difference', df_type)){
          min_year=suppressWarnings(min(as.numeric(year_set_t), na.rm=T))
          year_set_t = year_set_t[!grepl(toString(min_year), year_set_t)]   # exclude first year for integrated series
        }
        
        for (yr in year_set_t){
          
          # concatenated slices (country x year_variable)
          if (yr == 'Concat'){
            
            slice = df_final %>%
              filter(data == df_type) %>%
              filter(method == recov_met) %>%
              filter(year != 'Avg') %>%
              mutate(variable = paste0(year, '_', variable)) %>%
              select(-data, -method, -year) %>%
              spread(variable, val)
            
            # single slice
          } else {
            slice = df_final %>%
              filter(data == df_type) %>%
              filter(method == recov_met) %>%
              filter(year == yr) %>%
              select(-data, -method, -year) %>%
              spread(variable, val)
          }
          rownames(slice) = slice$country
          slice = slice %>% select(-country)
          
          if (nrow(slice) == 0){
            cat('\n *****  skipped year:', yr, 'for data:', df_type, 'recov_met:', recov_met)
          } else {
            
            # standardize input
            slice_scaled = clever_scale(input_df = slice, exclude_var = dummy_var)

            # PCA
            if(pc_met == 'PCA'){
              res_pca = pca_fun(slice_scaled, k_fold, cv = cv_eval, method_title = paste0('PCA - ', recov_met, ' - ', yr))
              res_pca$pca$pca_input = slice_scaled
            }
            
            # Robust PCA - https://arxiv.org/abs/0912.3599
            if(pc_met == 'RobPCA'){
              robpca = rpca(as.matrix(slice_scaled), trace = F, max.iter = 10000)
              if (robpca$convergence$converged == F){cat('\n #### RobPCA did not converge - year:', yr, '- data:', df_type, '- method:', recov_met)}
              L = robpca$L; colnames(L) = colnames(slice_scaled); rownames(L) = rownames(slice_scaled)
              res_pca = pca_fun(L, k_fold, cv = cv_eval, method_title = paste0('RobPCA - ', recov_met, ' - ', yr))
              res_pca$pca$pca_input = L
              res_pca$pca$add_components = robpca$S   # used to reconstruct the original data for any number of pc - method specific
            }
            
            # Robust Sparse PCA - https://cran.r-project.org/web/packages/sparsepca/sparsepca.pdf - https://github.com/erichson/spca
            if(pc_met == 'RobSparPCA'){
              # tune sparsity parameter in spca_fun to increase loading sparsity
              res_pca = spca_fun(slice_scaled, k_fold, cv = cv_eval, method_title = paste0('RobSparPCA - ', recov_met, ' - ', yr))
              if (length(res_pca$conv_err) > 0){cat('\n #### Sparse RobPCA did not converge with', res_pca$conv_err, 'elements exceeding toll - year:', yr, '- data:', df_type, '- method:', recov_met)}
              sparseness_count = sparseness_count %>% rbind( cbind(res_pca$sparseness_count) %>%
                                                               mutate(year = yr,
                                                                      method = recov_met,
                                                                      data = df_type))
              res_pca$pca$pca_input = slice_scaled
              res_pca$pca$add_components = res_pca$sparse   # used to reconstruct the original data for any number of pc - method specific
            }
            
            # in res_pca$pca "x" is the score matrix - add also original matrix for later calculation of R^2 given the number of PC
            res_pca$pca$orig_data = slice
            res_pca$pca$orig_data_scaled = slice_scaled
            
            # save results
            if (yr == 'Concat'){
              res_PCA_concat[[df_type]][[recov_met]][[pc_met]] = res_pca
              res_PCA_concat[['max_variance_list']] = c(res_PCA_concat[['max_variance_list']], max(res_pca$importance_table$`Proportion of Variance`))
              res_PCA_concat[['res_loadings']] = res_PCA_concat[['res_loadings']] %>% rbind(
                cbind(PCA = pc_met, res_pca$load_table) %>%
                  mutate(year = yr,
                         method = recov_met,
                         data = df_type)
              )
            } else {
              res_PCA_list[[df_type]][[recov_met]][[toString(yr)]][[pc_met]] = res_pca
              
              max_variance_PCA_list = c(max_variance_PCA_list,max(res_pca$importance_table$`Proportion of Variance`))
              res_PCA_loadings = res_PCA_loadings %>% rbind(
                cbind(PCA = pc_met, res_pca$load_table) %>%
                  mutate(year = yr,
                         method = recov_met,
                         data = df_type)
              )
              res_PCA_importance = res_PCA_importance %>% rbind(
                cbind(PCA = pc_met, res_pca$importance_table, t(res_pca$pca$PC_opt), PC_opt = min(res_pca$pca$PC_opt)) %>%
                  mutate(year = yr,
                         method = recov_met,
                         data = df_type)
              )
            }
          } # year check
          
        } # yr
      } # recov_met
    } # df_type
  } # pc_met
  cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
  sink()
  
  # sparseness_stats = sparseness_count %>%
  #   mutate(PC = paste0('PC', formatC(PC, width = 2, format = "d", flag = "0"))) %>%
  #   spread(PC, ZERO_ELEM) %>%
  #   arrange(desc(data), method, year)
  # write.table(sparseness_stats, "./Stats/2_PCA_Sparseness_count.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  saveRDS(res_PCA_list, './Checkpoints/res_PCA_list.rds')
  saveRDS(res_PCA_loadings, './Checkpoints/res_PCA_loadings.rds')
  saveRDS(res_PCA_importance, './Checkpoints/res_PCA_importance.rds')
  saveRDS(max_variance_PCA_list, './Checkpoints/max_variance_PCA_list.rds')
  saveRDS(res_PCA_concat, './Checkpoints/res_PCA_concat.rds')
} # skip if reloading
res_PCA_list = readRDS('./Checkpoints/res_PCA_list.rds')
res_PCA_loadings = readRDS('./Checkpoints/res_PCA_loadings.rds')
res_PCA_importance = readRDS('./Checkpoints/res_PCA_importance.rds')
res_PCA_concat = readRDS('./Checkpoints/res_PCA_concat.rds')
max_variance_PCA_list = readRDS('./Checkpoints/max_variance_PCA_list.rds')


# save Explained Variance (also 95h and 99th percentile) for PCA
{
  res_PCA_stats=c()
  for (pc_met in names(res_PCA_list[[1]][[1]][[1]])){
    for (df_type in names(res_PCA_list)){
      for (recov_met in names(res_PCA_list[[1]])){
        for (yr in names(res_PCA_list[[1]][[1]])){
          
          ref_list=res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]]
          
          # take variance on full dataset
          cumul_var = ref_list$importance_table %>%
            select(PC, `Cumulative Proportion`) %>%
            mutate(PC = as.numeric(gsub('PC', '', PC))) %>%
            setNames(c('PC', 'Explain_var_Loadings'))
          
          # evaluate reconstructed X with different number of PC
          loadings = ref_list$pca$rotation
          X_pca = ref_list$pca$pca_input
          X_scaled = ref_list$pca$orig_data_scaled
          if (pc_met != 'PCA'){
            sparse = ref_list$pca$add_components
          } else {
            sparse = matrix(0, ncol = ncol(X_scaled), nrow = nrow(X_scaled))
          }
          
          expl_var_tab = c()
          for (pc in 1:ncol(loadings)){
            scores_pc = X_pca %*% loadings[, 1:pc] # you can also simply take first p columns of ref_list$pca$x
            X_reconstr = scores_pc %*% t(loadings[,1:pc]) + sparse
            
            TSS_val = (X_scaled - mean(X_scaled)) ^ 2
            TSS = sum(TSS_val)
            RSS_val = (X_scaled - X_reconstr) ^ 2
            RSS = sum(RSS_val)
            ind_95 = RSS_val <= quantile(RSS_val, 0.95)
            RSS_95 = sum(RSS_val[ind_95])
            TSS_95 = sum(TSS_val[ind_95])
            ind_99 = RSS_val <= quantile(RSS_val, 0.99)
            RSS_99 = sum(RSS_val[ind_99])
            TSS_99 = sum(TSS_val[ind_99])
            Explain_var = 1 - RSS / TSS
            Explain_var_95 = 1 - RSS_95 / TSS_95
            Explain_var_99 = 1 - RSS_99 / TSS_99
            
            expl_var_tab = expl_var_tab %>% rbind(
              c(PC = pc, Explain_var = Explain_var, Explain_var_99 = Explain_var_99, Explain_var_95 = Explain_var_95)
            )
          }
          expl_var_tab = expl_var_tab %>% as.data.frame(stringsAsFactors=F)
          cumul_var = cumul_var %>%
            left_join(expl_var_tab, by = "PC")
          
          
          res_PCA_stats = res_PCA_stats %>% bind_rows(
            cbind(data.frame(PCA=pc_met, method=recov_met, data=df_type, year=yr, cumul_var, stringsAsFactors=F))
          )
        } # yr
      } # recov_met
    } # df_type
  } # pc_met
  saveRDS(res_PCA_stats, './Checkpoints/res_PCA_stats.rds')
  write.table(res_PCA_stats, "./Stats/2_PCA_Stats_fitting.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
}
res_PCA_stats = readRDS('./Checkpoints/res_PCA_stats.rds')






PCA_to_keep = 'RobPCA'
df_type_to_keep = c("RestrictedDifference", "Restricted")


res_PCA_loadings = res_PCA_loadings %>%
  filter(PCA %in% PCA_to_keep) %>%
  filter(data %in% df_type_to_keep) %>%
  mutate(PCA = as.character(PCA))     # original column is factor
res_PCA_importance = res_PCA_importance %>%
  filter(PCA %in% PCA_to_keep) %>%
  filter(data %in% df_type_to_keep) %>%
  mutate(PCA = as.character(PCA))     # original column is factor



### plot all scree plot
ordered_df_type = c("Restricted", "RestrictedDifference") # "Original", "Difference", 
{
  # split by pc_met
  
  avail_years = c()
  for (l1 in res_PCA_list){
    for (l2 in l1){
      avail_years = c(avail_years, names(l2)) %>% unique() %>% sort()
    }
  }
  
  for (pc_met in names(res_PCA_list[[1]][[1]][[1]])){
    n_col = length(avail_years)     # number of years
    
    max_variance = c()
    for (l1 in res_PCA_list){
      for (l2 in l1){
        for (l3 in l2){
          max_variance = c(max_variance, l3[[pc_met]][['importance_table']] %>% pull(`Proportion of Variance`) %>% max()) %>% max()
        }
      }
    }
    
    row_list = list()
    for (yr in avail_years){
      i = 1
      for (recov_met in names(res_PCA_list[[1]])){
        for (df_type in ordered_df_type){
          
          if (!is.null(res_PCA_list[[df_type]][[recov_met]][[yr]])){   # check for empty year (for Difference df_type)
            row_list[[yr]][[i]] = ggplotGrob(
              res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]][['scree_plot']]+
                theme(axis.title.y=element_blank(),
                      axis.title.x=element_blank()) +
                ggtitle(paste0(recov_met, '\n', df_type, ' - ', yr)) +
                ylim(0, max_variance * 100 + 10)
            )
          } else {
            row_list[[yr]][[i]] = ggplotGrob(ggplot() + theme(panel.background = element_blank()))
          }
          i = i + 1
        } # df_type
      } # recov_met
    } # yr
    
    col_list = list()
    for (i in c(1:n_col)){
      col_list[[i]] = do.call(rbind, c(row_list[[i]], size="last"))
    }
    g = do.call(cbind, c(col_list, size="last"))
    png(paste0('./Stats/2_PCA_Scree_plot_', pc_met, '.png'), width = 30, height = 20, units = 'in', res=300)
    grid.draw(g)
    dev.off()
  } # pc_met
  
  # split by df_type
  
  # for (df_type in names(res_PCA_list)){
  #   n_row = length(names(res_PCA_list[[df_type]])) * length(names(res_PCA_list[[df_type]][[1]][[1]]))   # #_recov_met * #_PCA_meth
  #   n_col = length(names(res_PCA_list[[df_type]][[1]]))     # number of years
  #   
  #   row_list = list()
  #   for (yr in names(res_PCA_list[[df_type]][[1]])){
  #     i = 1
  #     for (recov_met in names(res_PCA_list[[df_type]])){
  #       
  #       for (pc_met in names(res_PCA_list[[df_type]][[1]][[1]])){
  #         row_list[[yr]][[i]] = ggplotGrob(
  #           res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]][['scree_plot']]+
  #             theme(axis.title.y=element_blank(),
  #                   axis.title.x=element_blank()) +
  #             ggtitle(paste0(pc_met, ' - ', recov_met, ' - ', yr)) +
  #             ylim(0, max(max_variance_PCA_list) * 100 + 10)
  #         )
  #         i = i + 1
  #       } # pc_met
  #     } # recov_met
  #   } # yr
  #   
  #   col_list = list()
  #   for (i in c(1:n_col)){
  #     col_list[[i]] = do.call(rbind, c(row_list[[i]], size="last"))
  #   }
  #   g = do.call(cbind, c(col_list, size="last"))
  #   png(paste0('./Stats/2_PCA_Scree_plot_', df_type, '.png'), width = 30, height = 20, units = 'in', res=300)
  #   grid.draw(g)
  #   dev.off()
  # } # df_type
}


### evaluate average (over years) Cumulative Explained Variance for df_type+pc_met+recov_met
PC_to_compare = 2
{
  summary_Exp_Var = res_PCA_importance %>%
    filter(PC == paste0("PC", PC_to_compare)) %>%
    group_by(PCA, data, method) %>%
    summarise(Average_Cum_Exp_Var = round(mean(`Cumulative Proportion`)*100, 1), .groups = "drop") %>%
    arrange(desc(Average_Cum_Exp_Var))
  write.table(summary_Exp_Var, "./Stats/2_PCA_Scree_plot_Avg_Cum_Exp_Var.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
}


### plot loadings comparison for each year
# for each year, all PC are inverted (change sign) such that leading variable [leading_var] has constant sign [leading_sign] over time
leading_var = "growth"
leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
signif_thresh_plot = 0.2  # significance threshold for loading (just for plotting lines)
PC_to_compare = 2
var_split_length = 9 # split variables into var_split_length chunks over multiple lines for better visualisation
# split by variables
{
  y_range = range((res_PCA_loadings %>% filter(PC <= PC_to_compare))$loading); y_range[2] = y_range[2] * 1.1 # extra space for avg_expl_var_lab
  var_set_all = c(leading_var, setdiff(unique(res_PCA_loadings$variable), leading_var))
  var_set_split = split_var(var_set_all, var_split_length)
  n_year = uniqueN(res_PCA_loadings$year) - ifelse('Avg' %in% res_PCA_loadings$year, 1, 0)
  for (df_type in unique(res_PCA_loadings$data)){
    for (recov_met in unique(res_PCA_loadings$method)){
      n_row = length(unique(res_PCA_loadings$PCA))  * length(var_set_split)  # number of PCA met * number of variables (split in multiple rows)
      n_col = max(unlist(lapply(var_set_split, length)))      # max number of variables per column
      
      # adjust sign according to leading variable
      res_PCA_loadings_adj = sort_loading_PCA(res_PCA_loadings, leading_var, leading_sign) %>%
        filter(data == df_type) %>%
        filter(method == recov_met)
      row_list = list()
      
      for (var_s in 1:length(var_set_split)){ # variable set
        var_set = var_set_split[[var_s]]
        
        for (j_var in 1:var_split_length){ # used to include empty plot in row_list (last row may have less elements than the others)
          
          var = var_set[j_var]
          i = 1
          for (pc_met in unique(res_PCA_loadings$PCA)){
            
            if (j_var <= length(var_set)){
              if (i == 1){
                tit_lab = paste0("<span style='font-size:21'><p><b>", var, "</b></p><span style='font-size:15'><p>", pc_met, "</p>")
                tit_lab <- rich_text_grob(tit_lab, x = unit(3, "lines"), y = unit(2, "lines"))
              } else {
                tit_lab = paste0("<span style='font-size:15'>", pc_met)
                tit_lab <- rich_text_grob(tit_lab, x = unit(0, "lines"), y = unit(1, "lines"))
              }
              if (i == 1){grad = colorRampPalette(c('#f0f0f0', 'black'))}  # greyscale 
              if (i == 2){grad = colorRampPalette(c('#deebf7', '#3182bd'))}  # bluescale  
              if (i == 3){grad = colorRampPalette(c('#fee6ce', '#e6550d'))}  # redscale  
              if (i == 4){grad = colorRampPalette(c('#e5f5e0', '#31a354'))}  # greenscale 
              avg_expl_var = res_PCA_importance %>%
                filter(PCA == pc_met & method == recov_met & data == df_type & PC %in% paste0('PC', c(1:PC_to_compare))) %>%
                group_by(PC) %>%
                summarise(AVG = mean(`Proportion of Variance` * 100),
                          STD = sd(`Proportion of Variance` * 100), .groups = 'drop') %>%
                arrange(PC)
              avg_expl_var_lab = paste0('PC', c(1:PC_to_compare), '\nAvg Expl Var:\n', round(avg_expl_var$AVG), ' ± ', round(avg_expl_var$STD), '%')
              
              p = ggplot(res_PCA_loadings_adj %>%
                           filter(data == df_type) %>%
                           filter(variable == var) %>%
                           filter(PCA == pc_met) %>%
                           filter(method == recov_met) %>%
                           filter(PC <= PC_to_compare) %>%
                           mutate(year = as.factor(year), PC = paste0('PC ', PC)),
                         aes(fill=year, y=loading, x=PC)) + 
                geom_bar(position="dodge", stat="identity") +
                # facet_wrap(~PC,scales = "free_x") + 
                scale_fill_manual(values = c(rev(grad(n_year)), 'darkgoldenrod1')) +
                ylim(y_range[1], y_range[2]) +
                geom_vline(xintercept = c(1:(PC_to_compare - 1)) + 0.5) +
                geom_hline(aes(yintercept = -signif_thresh_plot), color = 'red', linetype = 'twodash', size = 0.4) +
                geom_hline(aes(yintercept = signif_thresh_plot), color = 'red', linetype = 'twodash', size = 0.4) +
                theme(axis.title.y=element_blank(),
                      axis.title.x=element_blank(),
                      legend.position = "left",
                      legend.text=element_text(size=20),
                      legend.title=element_text(size=20),
                      legend.key.size = unit(0.35, "cm"),
                      plot.margin = ggplot2::margin(0.7, 0, 0.7, 0.3, "cm"),
                      panel.background = element_rect(fill = "white", colour = "black"),
                      panel.grid.major.y = element_line(colour = "black", linetype = 'dashed', size = 0.2),
                      # panel.grid.minor.y = element_line(colour = "black", linetype = 'dashed', size = 0.2),
                      # strip.text.x = element_blank(),
                      axis.text.y=element_text(size=20),
                      axis.text.x=element_blank(),
                      axis.ticks.x = element_blank()) +
                annotate("text", x = c(1:PC_to_compare), y = y_range[2], label = avg_expl_var_lab, hjust = 0.5, vjust = 1, size = 5)
              # show x label only on bottom row and odd bottom-1 row (if any)
              # if ((length(unique(res_PCA_loadings$PCA))*(var_s - 1) + i == n_row) |
              #     (j_var == n_col & length(var_set_all) %% var_split_length != 0 & var_s == (length(var_set_split) - 1) & i == (n_row / length(var_set_split)))){
              #   p = p + theme(axis.text.x = element_text(color = 'black'))
              # }
              # hide legend
              if (var != var_set[1]){
                p = p + theme(axis.text.y=element_blank(),
                              legend.position = "none")
              }
              p = ggplotGrob(p)
              p$grobs[[16]] <- tit_lab
            } else {  # empty plot
              p = ggplotGrob(ggplot() + theme(panel.background = element_blank()))
            }
            p = gtable_add_padding(p, unit(c(0,0,1.5,0), "cm")) # t,r,b,l
            
            row_list[[toString(j_var)]][[length(unique(res_PCA_loadings$PCA))*(var_s - 1) + i]] = p
            i = i + 1
          } # pc_met
        } # j_var
      } #var_s
      
      col_list = list()
      for (j in c(1:n_col)){
        # fill different gtable with empty cols/rows
        max_row = max(unlist(lapply(row_list[[j]], nrow)))
        max_col = max(unlist(lapply(row_list[[j]], ncol)))
        for (i_c in 1:length(row_list[[j]])){
          if(nrow(row_list[[j]][[i_c]]) < max_row){
            while (max_row - nrow(row_list[[j]][[i_c]]) > 0){
              row_list[[j]][[i_c]] = gtable_add_rows(row_list[[j]][[i_c]], unit(1, "null"))
            }
          }
          if(ncol(row_list[[j]][[i_c]]) < max_col){
            while (max_col - ncol(row_list[[j]][[i_c]]) > 0){
              row_list[[j]][[i_c]] = gtable_add_cols(row_list[[j]][[i_c]], unit(1, "null"))
            }
          }
        }
        # create list of columns
        col_list[[j]] = do.call(rbind, c(row_list[[j]], size="last"))
      }
      g = do.call(cbind, c(col_list, size="last"))
      g = gtable_add_padding(g, unit(2, "cm")) # t,r,b,l
      png(paste0('./Stats/2_PCA_Loading_plot_by_variable_', df_type, '_', recov_met, '.png'), width = 32, height = 13, units = 'in', res=300)
      grid.draw(g)
      dev.off()
      
    } # recov_met
  } # df_type
}
# all variables together
signif_thresh_plot = 0.2 # used to highlight significative loadings
leading_var = "growth"
leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
PC_to_compare = 2
{
  add_recov = c()#ifelse(length(unique(res_PCA_loadings$method)) >= 2, 'DIFFERENCE', c())
  y_range = range((res_PCA_loadings %>% filter(PC <= PC_to_compare))$loading); y_range[2] = y_range[2] * 1.1 # extra space for avg_expl_var_lab
  n_year = uniqueN(res_PCA_loadings$year) - ifelse('Avg' %in% res_PCA_loadings$year, 1, 0)
  for (df_type in unique(res_PCA_loadings$data)){
    for (recov_met in c(unique(res_PCA_loadings$method), add_recov)){
      
      # adjust sign according to leading variable
      if (recov_met != 'DIFFERENCE'){
        res_PCA_loadings_adj = sort_loading_PCA(res_PCA_loadings, leading_var, leading_sign) %>%
          filter(data == df_type) %>%
          filter(method == recov_met) %>%
          mutate(SIGNIF_LAB = '',
                 SIGN = 1)
        # create main title
        avg_expl_var = res_PCA_importance %>%
          filter(method == recov_met & data == df_type & PC %in% paste0('PC', c(1:PC_to_compare))) %>%
          group_by(PCA, PC) %>%
          summarise(AVG = mean(`Proportion of Variance` * 100),
                    STD = sd(`Proportion of Variance` * 100), .groups = 'drop') %>%
          mutate(LAB = paste0(PCA, ' ', round(AVG), ' ± ', round(STD), '%')) %>%
          group_by(PC) %>%
          summarise(LAB = paste0(LAB, collapse = '  |  '), .groups = 'drop') %>%
          mutate(LAB = paste0('Average Explained Variance: ', LAB))
        
        # evaluate difference 
      } else {
        diff_met = unique(res_PCA_loadings$method)[1:2]
        cat('\n ----', df_type, ': Difference evaluated between', diff_met[1], 'and', diff_met[2],'\n')
        res_PCA_loadings_adj = sort_loading_PCA(res_PCA_loadings, leading_var, leading_sign) %>%
          filter(data == df_type) %>%
          filter(method %in% diff_met) %>%
          select(PCA, variable, PC, year, data, method, loading) %>%
          spread(method, loading) %>%
          mutate(ABS_DIFF = abs(abs((!!as.name(diff_met[1]))) - abs((!!as.name(diff_met[2])))),
                 SIGN = sign((!!as.name(diff_met[1]))) * sign((!!as.name(diff_met[2])))) %>%
          mutate(SIGN = ifelse(SIGN == 0, 1, SIGN)) %>%
          mutate(loading = ABS_DIFF * SIGN,
                 method = recov_met) %>%
          mutate((!!as.name(paste0('SIGNIF_', diff_met[1]))) := ifelse(abs((!!as.name(diff_met[1]))) >= signif_thresh_plot, '*', ''),
                 (!!as.name(paste0('SIGNIF_', diff_met[2]))) := ifelse(abs((!!as.name(diff_met[2]))) >= signif_thresh_plot, '°', '')) %>%
          mutate(SIGNIF_LAB = paste0((!!as.name(paste0('SIGNIF_', diff_met[1]))), (!!as.name(paste0('SIGNIF_', diff_met[2])))))
        avg_expl_var = data.frame(PC = paste0('PC', 1:PC_to_compare),
                                  LAB = paste0('Magnitude is Abs[ Abs(', diff_met[1],') - (Abs(', diff_met[2],
                                               ') ] and sign is + if loadings have same sign, - otherwise  |  * ', diff_met[1], '  ° ', diff_met[2],
                                               ' loading above ', signif_thresh_plot))
      }
      
      row_list = list()
      for (pc in 1:PC_to_compare){
        
        if (pc == 1){grad = colorRampPalette(c('#f0f0f0', 'black'))}  # greyscale
        if (pc == 2){grad = colorRampPalette(c('#deebf7', '#3182bd'))}  # bluescale
        if (pc == 3){grad = colorRampPalette(c('#fee6ce', '#e6550d'))}  # redscale
        if (pc == 4){grad = colorRampPalette(c('#e5f5e0', '#31a354'))}  # greenscale
        tit_lab = paste0("<span style='font-size:25'><b>", paste0('PC ', pc, '  '), "</b><span style='font-size:20'>",
                         (avg_expl_var %>% filter(PC == paste0('PC', pc)))$LAB)
        tit_lab <- rich_text_grob(tit_lab,
                                  x = unit(28, "lines"),#-0.62 * uniqueN(res_PCA_loadings_adj$variable) * uniqueN(res_PCA_loadings_adj$year), "lines"), # for 17*8 should be 95
                                  y = unit(2, "lines")) # for 2 PC should be 33    83 -25 * PC_to_compare
        
        p = ggplot(res_PCA_loadings_adj %>%
                     filter(data == df_type) %>%
                     filter(method == recov_met) %>%
                     filter(PC == pc) %>%
                     mutate(variable = split_string(variable, 15)) %>%
                     mutate(year = as.factor(year), PC = paste0('PC ', PC)),
                   aes(fill=year, y=loading, x=variable)) + 
          geom_bar(position="dodge", stat="identity") +
          facet_wrap(~ PCA, ncol = 1, dir = 'v', scales = 'free_y', strip.position = 'left') +
          scale_fill_manual(values = c(rev(grad(n_year)), 'darkgoldenrod1')) +
          ylim(y_range[1], y_range[2]) +
          geom_vline(xintercept = c(1:(uniqueN(res_PCA_loadings_adj$variable) - 1)) + 0.5) +
          # geom_hline(aes(yintercept = -signif_thresh_plot), color = 'red', linetype = 'twodash', size = 0.4) +
          # geom_hline(aes(yintercept = signif_thresh_plot), color = 'red', linetype = 'twodash', size = 0.4) +
          theme(axis.title.y=element_blank(),
                axis.title.x=element_blank(),
                axis.text.x = element_text(angle = 90,  size = 17, vjust=0.3),
                axis.text.y=element_text(size=15),
                legend.text=element_text(size=20),
                legend.title=element_text(size=20),
                strip.text.y = element_text(size = 18,face="bold"),
                plot.margin = ggplot2::margin(2, 0, 0.7, 0.3, "cm"),
                panel.background = element_rect(fill = "white", colour = "black"),
                panel.grid.major.y = element_line(colour = "black", linetype = 'dashed', size = 0.2))+
          scale_x_discrete(position = "bottom") +
          geom_text(aes(label=SIGNIF_LAB, vjust=0.5, hjust = ifelse(SIGN >= 0, -0.2, 1.2)), color="black", size=3.5, position=position_dodge(.9), angle = 90)
        
        p = ggplotGrob(p)
        p$grobs[[16]] <- tit_lab
        p = gtable_add_padding(p, unit(c(0,0,0.35,0), "cm")) # t,r,b,l
        row_list[[pc]] = p
      } # pc
      g = do.call(rbind, c(row_list, size="last"))
      g = gtable_add_padding(g, unit(2, "cm")) # t,r,b,l
      png(paste0('./Stats/2_PCA_Loading_plot_all_', df_type, '_', recov_met, '.png'), width = 22, height = 17, units = 'in', res=300)
      grid.draw(g)
      dev.off()
      
    } # recov_met
  } # df_type
}




recov_met_set = c('SOFT_IMPUTE')

### evaluate PCA final index based on PC (=scores)
PC_to_keep = 2
index_1_thresh = 0  # threshold to split index 1 (first PC scores)
index_2_thresh = 0  # threshold to split index 2 (second PC scores) - if PC_to_keep == 2
load_thresh = 0  # loadings with magnitude below threshold are set to 0 and scores are evaluated always as X * loadings
leading_var = "growth"
leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
{
  res_index_PCA = list()
  for (df_type in unique(res_PCA_loadings$data)){
    for (recov_met in recov_met_set){   #unique(res_PCA_loadings$method)){
      for (pc_met in unique(res_PCA_loadings$PCA)){
        
        cat('\n evaluating data:', df_type, '- recov_met:', recov_met, '- pc_met:', pc_met)
        res_index_PCA = evaluate_index_PCA(res_index_PCA, res_PCA_list, res_PCA_loadings, res_PCA_importance, PC_to_keep, index_1_thresh, index_2_thresh,
                                           load_thresh, leading_var, leading_sign, df_type, recov_met, pc_met)
        
      } # pc_met
    } # recov_met
  } # df_type
  
  saveRDS(res_index_PCA, './Checkpoints/res_index_PCA.rds')
}
res_index_PCA = readRDS('./Checkpoints/res_index_PCA.rds')
















df_type_to_keep = c("RestrictedDifference", "Restricted")
recov_met_to_keep = c("SOFT_IMPUTE")


### Dynamic Factor Model
n_factor_tot = 1  # number of factors to evaluate (all incremental sets up to n_factor_tot are evaluated)
VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
univ_reload = T  # reload previous evaluation (if available) for evaluate_DFM_univar
univ_save = F  # save intermediate evaluation
total_save = T  # save RDS for univariate + multivariate for each pair df_type + recov_met
{
  res_DFM_factors = res_DFM_loadings = res_DFM_stats = res_DFM_MAPE = c()
  sink(paste0('./Log/DFM_fitting_log_', format(Sys.time(), "%Y-%m-%dT%H%M"), '.txt'), append = F, split = T, type = 'output')
  cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
  for (df_type in df_type_to_keep){ #unique(df_final$data)){
    for (recov_met in recov_met_to_keep){   # unique(df_final$method)){
      
      cat('\n\n\n-----------------------################# evaluating data:', df_type, '- recov_met:', recov_met, '#################-----------------------\n')
      
      RDS_lab = paste0(gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '_', df_type, '_', recov_met, '.rds')
      res_DFM_list = list()
      
      # selecting the data set
      df_SP_orig = df_final %>%
        filter(year != 'Avg') %>%
        filter(method == recov_met) %>%
        filter(data == df_type) %>%
        select(-method, -data) %>%
        mutate(year = as.numeric(year)) %>%
        spread(variable, val) %>% # spread uses alphabetical ordering for columns
        select(-starts_with('GEO')) %>%
        arrange(country, year)
      
      # ID must be integer
      country_code = data.frame(country = unique(df_SP_orig$country), ID = 1:uniqueN(df_SP_orig$country), stringsAsFactors = F)
      
      df_SP = df_SP_orig %>%
        left_join(country_code, by = "country") %>%
        rename(TIME = year) %>%
        select(-country) %>%
        select(ID, everything())
      
      # variable names must be numeric
      variable_code = data.frame(variable = colnames(df_SP)[-c(1,2)], CODE = as.character(1:(ncol(df_SP) - 2)), stringsAsFactors = F)
      colnames(df_SP)[-c(1,2)] = 1:nrow(variable_code)
      
      # reload previous evaluation is univ_reload = T
      reload_err = try(res_DFM_list_reload <- suppressWarnings(readRDS(paste0('./Checkpoints/DFM/', RDS_lab))), silent = T)
      if (class(reload_err) == "try-error"){
        res_DFM_list_reload = NULL
      }
      
      # evaluate model for each set of factors (1, 2, 3, ..., n_factors)
      for (n_factor in 1:n_factor_tot){
        
        cat('\n\n\n  --------------------------------- Testing', n_factor, ifelse(n_factor == 1, 'factor', 'factors') ,'model   ---------------------------------\n\n')
        
        for (country_i in 1:nrow(country_code)){
          
          # Dynamic Factor Model - evaluate for single country (data are standardized)
          DFM_uni = evaluate_DFM_univar(res_DFM_factors, res_DFM_loadings, res_DFM_stats, res_DFM_list, res_DFM_list_reload,
                                        variable_code, country_code, df_SP, n_factor, country_i, df_type, recov_met,
                                        univ_reload = univ_reload,
                                        dfm_eval = 'BFGS',  # 'BFGS' or 'kem' - kem seems not to converge due to maxiter
                                        dfm_max_iter = 3000)
          res_DFM_stats = DFM_uni$res_DFM_stats
          res_DFM_list = DFM_uni$res_DFM_list
          res_DFM_factors = DFM_uni$res_DFM_factors
          res_DFM_loadings = DFM_uni$res_DFM_loadings
          
          if (univ_save){
            saveRDS(res_DFM_list, paste0('./Checkpoints/DFM/', RDS_lab))
          }
          
        } # country_i
        
        # adjust DFM taking into account all countries together (factors are standardized)
        cat('\n\n\n  --------------------------------- Filtering all countries together with', n_factor, ifelse(n_factor == 1, 'factor', 'factors') ,'models   ---------------------------------\n\n')
        
        DFM_multi = evaluate_DFM_multivar(res_DFM_factors, res_DFM_loadings, res_DFM_stats, res_DFM_list, res_DFM_MAPE,
                                          df_SP_orig, df_type, recov_met, n_factor,
                                          VAR_alpha = VAR_alpha,
                                          kalm_Q_hat_mode = kalm_Q_hat_mode,
                                          multiv_reload_Kalm = T, res_DFM_list_reload = res_DFM_list_reload)
        res_DFM_stats = DFM_multi$res_DFM_stats
        res_DFM_list = DFM_multi$res_DFM_list
        res_DFM_factors = DFM_multi$res_DFM_factors
        res_DFM_loadings = DFM_multi$res_DFM_loadings
        res_DFM_MAPE = DFM_multi$res_DFM_MAPE
      } # n_factor
      
      if (total_save){
        cat('\n\n** saving... ')
        saveRDS(res_DFM_list, paste0('./Checkpoints/DFM/', RDS_lab))
        cat('Done')
      }
      
    } # recov_met
  } # df_type
  cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
  sink()
  
  # evaluate final metrics to choose the optimal number of factors (AIC or Explain Variance, i.e. R^2) -> Explain variance on 95th percentile according to RSS
  # univariate DFM is forced to have the same number of factors of corresponding multivariate DFM
  
  res_DFM_best_model = res_DFM_stats %>%
    group_by(data, method, DFM, Total_Factors) %>%
    summarise(Total_AICc = sum(AICc, na.rm = T),
              Explain_var = 1 - sum(y_RSS, na.rm = T)/ sum(y_TSS, na.rm = T),
              Explain_var_95 = 1 - sum(y_RSS_95, na.rm = T)/ sum(y_TSS_95, na.rm = T),
              Explain_var_99 = 1 - sum(y_RSS_99, na.rm = T)/ sum(y_TSS_99, na.rm = T),
              Total_LogLik = sum(LogLik, na.rm = T),
              Total_y_RSS = sum(y_RSS, na.rm = T),
              Total_y_RSS_95 = sum(y_RSS_95, na.rm = T),
              Total_y_RSS_99 = sum(y_RSS_99, na.rm = T),
              Total_y_TSS = sum(y_TSS, na.rm = T),
              Total_y_TSS_95 = sum(y_TSS_95, na.rm = T),
              Total_y_TSS_99 = sum(y_TSS_99, na.rm = T),
              Null_values = sum(is.na(AICc)), .groups = 'drop') %>%
    ungroup() %>%
    group_by(data, method, DFM) %>%
    arrange(desc(Explain_var_95)) %>%
    mutate(Best_model = ifelse(row_number()==1, 'YES', '')) %>% # force DFM_univar factors
    group_by(data, method) %>%
    mutate(best_DFM_factor = Total_Factors[Best_model == 'YES' & DFM == 'DFM_multivar']) %>%
    mutate(Best_model = ifelse(Total_Factors == best_DFM_factor, 'YES', '')) %>%
    select(-best_DFM_factor) %>%
    as.data.frame()
  
  cat('\n\n\n ----- Best number of factors:\n')
  print(res_DFM_best_model %>% filter(Best_model == 'YES') %>% select(data, method, DFM, Total_Factors, Explain_var_95) %>% arrange(DFM, data, method))
  
  # save results
  write.table(res_DFM_stats, paste0('./Stats/2_DFM_Stats_fitting_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".", na = "")
  write.table(res_DFM_best_model, paste0('./Stats/2_DFM_Stats_factors_selection_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  write.table(res_DFM_MAPE, paste0('./Stats/2_DFM_Stats_MAPE_list_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  saveRDS(res_DFM_factors, paste0('./Checkpoints/res_DFM_factors_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  saveRDS(res_DFM_loadings, paste0('./Checkpoints/res_DFM_loadings_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  saveRDS(res_DFM_stats, paste0('./Checkpoints/res_DFM_stats_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  saveRDS(res_DFM_best_model, paste0('./Checkpoints/res_DFM_best_model_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  saveRDS(res_DFM_MAPE, paste0('./Checkpoints/res_DFM_MAPE_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
}
res_DFM_factors = readRDS(paste0('./Checkpoints/res_DFM_factors_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
res_DFM_loadings = readRDS(paste0('./Checkpoints/res_DFM_loadings_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
res_DFM_stats = readRDS(paste0('./Checkpoints/res_DFM_stats_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
res_DFM_best_model = readRDS(paste0('./Checkpoints/res_DFM_best_model_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
res_DFM_MAPE = readRDS(paste0('./Checkpoints/res_DFM_MAPE_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))





dfm_met_to_plot = c("DFM_multivar")


### evaluate DFM final index based on Factors
# only 1 and 2 factors plot are showed
index_1_thresh = 0  # threshold to split index 1 (first factor)
index_2_thresh = 0  # threshold to split index 2 (second factor)
load_thresh = 0  # loadings with magnitude below threshold are set to 0
leading_var = "growth"
leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
{
  res_index_DFM = list()
  for (df_type in setdiff(unique(res_DFM_loadings$data), c('ExtendedDifference', 'ExtendedOriginal'))){
    for (recov_met in unique(res_DFM_loadings$method)){
      
      RDS_lab = paste0(gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '_', df_type, '_', recov_met, '.rds')
      res_DFM_list = readRDS(paste0('./Checkpoints/DFM/', RDS_lab))
      
      for (dfm_met in dfm_met_to_plot){
        
        res_index_DFM = evaluate_index_DFM(res_index_DFM, res_DFM_best_model, res_DFM_factors, res_DFM_list, res_DFM_loadings, index_1_thresh, index_2_thresh,
                                           load_thresh, leading_var, leading_sign, df_type, recov_met, dfm_met, expl_var_to_show=95)
        
      } # dfm_met
    } # recov_met
  } # df_type
  
  saveRDS(res_index_DFM, './Checkpoints/res_index_DFM.rds')
}
res_index_DFM = readRDS('./Checkpoints/res_index_DFM.rds')



### plot heatmap for matrix A of factors interaction between countries
VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
{
  for (df_type in setdiff(unique(res_DFM_loadings$data), c('ExtendedDifference', 'ExtendedOriginal'))){
    for (recov_met in unique(res_DFM_loadings$method)){
      
      RDS_lab = paste0(gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '_', df_type, '_', recov_met, '.rds')
      res_DFM_list = readRDS(paste0('./Checkpoints/DFM/', RDS_lab))
      
      for (dfm_met in dfm_met_to_plot){
        
        best_model = res_DFM_best_model %>%
          filter(data == df_type) %>%
          filter(method == recov_met) %>%
          filter(DFM == dfm_met) %>%
          filter(Best_model == 'YES')
        n_factor_best = best_model$Total_Factors
        expl_var = best_model$Explain_var
        
        if (dfm_met == 'DFM_univar'){
          mat_names = names(res_DFM_list[[df_type]][[recov_met]][[dfm_met]][[paste0(n_factor_best, '_factors')]])
          plot_mat = list()
          pm = 1
          mat_col = c()
          for (nam in mat_names){
            plot_mat[[pm]] = res_DFM_list[[df_type]][[recov_met]][[dfm_met]][[paste0(n_factor_best, '_factors')]][[nam]][['A']]
            pm = pm + 1
            mat_col = c(mat_col, paste0(nam,paste0('_', 1:2)))
          }
          plot_mat = bdiag(plot_mat) %>%
            as.matrix() %>%
            `colnames<-`(mat_col) %>%
            `rownames<-`(mat_col) %>%
            as.table() %>%
            as.data.frame() %>%
            rename(Factor = Freq)
        } else {
          plot_mat = res_DFM_list[[df_type]][[recov_met]][[dfm_met]][[paste0(n_factor_best, '_factors')]][['All']][['A_hat']] %>%
            as.table() %>%
            as.data.frame() %>%
            rename(Factor = Freq)
        }
        
        png(paste0('./Results/1b_', df_type, '_', recov_met, '_', dfm_met, '_factors_interaction.png'), width = 22, height = 22, units = 'in', res=300) 
        plot(ggplot(plot_mat, aes(x = Var1, y = Var2, fill = Factor)) +
               geom_tile(colour="grey",size=0.25) +
               scale_fill_gradient2() +
               theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
               labs(x = '', y = '') +
               ggtitle(paste0('Interaction between factors - Explained Variance: ', round(expl_var * 100), '%')))
        dev.off()
        
        write.table(plot_mat %>%
                      filter(Factor != 0) %>%
                      arrange(desc(abs(Factor))),
                    paste0('./Results/1b_', df_type, '_', recov_met, '_', dfm_met, '_factors_interaction.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
        
      } # dfm_met
    } # recov_met
  } # df_type
  
}





df_type_to_keep = c("RestrictedDifference", "Restricted")
recov_met_to_keep = c("SOFT_IMPUTE")
PCA_to_keep = c('RobPCA')
DFM_to_keep = c('DFM_multivar')


### plot comparison of index evolution over time for all methods
res_PCA_loadings = readRDS('./Checkpoints/res_PCA_loadings.rds')
res_index_PCA = readRDS('./Checkpoints/res_index_PCA.rds')
VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
res_DFM_factors = readRDS(paste0('./Checkpoints/res_DFM_factors_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
res_DFM_loadings = readRDS(paste0('./Checkpoints/res_DFM_loadings_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
res_DFM_best_model = readRDS(paste0('./Checkpoints/res_DFM_best_model_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
res_index_DFM = readRDS('./Checkpoints/res_index_DFM.rds')
algo_to_plot = c('RobPCA', 'DFM_multivar')
{
  # create single list of evaluated index (raw scores)
  res_ALL_scores = c()
  # from PCA
  for (df_type in df_type_to_keep){
    for (recov_met in recov_met_to_keep){
      for (pc_met in PCA_to_keep){
        
        pc_mat = res_index_PCA[[df_type]][[recov_met]][[pc_met]][['scores_raw']]
        new_col = paste0(colnames(pc_mat), pc_mat[1,])
        res_ALL_scores = res_ALL_scores %>% bind_rows(
          pc_mat %>%
            setNames(new_col) %>%
            filter(row_number() != 1) %>%
            `rownames<-`(pc_mat$COUNTRY[-1]) %>%
            select(-'COUNTRYINDEX') %>%
            select(-starts_with('Avg')) %>%
            as.matrix() %>%
            as.table() %>%
            as.data.frame(stringsAsFactors = F) %>%
            mutate(year = as.numeric(substr(Var2, 1, 4)),
                   factor = as.numeric(substr(Var2, 7, 7)),
                   index = as.numeric(as.character(Freq)),
                   method = recov_met,
                   data = df_type,
                   algo = pc_met,
                   family = 'PCA') %>%
            rename(country = Var1) %>%
            select(-Var2, -Freq) %>%
            left_join(res_index_PCA[[df_type]][[recov_met]][[pc_met]][['Expl_Var']] %>%
                        filter(year != 'Avg') %>%
                        mutate(year = as.numeric(year),
                               Expl_Var = as.numeric(Expl_Var)) %>%
                        rename(Explain_var = Expl_Var) %>%
                        select(-PC), by = "year") %>%
            mutate(Explain_var_95 = Explain_var,
                   Explain_var_99 = Explain_var)
        )
        
      } # pc_met
    } # recov_met
  } # df_type
  
  # from DFM
  for (df_type in df_type_to_keep){
    for (recov_met in recov_met_to_keep){
      for (dfm_met in DFM_to_keep){
        
        best_model = res_DFM_best_model %>%
          filter(data == df_type) %>%
          filter(method == recov_met) %>%
          filter(DFM == dfm_met) %>%
          filter(Best_model == 'YES')
        n_factor_best = best_model$Total_Factors
        
        res_ALL_scores = res_ALL_scores %>% bind_rows(
          res_DFM_factors %>%
            filter(data == df_type) %>%
            filter(method == recov_met) %>%
            filter(DFM == dfm_met) %>%
            filter(Total_Factors == n_factor_best) %>%
            rename(algo = DFM,
                   factor = Factor,
                   index = val) %>%
            mutate(family = 'DFM') %>%
            select(-Total_Factors, -Var_removed) %>%
            cbind(best_model %>% select(starts_with('Explain')))
        )
      } # dfm_met
    } # recov_met
  } # df_type
  
  # plot for each country
  for (df_type in unique(res_ALL_scores$data)){
    for (recov_met in unique(res_ALL_scores$method)){
      
      comp_data = res_ALL_scores %>%
        filter(data == df_type) %>%
        filter(method == recov_met)  %>%
        filter(algo %in% algo_to_plot) %>%
        mutate(factor = paste0('Index ', factor),
               color_id = paste0(family, '_', algo))
      # filter(country %in% country_list[((con - 1) * n_col + 1):min(c(con * n_col, length(country_list)))])
      
      # create color palette for each family
      grad_bl = c('blue3', 'deepskyblue', 'deepskyblue3')  # bluescale
      grad_rd = c('brown', 'brown1', 'darkorange')  # redscale
      grad_gr = c('chartreuse4', 'chartreuse3', 'chartreuse')  # greenscale
      color_palette_list = comp_data %>%
        select(family, color_id) %>%
        unique() %>%
        arrange(color_id) %>%
        group_by(family) %>%
        summarise(rep = n(), .groups = 'drop')
      color_palette = c()
      for (pal in 1:nrow(color_palette_list)){
        if (pal == 1){color_palette = c(color_palette, grad_bl[1:color_palette_list$rep[pal]])}
        if (pal == 2){color_palette = c(color_palette, grad_rd[1:color_palette_list$rep[pal]])}
        if (pal == 3){color_palette = c(color_palette, grad_gr[1:color_palette_list$rep[pal]])}
      }
      
      png(paste0('./Results/2_', df_type, '_', recov_met, '_factors_evolution_over_time.png'), width = 22, height = 50, units = 'in', res=300) 
      plot(ggplot(comp_data,
                  aes(fill=algo, y=index, x=year)) + 
             geom_line(aes(colour = color_id, linetype = family), size = 1) +
             scale_linetype_manual(values=c('solid', 'twodash', 'dotted', 'dotdash'), name = 'Algo type') +
             scale_colour_manual(values=color_palette, name = "Algo", labels = (comp_data %>% select(color_id, algo) %>% unique() %>% arrange(color_id))$algo) +
             # facet_grid(country ~ factor, scales = "free_y", switch = 'y') +
             facet_wrap(country ~ factor, dir = 'v', scales = 'free', strip.position = 'top', ncol = 7) +
             scale_x_continuous(labels = as.character(sort(unique(comp_data$year))), breaks = unique(comp_data$year)) +
             theme(panel.background = element_rect(fill = "white", colour = "black"),
                   panel.grid = element_line(colour = "black", linetype = 'dashed', size = 0.2),
                   strip.text = element_text(size = 8)) +
             labs(x = '', y = ''))
      dev.off()
      
    } # recov_met
  } # df_type
  
  saveRDS(res_ALL_scores, './Checkpoints/res_ALL_scores.rds')
  
}
res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')



### plot variation coefficient over time
res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
algo_to_plot = c('RobPCA', 'DFM_multivar')
algo_to_plot_lab = c('RobPCA', 'DFM')
{
  for (df_type in unique(res_ALL_scores$data)){
    for (recov_met in unique(res_ALL_scores$method)){
      
      data = res_ALL_scores %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        filter(algo %in% algo_to_plot) %>%
        mutate(factor = paste('Index', factor)) %>%
        group_by(algo, factor, country) %>%
        # filter(factor == max(factor)) %>%
        mutate(index = index + min(index)) %>%  # in order to make all values positive andd avoid 0 mean
        summarise(Variation_coeff = sd(index) / abs(mean(index)), .groups = 'drop') %>%
        mutate(quantile = quantile(Variation_coeff, 0.90)) %>%
        ungroup() %>%
        left_join(data.frame(old = algo_to_plot, new = algo_to_plot_lab, stringsAsFactors = F), by = c('algo' = 'old')) %>%
        select(-algo) %>%
        rename(algo = new) %>%
        filter(Variation_coeff <= quantile)
      
      png(paste0('./Results/2_', df_type, '_', recov_met, '_factors_variation_coeff.png'), width = 10, height = 10, units = 'in', res=300) 
      plot(ggplot(data,
                  aes(x=algo, y=Variation_coeff, color=algo)) +
             geom_boxplot(lwd=1.5) +
             scale_x_discrete(limits=sort(algo_to_plot_lab)) +
             scale_color_manual(values=c("blue3", "brown", "chartreuse4")) +
             labs(y = "Coefficient of Variation") +
             facet_wrap(~ factor, dir = 'v', scales = 'free_x', strip.position = 'top', ncol = 2) +
             theme(legend.position = "none",
                   axis.title.x = element_blank(),
                   text = element_text(size=22))
      )
      dev.off()
    } # recov_met
  } # df_type
}



### plot distributions of index
{
  res_index_distribution = c()
  for (df_type in unique(res_ALL_scores$data)){
    for (recov_met in unique(res_ALL_scores$method)){
      
      symmetric_log = function(x, base = 10){
        if (x == 0){out = 0
        } else if (x >= -1 & x <= 1){out = x
        } else if (x > 1){out = log(x, base = base)
        } else if (x < -1){out = -log(-x, base = base)}
        
        return(out)
      }
      symmetric_log = Vectorize(symmetric_log)
      
      data = res_ALL_scores %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        mutate(index_log = symmetric_log(index))
      
      p = ggplot(data = data, aes(x = index_log)) +
        geom_density(aes(fill = as.factor(paste0('Index ',factor))), alpha = 0.5) +
        facet_wrap(algo ~ year, ncol = uniqueN(data$year), dir = 'h', scales = 'free_y') +
        ggtitle('Index scale is linear in [-1, 1] and logarithmic otherwise\n') +
        labs(fill = 'Index') +
        theme(
          plot.title = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size=12, face = "bold"))
      
      res_index_distribution = res_index_distribution %>% bind_rows(
        data.frame(data = df_type,
                   method = recov_met,
                   min = min(data$index),
                   perc05 = quantile(data$index, 0.05),
                   perc95 = quantile(data$index, 0.95),
                   max = max(data$index), stringsAsFactors = F)
      )
      
      png(paste0('./Results/2_', df_type, '_', recov_met, '_index_distributions.png'), width = 22, height = 17, units = 'in', res=300)
      plot(p)
      dev.off()
      
    } # recov_met
  } # df_type
  
  write.table(res_index_distribution, './Results/2_index_distributions.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  print(res_index_distribution)
}



### add additional variables
NA_toll = 10  # max tolerance (%) of NAs
{
  ref_country_names = read.csv("./Data/Data_set.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)$country %>% unique()
  
  # load dataset
  {
    # WEF_Global competitiveness data_2006-2017.csv
    add_WEF = read.csv("./Data/Additional_Variables/WEF_Global competitiveness data_2006-2017.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F, skip = 1) %>%
      mutate_all(funs(replace(., . == "", NA))) %>%
      gather('country', 'val', -c(Year, Series)) %>%
      spread(Series, val) %>%
      setnames('Year', 'year') %>%
      mutate(country = gsub('\\.\\.', ', ', country)) %>%
      mutate(country = gsub('\\.', ' ', country)) %>%
      select('country', everything()) %>%
      mutate_at(vars(-country), as.numeric)
    
    set_bef = setdiff(ref_country_names, unique(add_WEF$country))
    
    rep_tab = read.csv("./Data/Additional_Variables/WEF_Global competitiveness data_2006-2017 - match.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)
    for (i in 1:nrow(rep_tab)){
      add_WEF = add_WEF %>% mutate(country = gsub(rep_tab$ORIGINAL[i], rep_tab$REPLACE[i], country))
    }
    
    set_aft = setdiff(ref_country_names, unique(add_WEF$country))
    if (!length(set_aft) + nrow(rep_tab) == length(set_bef)){'\n ***** Error on WEF countries match'}
    
    # WB_WDI_2005-2017.csv
    add_WB_WDI = read.csv("./Data/Additional_Variables/WB_WDI_2005-2017.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
      select(-Country.Code) %>%
      gather('year', 'val', -c(Country.Name, Indicator.Name)) %>%
      mutate(year = as.numeric(gsub('X', '', year))) %>%
      spread(Indicator.Name, val) %>%
      setnames('Country.Name', 'country')
    
    set_bef = setdiff(ref_country_names, unique(add_WB_WDI$country))
    
    rep_tab = read.csv("./Data/Additional_Variables/WB_WDI_2005-2017 - match.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)
    for (i in 1:nrow(rep_tab)){
      add_WB_WDI = add_WB_WDI %>% mutate(country = gsub(rep_tab$ORIGINAL[i], rep_tab$REPLACE[i], country))
    }
    
    set_aft = setdiff(ref_country_names, unique(add_WB_WDI$country))
    if (!length(set_aft) + nrow(rep_tab) == length(set_bef)){'\n ***** Error on WB_WDI countries match'}
    
    # WB_Financial structure and development data_2005-2017.csv
    add_WB_Fin = read.csv("./Data/Additional_Variables/WB_Financial structure and development data_2005-2017.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
      mutate_all(funs(replace(., . == "", NA))) %>%
      select(-WB.COUNTRY.CODE, -WB.REGION, -WB.INCOME.GROUP) %>%
      setNames(tolower(names(.))) %>%
      mutate_at(vars(-country), as.numeric)
    
    set_bef = setdiff(ref_country_names, unique(add_WB_Fin$country))
    
    rep_tab = read.csv("./Data/Additional_Variables/WB_Financial structure and development data_2005-2017 - match.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)
    for (i in 1:nrow(rep_tab)){
      add_WB_Fin = add_WB_Fin %>% mutate(country = gsub(rep_tab$ORIGINAL[i], rep_tab$REPLACE[i], country))
    }
    
    set_aft = setdiff(ref_country_names, unique(add_WB_Fin$country))
    if (!length(set_aft) + nrow(rep_tab) == length(set_bef)){'\n ***** Error on WB_Fin countries match'}
    
    # IMF_fiscal position_2006-2018.csv
    add_IMF_Fis = read.csv("./Data/Additional_Variables/IMF_fiscal position_2006-2018.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
      select(-Country.Code) %>%
      setnames('Country.Name', 'country') %>%
      setNames(tolower(names(.)))
    
    # IMF Financial market development data_2005-2016.csv
    add_IMF_FMa = read.csv("./Data/Additional_Variables/IMF Financial market development data_2005-2016.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
      select(-Country.Code) %>%
      mutate_all(funs(replace(., . == 0, NA))) %>%
      setnames(c('Country.Name', 'Time.Period'), c('country', 'year')) %>%
      setNames(tolower(names(.)))
    
    # IMF Financial access data_2004-2016.csv
    add_IMF_FAc = read.csv("./Data/Additional_Variables/IMF Financial access data_2004-2016.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
      select(-c_code1)
  }
  
  # merge data
  df_add = df %>% select(country, year)
  df_add = add_variables(df_add, df_test = add_WEF, df_test_lab = 'add_WEF', NA_toll =  NA_toll, summary_fold = './Stats/')
  df_add = add_variables(df_add, df_test = add_WB_WDI, df_test_lab = 'add_WB_WDI', NA_toll =  NA_toll, summary_fold = './Stats/')
  df_add = add_variables(df_add, df_test = add_WB_Fin, df_test_lab = 'add_WB_Fin', NA_toll =  NA_toll, summary_fold = './Stats/')
  df_add = add_variables(df_add, df_test = add_IMF_Fis, df_test_lab = 'add_IMF_Fis', NA_toll =  NA_toll, summary_fold = './Stats/')
  df_add = add_variables(df_add, df_test = add_IMF_FMa, df_test_lab = 'add_IMF_FMa', NA_toll =  NA_toll, summary_fold = './Stats/')
  df_add = add_variables(df_add, df_test = add_IMF_FAc, df_test_lab = 'add_IMF_FAc', NA_toll =  NA_toll, summary_fold = './Stats/')
  
  # stats
  add_summary = variable_stats(df_add)
  
  write.table(add_summary, "./Stats/3_Additional_variable_summary.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  # save
  saveRDS(df_add, './Checkpoints/df_add.rds')
}
df_add = readRDS('./Checkpoints/df_add.rds')



### evaluate sensitivity on index threshold - execute regression task on test_variable_set with both original regressors and binarized indices
df_set = c('Original', 'Difference')
recov_set = c('TENSOR_BF')
fitting_set = c('RobPCA', 'DFM_multivar')
index_1_set = c(-1.5, -1 , -0.5, 0, 0.5, 1, 1.5)
index_2_set = c(-1.5, -1 , -0.5, 0, 0.5, 1, 1.5)
test_variable_set = c('var_add_WB_WDI_Bank_nonperforming_loans_to_total_gross_loans_PERC',
                      'var_add_WB_WDI_GDP_per_capita_current_USDOLL',
                      'var_add_WB_WDI_GDP_per_capita_growth_annual_PERC',
                      'var_add_WB_WDI_GDP_per_capita_PPP_current_international_DOLL',
                      'var_add_WB_WDI_Domestic_credit_to_private_sector_PERC_of_GDP')
flag_tuning = T  # tune models parameters. If FALSE, saved parameters will be reloaded
force_tuning = F  # force tuning even if previous tuned parameters are already stored
tuning_strategy = 'bayes'  # 'bayes' or 'grid'
save_all = T  # save parameters and stats
inn_cross_val_fold = 5  # inner cross-validation fold
out_cross_val_fold = 5  # outer cross-validation fold
tuning_criterion = 'rmse'  # tuning criterion and performance measure
algo_set = c('RandomForest', 'GBM')
df_final = readRDS('./Checkpoints/df_final.rds')
df_add = readRDS('./Checkpoints/df_add.rds')
res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
{
  summary(df_add %>% select(test_variable_set))
  res_thresh_sensitivity_list = res_thresh_sensitivity_best = res_thresh_sensitivity_residual = c()
  comb_tot = length(index_1_set) * length(index_2_set)
  sink(paste0('./Log/Threshold_sensitivity_log_', format(Sys.time(), "%Y-%m-%dT%H%M"), '.txt'), append = F, split = T, type = 'output')
  cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
  for (df_type in df_set){
    for (recov_met in recov_set){
      
      cat('\n\n\n-------------------------------------------------------------------------------------------------------------------------')
      cat('\n                                      data:', df_type, '    method:', recov_met)
      cat('\n-------------------------------------------------------------------------------------------------------------------------\n\n')
      
      for (fit_met in fitting_set){
        
        cat('\n\n\n ----#####  fitting method:  ', fit_met, '  #####----')
        
        var_count = 1
        for (var_target in test_variable_set){
          
          for (algo_type in algo_set){
            
            # reload settings and parameters
            RDS_lab = paste0('./Checkpoints/Threshold_sensitivity/', df_type, '_', recov_met, '_', fit_met, '_', var_target, '_', algo_type, '.rds')
            reload_out = tryCatch({
              readRDS(RDS_lab)
            }, warning = function(e) {
              'warn'
            }, error = function(e) {
              'error'
            }, finally = {
            }, quiet = TRUE)
            if (typeof(reload_out) == "character"){reload_out = list()}
            
            cat('\n\n     ++++ regressing', var_target, ' -', var_count, '/', length(test_variable_set), '  with', algo_type)
            if (algo_type == 'RandomForest'){
              algo_type_work = 'CITree'
            } else {
              algo_type_work = algo_type
            }
            
            # fit model with original regressors
            cat('\n\n       °°°° with original regressors')
            # define dataset df_work_orig
            {
              # check for non-missing for all variables
              data_match = res_ALL_scores %>%
                filter(data == df_type) %>%
                select(country, year) %>%
                unique() %>%
                left_join(df_add, by = c("country", "year")) %>%
                select(c('country', 'year', test_variable_set)) %>%
                group_by(country, year) %>%
                summarize_all(function(x) sum(is.na(x))) %>%
                ungroup() %>%
                mutate(NA_SUM = rowSums(select(., -c(1,2)))) %>%
                mutate(KEEP = ifelse(NA_SUM == 0, 1, 0)) %>%
                select(country, year, KEEP)
              cat('  - removed', sum(data_match$KEEP == 0), 'observations from', nrow(data_match), '-', sum(data_match$KEEP == 1), 'remaining\n')
              
              data_work = res_ALL_scores %>%
                filter(data == df_type) %>%
                filter(method == recov_met) %>%
                filter(algo == fit_met) %>%
                select(country, year) %>%
                unique() %>%
                left_join(df_add, by = c("country", "year")) %>%
                select(c('country', 'year', test_variable_set)) %>%
                left_join(data_match, by = c("country", "year")) %>%
                filter(KEEP == 1) %>%
                select(-KEEP)
              if (sum(is.na(data_work)) != 0){cat('\n ############ data_work: observations with missing not removed')}
              
              df_recover = df_final %>%
                filter(data == df_type) %>%
                filter(method == recov_met) %>%
                filter(year != 'Avg') %>%
                mutate(year = as.numeric(year)) %>%
                select(-method, -data) %>%
                spread(variable, val)
              
              # df with original variables - TARGET is standardized
              df_work_orig = data_work %>%
                select_(.dots = c('country', 'year', var_target)) %>%
                rename_(.dots = setNames(var_target, 'TARGET')) %>%
                left_join(df_recover, by = c("country", "year")) %>%
                mutate(TARGET = scale(TARGET))
              if (sum(is.na(df_work_orig)) != 0){cat('\n ############ df_work_orig: observations with missing not removed')}
            }
            
            regr_orig = threshold_sensitivity_fit(flag_tuning, force_tuning, save_all, inn_cross_val_fold, out_cross_val_fold,
                                                  tuning_criterion, tuning_strategy, df_type, recov_met, fit_met, algo_type, algo_type_work = algo_type,
                                                  var_target, reload_out, index_1_set, index_2_set,
                                                  res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
                                                  regr_to_test = 'original', 
                                                  index_1_thresh = NULL, 
                                                  index_2_thresh = NULL, 
                                                  df_work_orig = df_work_orig, 
                                                  df_work_index = NULL,
                                                  df_work_rank = NULL,
                                                  df_work_raw = NULL)
            
            res_thresh_sensitivity_list = regr_orig$res_thresh_sensitivity_list
            res_thresh_sensitivity_best = regr_orig$res_thresh_sensitivity_best
            res_thresh_sensitivity_residual = regr_orig$res_thresh_sensitivity_residual
            reload_out = regr_orig$reload_out
            
            # fit model with index regressors
            comb_count = 1
            for (index_1_thresh in index_1_set){
              for (index_2_thresh in index_2_set){
                
                cat('\n\n       °°°° with index regressors, thresholds:', index_1_thresh, 'and', index_2_thresh, ' -', round(comb_count / comb_tot * 100, 2), '%')
                # define dataset df_work_index
                {
                  # df with binary index - TARGET is standardized
                  df_work_index = data_work %>%
                    select_(.dots = c('country', 'year', var_target)) %>%
                    rename_(.dots = setNames(var_target, 'TARGET')) %>%
                    left_join(res_ALL_scores %>%
                                filter(data == df_type) %>%
                                filter(method == recov_met) %>%
                                filter(algo == fit_met) %>%
                                mutate(factor = paste0('Index_', factor)) %>%
                                spread(factor, index) %>%
                                select(-method, -data, -algo, -family, -starts_with('Explain')) %>%
                                mutate(Index_1 = ifelse(Index_1 >= index_1_thresh, 1, 0),
                                       Index_2 = ifelse(Index_2 >= index_2_thresh, 1, 0)),
                              by = c("country", "year")) %>%
                    mutate(TARGET = scale(TARGET))
                  if (sum(is.na(df_work_index)) != 0){cat('\n ############ df_work_index: observations with missing not removed')}
                }
                
                regr_index = threshold_sensitivity_fit(flag_tuning, force_tuning, save_all, inn_cross_val_fold, out_cross_val_fold,
                                                       tuning_criterion, tuning_strategy, df_type, recov_met, fit_met, algo_type, algo_type_work,
                                                       var_target, reload_out, index_1_set, index_2_set,
                                                       res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
                                                       regr_to_test = 'index', 
                                                       index_1_thresh = index_1_thresh, 
                                                       index_2_thresh = index_2_thresh, 
                                                       df_work_orig = NULL, 
                                                       df_work_index = df_work_index,
                                                       df_work_rank = NULL,
                                                       df_work_raw = NULL)
                
                res_thresh_sensitivity_list = regr_index$res_thresh_sensitivity_list
                res_thresh_sensitivity_best = regr_index$res_thresh_sensitivity_best
                res_thresh_sensitivity_residual = regr_index$res_thresh_sensitivity_residual
                reload_out = regr_index$reload_out
                comb_count = comb_count + 1
                
              } # index_2_thresh
            } # index_1_thresh
            
            # fit model with ranked index regressors
            cat('\n\n       °°°° with ranked index regressors')
            # define dataset df_work_rank 
            {
              # df with binary index - TARGET is standardized
              df_work_rank = data_work %>%
                select_(.dots = c('country', 'year', var_target)) %>%
                rename_(.dots = setNames(var_target, 'TARGET')) %>%
                left_join(res_ALL_scores %>%
                            filter(data == df_type) %>%
                            filter(method == recov_met) %>%
                            filter(algo == fit_met) %>%
                            mutate(factor = paste0('Index_', factor)) %>%
                            spread(factor, index) %>%
                            select(-method, -data, -algo, -family, -starts_with('Explain')) %>%
                            mutate(Index_1 = cut(Index_1, breaks = c(-Inf,index_1_set,Inf)),
                                   Index_2 = cut(Index_2, breaks = c(-Inf,index_2_set,Inf))),
                          by = c("country", "year")) %>%
                mutate(TARGET = scale(TARGET))
              if (sum(is.na(df_work_rank)) != 0){cat('\n ############ df_work_rank: observations with missing not removed')}
            }
            
            regr_rank = threshold_sensitivity_fit(flag_tuning, force_tuning, save_all, inn_cross_val_fold, out_cross_val_fold,
                                                  tuning_criterion, tuning_strategy, df_type, recov_met, fit_met, algo_type, algo_type_work = algo_type,
                                                  var_target, reload_out, index_1_set, index_2_set,
                                                  res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
                                                  regr_to_test = 'rank_index', 
                                                  index_1_thresh = index_1_thresh, 
                                                  index_2_thresh = index_2_thresh, 
                                                  df_work_orig = NULL, 
                                                  df_work_index = NULL,
                                                  df_work_rank = df_work_rank,
                                                  df_work_raw = NULL)
            
            res_thresh_sensitivity_list = regr_rank$res_thresh_sensitivity_list
            res_thresh_sensitivity_best = regr_rank$res_thresh_sensitivity_best
            res_thresh_sensitivity_residual = regr_rank$res_thresh_sensitivity_residual
            reload_out = regr_rank$reload_out
            
            # fit model with raw index regressors
            cat('\n\n       °°°° with raw index regressors')
            # define dataset df_work_raw 
            {
              # df with binary index - TARGET is standardized
              df_work_raw = data_work %>%
                select_(.dots = c('country', 'year', var_target)) %>%
                rename_(.dots = setNames(var_target, 'TARGET')) %>%
                left_join(res_ALL_scores %>%
                            filter(data == df_type) %>%
                            filter(method == recov_met) %>%
                            filter(algo == fit_met) %>%
                            mutate(factor = paste0('Index_', factor)) %>%
                            spread(factor, index) %>%
                            select(-method, -data, -algo, -family, -starts_with('Explain')),
                          by = c("country", "year")) %>%
                mutate(TARGET = scale(TARGET))
              if (sum(is.na(df_work_raw)) != 0){cat('\n ############ df_work_raw: observations with missing not removed')}
            }
            
            regr_raw = threshold_sensitivity_fit(flag_tuning, force_tuning, save_all, inn_cross_val_fold, out_cross_val_fold,
                                                 tuning_criterion, tuning_strategy, df_type, recov_met, fit_met, algo_type, algo_type_work = algo_type,
                                                 var_target, reload_out, index_1_set, index_2_set,
                                                 res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
                                                 regr_to_test = 'raw_index', 
                                                 index_1_thresh = index_1_thresh, 
                                                 index_2_thresh = index_2_thresh, 
                                                 df_work_orig = NULL, 
                                                 df_work_index = NULL,
                                                 df_work_rank = NULL,
                                                 df_work_raw = df_work_raw)
            
            res_thresh_sensitivity_list = regr_raw$res_thresh_sensitivity_list
            res_thresh_sensitivity_best = regr_raw$res_thresh_sensitivity_best
            res_thresh_sensitivity_residual = regr_raw$res_thresh_sensitivity_residual
            reload_out = regr_raw$reload_out
            
            
            if (save_all){
              cat('\n --- saving RDS')
              saveRDS(reload_out, RDS_lab)
            }
            
          } # algo_type
          var_count = var_count + 1
        } # var_target
      } # fit_met
    } # recov_met
  } # df_type
  cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
  sink()
  
  if (save_all){
    if (nrow(unique(res_thresh_sensitivity_residual)) != nrow(res_thresh_sensitivity_residual)){cat('\n\n\n\n\n ##################### duplicates in res_thresh_sensitivity_residual\n\n\n\n')}
    saveRDS(res_thresh_sensitivity_residual, './Checkpoints/res_thresh_sensitivity_residual.rds')
    saveRDS(res_thresh_sensitivity_best, './Checkpoints/res_thresh_sensitivity_best.rds')
    write.table(res_thresh_sensitivity_list, "./Results/3_threshold_sensitivity_performance_list.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    write.table(res_thresh_sensitivity_best, "./Results/3_threshold_sensitivity_performance_best.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  }
}



### evaluate outlier stability to different index thresholds
res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
res_thresh_sensitivity_residual = readRDS('./Checkpoints/res_thresh_sensitivity_residual.rds')
{
  res_outlier_stability = c()
  for (df_type in unique(res_thresh_sensitivity_residual$data)){
    for (recov_met in unique(res_thresh_sensitivity_residual$method)){
      for (fit_met in unique(res_thresh_sensitivity_residual$fit_method)){
        for (var_target in unique(res_thresh_sensitivity_residual$target_var)){
          
          data_t1 = res_thresh_sensitivity_residual %>%
            filter(data == df_type) %>%
            filter(method == recov_met) %>%
            filter(fit_method == fit_met) %>%
            filter(target_var == var_target)
          
          
          ID_code = data_t1 %>%
            select(country, year) %>%
            unique() %>%
            mutate(ID = paste0(country, year)) %>%
            mutate(ID_NUM = 1:n()) %>%
            select(ID, ID_NUM)
          
          # useless because it cannot be compared with R^2 of target variable
          # Expl_var = res_ALL_scores %>%
          #   filter(data == df_type) %>%
          #   filter(method == recov_met) %>%
          #   filter(algo == fit_met) %>%
          #   select(starts_with('Explain')) %>%
          #   unique() %>%
          #   mutate(Explain_var_avg = mean(Explain_var),
          #          Explain_var_std = sd(Explain_var),
          #          Explain_var_95_avg = mean(Explain_var_95),
          #          Explain_var_95_std = sd(Explain_var_95),
          #          Explain_var_99_avg = mean(Explain_var_99),
          #          Explain_var_99_std = sd(Explain_var_99)) %>%
          #   select(matches('avg|std')) %>%
          #   unique()
          
          outlier_list = c()
          for (reg_type in unique(data_t1$regressor_type)){
            
            data_t2 = data_t1 %>%
              filter(regressor_type == reg_type)
            
            for (algo_m in unique(data_t2$algo_main)){
              for (perf in unique(data_t2$perform)){
                
                data_t3 = data_t2 %>%
                  filter(algo_main == algo_m) %>%
                  filter(perform == perf)
                
                
                for (algo_type in unique(data_t3$algo)){
                  
                  data = data_t3 %>%
                    filter(algo == algo_type)
                  tot_comb = nrow(data %>%
                                    select(index_1_thresh, index_2_thresh) %>%
                                    unique())
                  
                  for (ind1 in unique(data$index_1_thresh)){
                    for (ind2 in unique(data$index_2_thresh)){
                      
                      data_list = data %>%
                        filter(index_1_thresh == ind1) %>%
                        filter(index_2_thresh == ind2) %>%
                        select(country, year, residual, truth) %>%
                        mutate(ID = paste0(country, year),
                               abs_err = abs(residual),
                               abs_perc_err = abs(residual / truth),
                               TSS = (truth - mean(truth)) ^ 2,
                               RSS = residual ^ 2)
                      
                      # evaluate overall performance
                      RSS_val = data_list$RSS
                      TSS_val = data_list$TSS
                      RSS = sum(RSS_val)
                      TSS = sum(TSS_val)
                      ind_95 = RSS_val <= quantile(RSS_val, 0.95)
                      RSS_95 = sum(RSS_val[ind_95])
                      TSS_95 = sum(TSS_val[ind_95])
                      ind_99 = RSS_val <= quantile(RSS_val, 0.99)
                      RSS_99 = sum(RSS_val[ind_99])
                      TSS_99 = sum(TSS_val[ind_99])
                      Explain_var = 1 - RSS / TSS
                      Explain_var_95 = 1 - RSS_95 / TSS_95
                      Explain_var_99 = 1 - RSS_99 / TSS_99
                      
                      for (outlier_measure in c('abs_perc_err', 'abs_err')){
                        obs_list = ID_code %>%
                          left_join(data_list %>%
                                      select_(.dots = c(outlier_measure, 'ID')), by = "ID") %>%
                          arrange(ID_NUM) %>%
                          rename_('measure' = outlier_measure)
                        obs_list = obs_list %>%
                          left_join(rosnerTest(obs_list$measure, k = round(0.8 * nrow(obs_list)), warn = F)$all.stats %>%
                                      select(Obs.Num, Outlier),
                                    by = c('ID_NUM' = 'Obs.Num')) %>%
                          mutate(Outlier = replace(Outlier, is.na(Outlier), FALSE)) %>%
                          arrange(ID_NUM)
                        
                        outlier_list = outlier_list %>% bind_rows(
                          data.frame(algo_main = algo_m,
                                     reg_type = reg_type,
                                     algo_type = algo_type,
                                     algo_perf = perf,
                                     measure = outlier_measure,
                                     Ind1 = ind1, Ind2 = ind2,
                                     obs_list %>% select(ID_NUM, Outlier),
                                     tot_outlier = sum(obs_list$Outlier), stringsAsFactors = F,
                                     mean_measure = mean(obs_list$measure),
                                     mean_measure_95 = mean(sort(obs_list$measure)[1:sum(ind_95)]),
                                     mean_measure_99 = mean(sort(obs_list$measure)[1:sum(ind_99)]),
                                     Model_Explain_var = Explain_var,
                                     Model_Explain_var_95 = Explain_var_95,
                                     Model_Explain_var_99 = Explain_var_99)
                        )
                      } # outlier_measure
                    } # ind2
                  } # ind1
                } # algo_type
              } # perf
            } # algo_m
          } # reg_type
          
          stability_ref = outlier_list %>%
            group_by(algo_main, reg_type, algo_type, algo_perf, measure, ID_NUM) %>%
            summarise(occurrence = sum(Outlier)) %>%
            ungroup()
          # filter(occurrence != 0)
          
          res_outlier_stability = res_outlier_stability %>% bind_rows(
            cbind(data.frame(data = df_type,
                             method = recov_met,
                             fit_method = fit_met,
                             target_var = var_target, stringsAsFactors = F),
                  outlier_list %>%
                    left_join(stability_ref, by = c("algo_main", "reg_type", "algo_type", "algo_perf", "measure", "ID_NUM")) %>%
                    filter(!is.na(occurrence)) %>%
                    group_by(algo_main, reg_type, algo_type, algo_perf, measure, Ind1, Ind2) %>%
                    summarize(maxAll = sum(Outlier == T & occurrence == tot_comb), # common outliers for all thresholds
                              max1 = sum(Outlier == T & (occurrence == tot_comb - 1)), # common outliers for (all-1) thresholds
                              max2 = sum(Outlier == T & occurrence == tot_comb - 2),
                              max3 = sum(Outlier == T & occurrence == tot_comb - 3),
                              max4 = sum(Outlier == T & occurrence == tot_comb - 4),
                              tot_outlier = unique(tot_outlier),
                              mean_measure = unique(mean_measure),
                              mean_measure_99 = unique(mean_measure_99),
                              mean_measure_95 = unique(mean_measure_95),
                              Model_Explain_var = unique(Model_Explain_var),
                              Model_Explain_var_99 = unique(Model_Explain_var_99),
                              Model_Explain_var_95 = unique(Model_Explain_var_95)) %>%
                    ungroup()
                  #Expl_var %>% setNames(paste0('Theoretical_', names(.)))
            )
          )
        } # target_var
      } # fit_met
    } # recov_met
  } # df_type
  saveRDS(res_outlier_stability, './Checkpoints/res_outlier_stability.rds')
}  



### plot comparison of outlier stability and regression model performance (3D plot)
res_outlier_stability = readRDS('./Checkpoints/res_outlier_stability.rds')
res_thresh_sensitivity_best = readRDS('./Checkpoints/res_thresh_sensitivity_best.rds')
{  
  res_outlier_stability_work = res_outlier_stability %>%
    mutate(Ind1 = ifelse(reg_type != 'index', NA, Ind1),
           Ind2 = ifelse(reg_type != 'index', NA, Ind2)) %>%
    mutate(Ind1 = as.numeric(Ind1),
           Ind2 = as.numeric(Ind2),
           maxAll = maxAll / tot_outlier,
           max1 = max1 / tot_outlier,
           max2 = max2 / tot_outlier,
           max3 = max3 / tot_outlier,
           max4 = max4 / tot_outlier) %>%
    mutate(outlier_score = (maxAll + 0.8 * max1 + 0.7 * max2 + 0.6 * max3 + 0.5 * max4) * 100) %>%
    mutate(outlier_score = ifelse(is.na(outlier_score), 0, outlier_score)) %>%
    mutate(mean_measure = ifelse(measure == 'abs_perc_err', mean_measure * 100, mean_measure),
           mean_measure_95 = ifelse(measure == 'abs_perc_err', mean_measure_95 * 100, mean_measure_95),
           mean_measure_99 = ifelse(measure == 'abs_perc_err', mean_measure_99 * 100, mean_measure_99),
           Model_Explain_var = Model_Explain_var * 100,
           Model_Explain_var_99 = Model_Explain_var_99 * 100,
           Model_Explain_var_95 = Model_Explain_var_95 * 100)
  z_lim = range(res_outlier_stability_work$outlier_score)
  z_lim_meas = res_outlier_stability_work %>%
    group_by(measure) %>%
    summarise(max = max(mean_measure)) %>%
    ungroup()
  z_lim_R2 = c(min(res_outlier_stability_work %>% select(starts_with('Model_Expl'))), 100)
  max_tot_outlier = max(res_outlier_stability_work$tot_outlier)
  upper_bin_max_tot_out = 150  # threshold to group all outlier above the value
  
  for (df_type in unique(res_outlier_stability_work$data)){
    for (recov_met in unique(res_outlier_stability_work$method)){
      for(var_target in unique(res_outlier_stability_work$target_var)){
        
        cat('\nevaluating: ', df_type, recov_met, var_target)
        
        set1 = res_outlier_stability %>%
          filter(data == df_type) %>%
          filter(method == recov_met) %>%
          filter(target_var == var_target) %>%
          select(fit_method, algo_main) %>%
          unique()
        
        p_row = c()
        for (fit_met in unique(set1$fit_method)){
          for (algo_m in unique(set1$algo_main)){
            
            best_perf = res_thresh_sensitivity_best %>%
              filter(data == df_type) %>%
              filter(method == recov_met) %>%
              filter(target_var == var_target) %>%
              filter(fit_method == fit_met) %>%
              filter(algo_main == algo_m) %>%
              group_by(algo_main, regressor_type, PERF) %>%
              summarize(train_avg = round(mean(TRAIN_mean), 2),
                        train_std = round(sd(TRAIN_mean), 2),
                        test_avg = round(mean(TEST_mean), 2),
                        test_std = round(sd(TEST_mean), 2)) %>%
              ungroup() %>%
              mutate(train_lab = paste0(train_avg, ifelse(is.na(train_std), '', paste0('±', train_std))),
                     test_lab = paste0(test_avg, ifelse(is.na(test_std), '', paste0('±', test_std)))) %>%
              mutate(label = paste0(' -', ifelse(regressor_type == 'index', 'index (avg)', regressor_type), ' Train: ', train_lab, ' Test: ', test_lab)) %>%
              arrange(regressor_type)
            
            # row header
            image_lab = image_graph(res = 100, width = 570, height = 500, clip = F)
            plot(ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + xlim(0, 10) + ylim(2, 3.5) +
                   annotate("text", x = 0, y = 3, label = paste(fit_met,
                                                                '\n\nModel:', algo_m,
                                                                '\nPerformance:', unique(best_perf$PERF),
                                                                '\nRegressor:',
                                                                paste0(c('', best_perf$label), collapse = '\n')),
                            cex = 7,
                            hjust = 0, vjust = 0.7) + 
                   theme_bw() +
                   theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
                          plot.margin=unit(c(0,0.4,0,0.4),"cm"))
            )
            dev.off()
            
            # plot
            data = res_outlier_stability_work %>%
              filter(data == df_type) %>%
              filter(method == recov_met) %>%
              filter(fit_method == fit_met) %>%
              filter(target_var == var_target) %>%
              filter(algo_main == algo_m)
            
            p_hist = p_surf = c()
            for (meas in unique(data$measure)){
              
              # plot outlier distribution (histogram)
              data_t = data %>%
                filter(measure == meas)
              
              data_hist = data_t %>%
                filter(reg_type == 'index')
              
              z_mat = data_hist %>%
                select(Ind1, Ind2, outlier_score) %>%
                # mutate(outlier_score = ifelse(outlier_score <= 0, NA, outlier_score)) %>%
                spread(Ind2, outlier_score) %>%
                arrange(Ind1) %>%
                select(-Ind1) %>%
                as.matrix()
              
              z_color = data_hist %>%
                select(Ind1, Ind2, maxAll) %>%
                # mutate(maxAll = ifelse(maxAll <= 0, NA, maxAll)) %>%
                spread(Ind2, maxAll) %>%
                arrange(Ind1) %>%
                select(-Ind1) %>%
                as.matrix()
              
              floor_tot_outlier = data_hist %>%
                select(Ind1, Ind2, tot_outlier) %>%
                mutate(tot_outlier = ifelse(tot_outlier >= upper_bin_max_tot_out, upper_bin_max_tot_out * 2, tot_outlier)) %>%
                spread(Ind2, tot_outlier) %>%
                arrange(Ind1) %>%
                select(-Ind1) %>%
                as.matrix()
              
              hist = image_graph(res = 100, width = 800, height = 600, clip = F)
              hist3D(z = z_mat, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                     border = "black",
                     xlab = "index 1", ylab = "index 2", zlab = "shared outlier (%)",
                     main = paste0('Outlier distribution for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                     shade = 0,
                     phi = 20,  theta = -50,
                     ticktype = "detailed",
                     nticks = uniqueN(data_hist$Ind1),
                     space = 0.65,
                     alpha = 0.8,
                     bty = 'b2',
                     zlim = z_lim,
                     col = ramp.col (col = c("red", "blue3"), n = 100),
                     colvar = z_color,
                     image = list(z = floor_tot_outlier, col = ramp.col (col = c("grey", "black"), n = 100, alpha = 0.5))
              )
              dev.off()
              hist = image_crop(hist, geometry_area(width = 500, height = 550, x_off = 120, y_off = 0))
              p_hist = c(p_hist, hist)
              
              
              # plot index stability according to meas (surface)
              z_mat = data_hist %>%
                select(Ind1, Ind2, mean_measure) %>%
                # mutate(mean_measure = ifelse(mean_measure <= 0, NA, mean_measure)) %>%
                spread(Ind2, mean_measure) %>%
                arrange(Ind1) %>%
                select(-Ind1) %>%
                as.matrix()
              
              z_mat_95 = data_hist %>%
                select(Ind1, Ind2, mean_measure_95) %>%
                # mutate(mean_measure_95 = ifelse(mean_measure_95 <= 0, NA, mean_measure_95)) %>%
                spread(Ind2, mean_measure_95) %>%
                arrange(Ind1) %>%
                select(-Ind1) %>%
                as.matrix()
              
              surf = image_graph(res = 100, width = 900, height = 600, clip = F)
              persp3D(z = z_mat, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                      facets = T, curtain = F,
                      xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error (%)', 'not found')),
                      main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                      col = 'darkslategrey', border = 'black',
                      shade = 0,
                      phi = 20,  theta = -50,
                      ticktype = "detailed",
                      nticks = uniqueN(data_hist$Ind1),
                      zlim = c(0, unlist(z_lim_meas %>% filter(measure == meas) %>% select(max))),
                      bty = 'b2'
              )
              # text3D(max(data_hist$Ind1), min(data_hist$Ind2), mean(z_mat, na.rm = T), paste0(' index ', round(mean(z_mat, na.rm = T), 2), '-', round(mean(z_mat_95, na.rm = T), 2)), colvar = NULL, add = T)
              # text3D(max(data_hist$Ind1), min(data_hist$Ind2), mean(z_mat, na.rm = T), ' index', colvar = NULL, add = T)
              persp3D(z = z_mat_95, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                      facets = T, curtain = F,
                      xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error (%)', 'not found')),
                      main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                      col = 'darkgrey', border = 'black',
                      shade = 0,
                      phi = 20,  theta = -50,
                      ticktype = "detailed",
                      nticks = uniqueN(data_hist$Ind1),
                      zlim = c(0, unlist(z_lim_meas %>% filter(measure == meas) %>% select(max))),
                      bty = 'b2',
                      add = T
              )
              # text3D(max(data_hist$Ind1), min(data_hist$Ind2), mean(z_mat_95, na.rm = T), ' 95% of index', colvar = NULL, add = T)
              label_list = data.frame(val = c(mean(z_mat, na.rm = T), mean(z_mat_95, na.rm = T)), label = c(' index', ' 95% of index'), stringsAsFactors = F)
              
              data_surf = data_t %>%
                filter(reg_type != 'index') %>%
                mutate(color = c('blue', 'coral', 'chartreuse'),
                       label = paste0(reg_type, ' ', round(mean_measure, 2), '-', round(mean_measure_95, 2)))
              for (rr in unique(data_surf$reg_type)){
                r_col = data_surf %>% filter(reg_type == rr) %>% select(color) %>% unlist() %>% setNames(NULL)
                r_val = data_surf %>% filter(reg_type == rr) %>% select(mean_measure) %>% unlist() %>% setNames(NULL)
                r_lab = data_surf %>% filter(reg_type == rr) %>% select(label) %>% unlist() %>% setNames(NULL)
                persp3D(z = z_mat_95, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                        facets = T, curtain = F,
                        xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error (%)', 'not found')),
                        main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                        col = 'darkgrey', border = 'black',
                        shade = 0,
                        phi = 20,  theta = -50,
                        ticktype = "detailed",
                        nticks = uniqueN(data_hist$Ind1),
                        zlim = c(0, unlist(z_lim_meas %>% filter(measure == meas) %>% select(max))),
                        bty = 'b2',
                        add = T,
                        image = list(side = r_val, facets = T, col = r_col, border = 'black', z = z_mat)
                )
                # text3D(max(data_hist$Ind1), min(data_hist$Ind2), r_val, r_lab, colvar = NULL, add = T)
                label_list = label_list %>% bind_rows(
                  data.frame(val = r_val, label = r_lab, stringsAsFactors = F)
                )
              }
              label_list = space_label(label_list, 0.07 * unlist(z_lim_meas %>% filter(measure == meas) %>% select(max)))
              text3D(rep(max(data_hist$Ind1), nrow(label_list)), rep(min(data_hist$Ind2), nrow(label_list)), label_list$val, label_list$label, colvar = NULL, add = T)
              dev.off()
              surf = image_crop(surf, geometry_area(width = 650, height = 550, x_off = 170, y_off = 0))
              p_surf = c(p_surf, surf)
            } # meas
            
            # plot Explained Variance (R^2)
            data_R = data %>%
              select(reg_type, Ind1, Ind2, starts_with('Model_Explain')) %>%
              unique()
            data_R_ind = data_R %>%
              filter(reg_type == 'index')
            
            z_mat = data_R_ind %>%
              select(Ind1, Ind2, Model_Explain_var) %>%
              # mutate(Model_Explain_var = ifelse(Model_Explain_var <= 0, NA, Model_Explain_var)) %>%
              spread(Ind2, Model_Explain_var) %>%
              arrange(Ind1) %>%
              select(-Ind1) %>%
              as.matrix()
            
            z_mat_95 = data_R_ind %>%
              select(Ind1, Ind2, Model_Explain_var_95) %>%
              # mutate(Model_Explain_var_95 = ifelse(Model_Explain_var_95 <= 0, NA, Model_Explain_var_95)) %>%
              spread(Ind2, Model_Explain_var_95) %>%
              arrange(Ind1) %>%
              select(-Ind1) %>%
              as.matrix()  
            
            p_R2 = image_graph(res = 100, width = 900, height = 600, clip = F)
            persp3D(z = z_mat, x = sort(unique(data_R_ind$Ind1)), y = sort(unique(data_R_ind$Ind2)),
                    facets = T, curtain = F,
                    xlab = "index 1", ylab = "index 2", zlab = 'R^2',
                    main = 'Explained Variance stability ',
                    col = 'darkslategrey', border = 'black',
                    shade = 0,
                    phi = 20,  theta = -50,
                    ticktype = "detailed",
                    nticks = uniqueN(data_R_ind$Ind1),
                    zlim = z_lim_R2,
                    bty = 'b2'
            )
            # text3D(max(data_R_ind$Ind1), min(data_R_ind$Ind2), mean(z_mat, na.rm = T), ' index', colvar = NULL, add = T)
            persp3D(z = z_mat_95, x = sort(unique(data_R_ind$Ind1)), y = sort(unique(data_R_ind$Ind2)),
                    facets = T, curtain = F,
                    xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error', 'not found')),
                    main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                    col = 'darkgrey', border = 'black',
                    shade = 0,
                    phi = 20,  theta = -50,
                    ticktype = "detailed",
                    nticks = uniqueN(data_R_ind$Ind1),
                    zlim = z_lim_R2,
                    bty = 'b2',
                    add = T
            )
            # text3D(max(data_R_ind$Ind1), min(data_R_ind$Ind2), mean(z_mat_95, na.rm = T), ' 95% of index', colvar = NULL, add = T)
            label_list = data.frame(val = c(mean(z_mat, na.rm = T), mean(z_mat_95, na.rm = T)), label = c(' index', ' 95% of index'), stringsAsFactors = F)
            
            data_R_other = data_R %>%
              filter(reg_type != 'index') %>%
              mutate(color = c('blue', 'coral', 'chartreuse'),
                     label = paste0(reg_type, ' ', round(Model_Explain_var, 2), '-', round(Model_Explain_var_95, 2)))
            for (rr in unique(data_R_other$reg_type)){
              r_col = data_R_other %>% filter(reg_type == rr) %>% select(color) %>% unlist() %>% setNames(NULL)
              r_val = data_R_other %>% filter(reg_type == rr) %>% select(Model_Explain_var) %>% unlist() %>% setNames(NULL)
              r_lab = data_R_other %>% filter(reg_type == rr) %>% select(label) %>% unlist() %>% setNames(NULL)
              persp3D(z = z_mat_95, x = sort(unique(data_R_ind$Ind1)), y = sort(unique(data_R_ind$Ind2)),
                      facets = T, curtain = F,
                      xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error', 'not found')),
                      main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                      col = 'darkgrey', border = 'black',
                      shade = 0,
                      phi = 20,  theta = -50,
                      ticktype = "detailed",
                      nticks = uniqueN(data_R_ind$Ind1),
                      zlim = z_lim_R2,
                      bty = 'b2',
                      add = T,
                      image = list(side = r_val, facets = T, col = r_col, border = 'black', z = z_mat)
              )
              # text3D(max(data_R_ind$Ind1), min(data_R_ind$Ind2), r_val, r_lab, colvar = NULL, add = T)
              label_list = label_list %>% bind_rows(
                data.frame(val = r_val, label = r_lab, stringsAsFactors = F)
              )
            }
            label_list = space_label(label_list, 0.09 * max(z_lim_R2))
            text3D(rep(max(data_hist$Ind1), nrow(label_list)), rep(min(data_hist$Ind2), nrow(label_list)), label_list$val, label_list$label, colvar = NULL, add = T)
            dev.off()
            p_R2 = image_crop(p_R2, geometry_area(width = 650, height = 550, x_off = 190, y_off = 0))
            
            
            # create legend
            bar_legend = image_graph(res = 100, width = 200, height = image_info(p_hist[[1]])$height / 2, clip = F)
            grid.draw(cowplot::get_legend(
              ggplot(data = data.frame(val = unique(res_outlier_stability_work$maxAll)) %>% mutate(x = 1:n()),
                     aes(x = x, y = val, fill = val)) +
                geom_point() +
                scale_fill_gradientn(colours=ramp.col (col = c("red", "blue3"), n = 100),
                                     breaks=c(0,1),labels=c("0% of combinations", "100% of combinations"),
                                     limits=c(0,1),
                                     name = 'Bar color:\nShare of outliers\n')
            ))
            dev.off()
            
            floor_legend= image_graph(res = 100, width = 200, height = image_info(p_hist[[1]])$height / 2, clip = F)
            grid.draw(cowplot::get_legend(
              ggplot(data = res_outlier_stability_work %>%
                       select(tot_outlier) %>%
                       mutate(tot_outlier = ifelse(tot_outlier >= upper_bin_max_tot_out, upper_bin_max_tot_out * 2, tot_outlier)) %>%
                       unique() %>%
                       mutate(x = 1:n()),
                     aes(x = x, y = tot_outlier, fill = tot_outlier)) +
                geom_point() +
                scale_fill_gradientn(colours=ramp.col (col = c("grey", "black"), n = 100, alpha = 0.5),
                                     breaks=c(0,1),labels=c(min(res_outlier_stability_work$tot_outlier), "300+"),
                                     limits=c(0,1),
                                     name = 'Floor color:\nNumber of total outliers\n')
            ))
            dev.off()
            
            regr_legend = image_graph(res = 100, width = 200, height = image_info(p_surf[[1]])$height, clip = F)
            grid.draw(cowplot::get_legend(
              ggplot(data = data.frame(reg = c('index', '95% of index', 'original (full-95%)', 'rank index (full-95%)', 'raw index (full-95%)')), aes(reg, fill = reg)) + 
                geom_bar() +
                scale_fill_manual(name = 'Regressor:', values = c('darkslategrey', 'darkgrey', 'blue', 'coral', 'chartreuse'))
            ))
            dev.off()
            
            # assemble row plot
            p_row = c(p_row, image_append(c(image_lab, p_hist, image_append(c(bar_legend, floor_legend), stack = T), p_surf, p_R2, regr_legend)))
            
          } # algo_m
        } # fit_met
        
        eval(parse(text=paste0('final_plot = image_append(c(', paste0('p_row[[', 1:length(p_row), ']]', collapse = ','), '), stack = T)')))
        
        png(paste0('./Results/3_threshold_sensitivity_stats_comparison_', df_type, '_', recov_met, '_', var_target, '.png'), width = 30, height = 4 * length(p_row), units = 'in', res=300)
        plot(final_plot)
        dev.off()
      } # var_target
    } # recov_met
  } # df_type
}



### plot binary index changes over time for each country
index_1_thresh = 0  # threshold to split index 1 (first factor)
index_2_thresh = 0  # threshold to split index 2 (second factor)
rectangle_threshold = 0.5  # % of total countries with sign change in each index
fit_met = c('RobPCA', 'DFM_multivar')
res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
{
  for (df_type in unique(res_ALL_scores$data)){
    for (recov_met in unique(res_ALL_scores$method)){
      
      tot_year = res_ALL_scores %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        filter(algo %in% fit_met) %>%
        select(year) %>%
        uniqueN()
      data = res_ALL_scores %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        filter(algo %in% fit_met) %>%
        mutate(index = ifelse(factor == 1, ifelse(index > index_1_thresh, 1, 0), index)) %>%
        mutate(index = ifelse(factor == 2, ifelse(index > index_2_thresh, 1, 0), index)) %>%
        mutate(factor = paste0('Index ', factor)) %>%
        # filter(factor == ind) %>%
        select(algo, factor, country, year, index) %>%
        group_by(algo, factor, country) %>%
        mutate(zero = sum(index == 0)) %>%
        ungroup() %>%
        mutate(prevail = ifelse(zero >= round(tot_year / 2), 0, 1)) %>%
        mutate(index_change = ifelse(index != prevail, 1, 0)) %>%
        filter(index_change == 1) %>%
        # filter(country %in% unique(data$country)[1:20]) %>%
        mutate_at(vars(-year), funs(as.factor(.)))
      
      d_rect = data %>%
        group_by(algo, factor, year) %>%
        summarise(COUNT = n()) %>%
        group_by(algo, factor) %>%
        mutate(MAX = max(COUNT)) %>%
        mutate(PLOT = ifelse(COUNT == MAX | COUNT >= rectangle_threshold * uniqueN(data$country), 1, 0)) %>%
        ungroup() %>%
        filter(PLOT == 1) %>%
        mutate(xmin = year - 0.5,
               xmax = year + 0.5,
               ymin = 0.5,
               ymax = uniqueN(data$country) + 0.5) %>%
        merge_rectagle()
      
      data = data %>%
        left_join(d_rect, by = c("algo", "factor", "year"))
      
      p = ggplot(data, aes(x = year, y = country)) +
        geom_rect(data = data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = 'grey', linetype = 2, color = 'black', size = 1, alpha = 0.01) +
        geom_point(aes(color = index), size = 7, shape = 16) +
        scale_colour_manual(values = c('deepskyblue4', 'chartreuse3'), name = 'Index Value') +
        facet_nested(. ~ algo + factor) +
        scale_x_continuous(breaks = unique(data$year)) +
        scale_y_discrete(limits = rev(levels(data$country))) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, face = 'bold'),
              axis.text.y = element_text(size = 12, face = 'bold', vjust = 0.3),
              axis.title = element_text(size = 15, face = 'bold'),
              panel.background = element_rect(fill = "white", colour = "black"),
              panel.grid.major = element_line(colour = "grey", linetype = 'dashed', size = 1),
              strip.text.x = element_text(size = 15, face = 'bold'),
              strip.background = element_rect(color = "black", size = 1),
              panel.border = element_rect(color = "black", fill = NA, size = 1),
              legend.text = element_text(size = 12, face = 'bold'),
              legend.title = element_text(size = 12, face = 'bold'),
              legend.key = element_rect(fill = "white"),
              plot.title = element_text(size = 25, face = 'bold'),
              plot.subtitle = element_text(size = 20)) +
        ggtitle('Change of index values over time', subtitle = 'Most frequent values are not showed, years with most changes are shaded in grey\n')
      
      png(paste0('./Results/3_', df_type, '_', recov_met, '_binary_index_change_over_time.png'), width = 20, height = 50, units = 'in', res=300)
      grid.draw(p)
      dev.off()
    } # recov_met
  } # df_type
}



### test ranking with other financial index
p_val_tol = 0.01 # p-val tolerance for correlation test
quantile_remov = 0.1 # remove quantile_remov from both side
algo_to_check = c('RobPCA', 'DFM_multivar')
index_sum_set = c('Average', 'Mahalanobis', 'Euclidean')#, 'Geometric')
{
  res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
  
  df_rank = read.csv("./Data/Additional_Variables/IMF Financial market development data_1980-2017.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
    select(Country.Name, year, starts_with('FD')) %>%
    rename(country = Country.Name) %>%
    left_join(read.csv("./Data/Additional_Variables/IMF Financial market development data_1980-2017 - match.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F),
              by = c('country' = 'ORIGINAL')) %>%
    mutate(country = ifelse(is.na(REPLACE), country, REPLACE)) %>%
    select(-REPLACE)
  
  quantile_remov_lab = paste0(round(quantile_remov * 100, 2), '% by both sides')
  res_ranking = c()
  for (df_type in c('Difference')){#unique(res_ALL_scores$data)){
    for (recov_met in c('TENSOR_BF')){#unique(res_ALL_scores$method)){
      for (algo_type in algo_to_check){
        for (index_sum in index_sum_set){
          
          index = res_ALL_scores %>%
            filter(data == df_type) %>%
            filter(method == recov_met) %>%
            filter(algo == algo_type) %>%
            select(country, year, factor, index) %>%
            mutate(factor = paste0('Index', factor)) %>%
            spread(factor, index)
          
          plot_final = c()
          
          # test different reference indexes
          for (ref_ind in setdiff(colnames(df_rank), c('country', 'year'))){
            
            plot_index = ref_duplicates = tot_obs = c()
            
            # evaluate ranking for each year
            for (yr in unique(index$year)){
              
              # the two factors are summarized into one using the Mahalanobis distance, where the covariance matrix is evaluated over all countries for each year
              # or simple average or Euclidean distance from (0,0)
              index_year = index %>%
                filter(year == yr)
              ma_centers = sapply(index_year %>% select(starts_with('Index')), mean)
              ma_covar = cov(index_year %>% select(starts_with('Index')))
              if (index_sum == 'Average'){
                index_year = index_year %>% rowwise() %>% mutate(distance = mean(!!quo(c(Index1, Index2))))
              } else if (index_sum == 'Mahalanobis'){
                index_year$distance = sqrt(mahalanobis(index_year %>% select(starts_with('Index')), ma_centers, ma_covar))
              } else if (index_sum == 'Euclidean'){
                index_year$distance = sqrt(index_year$Index1 ^ 2 + index_year$Index2 ^ 2)
              } else if (index_sum == 'Geometric'){
                index_year$distance = sqrt(index_year$Index1 * index_year$Index2)
              }
              
              index_ref_all = index_year %>%
                left_join(df_rank %>%
                            select(c('country', 'year', ref_ind)) %>%
                            rename(reference = !!as.name(ref_ind)), by = c("country", "year")) %>%
                filter(!is.na(reference))
              
              index_by_dist = index_by_refer = c()
              quantile_list_dist = quantile(index_ref_all$distance, c(0.25, 0.5, 0.75, quantile_remov, 1 - quantile_remov), na.rm = T, type = 7)
              quantile_list_refer = quantile(index_ref_all$reference, c(0.25, 0.5, 0.75, quantile_remov, 1 - quantile_remov), na.rm = T, type = 7)
              subset_lab = c()
              for (quantile_set in c('All data', quantile_remov_lab, '1st Quartile', '2nd Quartile', '3rd Quartile', '4th Quartile')){
                
                if (quantile_set == 'All data'){
                  index_by_dist = index_ref_all
                  index_by_refer = index_ref_all
                  plot_index = plot_index %>%   # used only to plot the index distribution
                    rbind(index_by_dist)
                } else if (quantile_set == '1st Quartile'){
                  index_by_dist = index_ref_all %>%
                    filter(distance <= quantile_list_dist[1])
                  index_by_refer = index_ref_all %>%
                    filter(reference <= quantile_list_refer[1])
                } else if (quantile_set == '2nd Quartile'){
                  index_by_dist = index_ref_all %>%
                    filter(distance <= quantile_list_dist[2]) %>%
                    filter(distance > quantile_list_dist[1])
                  index_by_refer = index_ref_all %>%
                    filter(reference <= quantile_list_refer[2]) %>%
                    filter(reference > quantile_list_refer[1])
                } else if (quantile_set == '3rd Quartile'){
                  index_by_dist = index_ref_all %>%
                    filter(distance <= quantile_list_dist[3]) %>%
                    filter(distance > quantile_list_dist[2])
                  index_by_refer = index_ref_all %>%
                    filter(reference <= quantile_list_refer[3]) %>%
                    filter(reference > quantile_list_refer[2])
                } else if (quantile_set == '4th Quartile'){
                  index_by_dist = index_ref_all %>%
                    filter(distance > quantile_list_dist[3])
                  index_by_refer = index_ref_all %>%
                    filter(reference > quantile_list_refer[3])
                } else if (quantile_set == quantile_remov_lab){
                  index_by_dist = index_ref_all %>%
                    filter(distance < quantile_list_dist[5]) %>%
                    filter(distance > quantile_list_dist[4])
                  index_by_refer = index_ref_all %>%
                    filter(reference < quantile_list_refer[5]) %>%
                    filter(reference > quantile_list_refer[4])
                }
                subset_lab = c(subset_lab, paste0(quantile_set, ' (', nrow(index_by_dist), ' obs)'))
                
                cor_kendall = suppressWarnings(cor.test(index_by_dist$distance, index_by_dist$reference,  method="kendall", exact = T)) # null hypotesis is 0 correlation
                cor_spearman = suppressWarnings(cor.test(index_by_dist$distance, index_by_dist$reference,  method="spearman", exact = T)) # null hypotesis is 0 correlation
                cor_somers = rcorr.cens(index_by_dist$distance, index_by_dist$reference, outx = TRUE)
                kruskal = kruskal.test(distance ~ reference, data = index_by_dist)  # null hypotesis is that the distributions are the same
                matching_countries = length(intersect(index_by_refer$country, index_by_dist$country)) / nrow(index_by_dist)
                
                ref_duplicates = c(ref_duplicates, nrow(index_by_dist) - uniqueN(index_by_dist$reference))
                tot_obs = c(tot_obs, nrow(index_by_dist))
                res_ranking = res_ranking %>%
                  rbind(data.frame(method = recov_met, data = df_type, algo = algo_type, index_summary = index_sum, year = yr,
                                   subset = subset_lab[length(subset_lab)], reference_index = ref_ind,
                                   tot_obs = nrow(index_by_dist),
                                   reference_duplicates = nrow(index_by_dist) - uniqueN(index_by_dist$reference),
                                   index_duplicates = nrow(index_by_dist) - uniqueN(index_by_dist$distance),
                                   kendall_corr = cor_kendall$estimate,
                                   kendall_pVal = cor_kendall$p.value,
                                   kendall_warn = ifelse(cor_kendall$p.value > p_val_tol, '*', ''),
                                   spearman_corr = cor_spearman$estimate,
                                   spearman_pVal = cor_spearman$p.value,
                                   spearman_warn = ifelse(cor_spearman$p.value > p_val_tol, '*', ''),
                                   somers_corr = cor_somers[2],
                                   somers_confidence = cor_somers[3],
                                   somers_warn = ifelse(abs(cor_somers[2]) * p_val_tol < cor_somers[3], '*', ''), # warn if confidence is greater than p_val_tol * somers_corr
                                   kruskalWallis_pVal_high_means_same = kruskal$p.value,
                                   kruskalWallis_warn = ifelse(kruskal$p.value <= p_val_tol, '*', ''),
                                   matching_countries = matching_countries,
                                   matching_countries_warn = '', stringsAsFactors = F)
                  )
              } # quantile_set
            } # year
            
            # plot index distribution
            plot_index = plot_index %>%
              rename(Index_distance = distance,
                     (!!as.name(ref_ind)) := reference) %>%
              gather('index', 'val', -c(year, country)) %>%
              mutate(index = gsub('Index_distance', 'Aggregated distance', index)) %>%
              mutate(index = gsub('Index1', 'Index 1', index)) %>%
              mutate(index = gsub('Index2', 'Index 2', index)) %>%
              mutate(country = as.factor(country),
                     year = as.factor(year),
                     index = factor(index, levels=c('Index 1', 'Index 2', 'Aggregated distance', ref_ind)))
            
            p_index = image_graph(res = 100, width = 1200, height = 600, clip = F)
            suppressMessages(plot(ggplot(plot_index, aes(x = val, y = year, fill = year)) +
                                    geom_density_ridges(scale = 5, alpha = 0.5) +
                                    facet_wrap(~index, ncol = 4, scales = 'free_x') +
                                    scale_fill_brewer(palette = "YlGnBu", guide = guide_legend(reverse = TRUE)) +
                                    ggtitle('Distribution of indexes',
                                            subtitle = paste0('Total observation: ', max(tot_obs), '  Duplicates in ', ref_ind, ': ', paste0(unique(range(ref_duplicates)), collapse = '-'))) +
                                    theme(axis.text = element_text(size = 18),
                                          axis.title=element_blank(),
                                          text = element_text(size=20),
                                          legend.text = element_text(size = 17),
                                          legend.title = element_text(size = 20))))
            dev.off()
            
            # plot correlation
            plot_corr = res_ranking %>%
              filter(method == recov_met) %>%
              filter(data == df_type) %>%
              filter(algo == algo_type) %>%
              filter(index_summary == index_sum) %>%
              filter(reference_index == ref_ind) %>%
              rename(kruskalWallis_corr = kruskalWallis_pVal_high_means_same) %>%
              select(-method, -data, -algo, -reference_index, -tot_obs, -reference_duplicates, -index_duplicates,
                     -ends_with('pVal'), -somers_confidence, -index_summary) %>%
              gather('corr', 'val', -c(subset, year, ends_with('warn'))) %>%
              mutate(warn = '')  %>%
              mutate(warn = ifelse(kendall_warn == '*' & corr == 'kendall_corr', '*', warn)) %>%
              mutate(warn = ifelse(spearman_warn == '*' & corr == 'spearman_corr', '*', warn)) %>%
              mutate(warn = ifelse(somers_warn == '*' & corr == 'somers_corr', '*', warn)) %>%
              mutate(warn = ifelse(kruskalWallis_warn == '*' & corr == 'kruskalWallis_corr', '*', warn)) %>%
              mutate(corr = gsub('_corr', '', corr)) %>%
              mutate(corr = gsub('kruskalWallis', 'kruskal wallis p-Val \n(high means same distribution)', corr)) %>%
              mutate(corr = gsub('matching_countries', '% of matching countries', corr)) %>%
              select(-ends_with('_warn')) %>%
              mutate(warn = as.factor(warn),
                     subset = factor(subset, levels = subset_lab),
                     corr = factor(corr, levels = c('kendall', 'spearman', 'somers', 'kruskal wallis p-Val \n(high means same distribution)', '% of matching countries')))
            
            p_corr = image_graph(res = 100, width = 1200 * 2, height = 600, clip = F)
            plot(ggplot(plot_corr %>%
                          mutate(split_by_plot = subset,
                                 split_by_line = corr), aes(x = year, y = val)) +
                   geom_hline(yintercept=0) +
                   geom_line(aes(colour = split_by_line), size = 1) +
                   geom_text(aes(label=warn, color=split_by_line), size=10, show.legend = FALSE) +
                   scale_x_continuous(breaks = unique(plot_corr$year)) +
                   # scale_color_manual(values=c('black', 'chocolate4', 'chocolate3', 'chocolate1', 'darkgoldenrod1', 'blue')) + # use if split_by_plot = corr
                   scale_color_manual(values=c('black', 'darkgoldenrod1', 'blue', 'red', 'green3')) + # use if split_by_plot = subset
                   ylim(-1 , 1) +
                   facet_grid(~split_by_plot, scales = "free_x") +
                   ggtitle('Rank correlations',
                           subtitle = paste0('* means p-Val > ', p_val_tol, '  - Subsets extracted by Aggregated distance')) +
                   theme(axis.text = element_text(size = 18),
                         axis.text.x = element_text(angle = 90),
                         axis.title=element_blank(),
                         text = element_text(size=20),
                         title = element_text(size=20),
                         panel.background = element_rect(fill = "white", colour = "black"),
                         panel.grid.major = element_line(colour = "grey", linetype = 'dashed', size = 1),
                         legend.text = element_text(size = 17),
                         legend.title = element_text(size = 20)) +
                   guides(colour = guide_legend(override.aes = list(size=1.5))) +
                   labs(color = "Set of data"))
            dev.off()
            
            # row label with reference index name
            image_lab = image_graph(res = 100, width = 200, height = 600, clip = F)
            plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            text(x = 0.5, y = 0.5, ref_ind, cex = 1.6, col = "black", srt = 90)
            box()
            dev.off()
            
            plot_final = c(plot_final, image_append(c(image_lab, p_index, p_corr), stack = F))
            
          } # ref_ind
          
          eval(parse(text=paste0('plot_final = image_append(c(', paste0('plot_final[[', 1:length(plot_final), ']]', collapse = ','), '), stack = T)')))
          
          png(paste0('./Results/4_ranking_power_', df_type, '_', recov_met, '_', algo_type, '_', index_sum, '.png'), width = 20, height = 5 * uniqueN(res_ranking$reference_index), units = 'in', res=300)
          plot(plot_final)
          dev.off()
        } # index_sum
      } # algo_type
    } # recov_met
  } # df_type
  
  write.table(res_ranking, './Results/4_ranking_power_summary.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
}



### tables and figures for Latex
{
  # variable stats
  {
    latex_stats = readRDS('./Checkpoints/stats_var.rds') %>%
      cbind(data.frame(Type = 'FSI', Frequency = 'Yearly', stringsAsFactors = F)) %>%
      select(Type, variable, Frequency, TOT, NAs, MIN, MAX, MEAN, STD, VAR_COEFF) %>%
      mutate(variable = gsub('FSI_', '', variable)) %>%
      mutate(variable = gsub('_', ' ', variable)) %>%
      arrange(variable) %>%
      rename(Variable = variable) %>%
      bind_rows(readRDS('./Checkpoints/add_centroid.rds') %>%
                  select(country, variable, val) %>%
                  unique() %>%
                  mutate(variable = gsub('GEO_lat', 'Latitude', variable)) %>%
                  mutate(variable = gsub('GEO_lon', 'Longitude', variable)) %>%
                  rename(Variable = variable) %>%
                  mutate(val = round(val, 2)) %>%
                  group_by(Variable) %>%
                  summarise(TOT = n(),
                            NAs = sum(is.na(val)),
                            MIN = min(val, na.rm = T),
                            MAX = max(val, na.rm = T),
                            MEAN = mean(val, na.rm = T),
                            STD = sd(val, na.rm = T),
                            VAR_COEFF = sd(val, na.rm = T) / abs(mean(val, na.rm = T))) %>%
                  ungroup() %>%
                  cbind(data.frame(Type = 'Geography', Frequency = 'Static', stringsAsFactors = F)) %>%
                  select(Type, Variable, Frequency, everything()) %>%
                  arrange(Variable)) %>%
      bind_rows(readRDS('./Checkpoints/add_Hofstede_or.rds') %>%
                  left_join(data.frame(variable = c("GEO_HOF_idv", "GEO_HOF_ivr", "GEO_HOF_ltowvs", "GEO_HOF_mas", "GEO_HOF_pdi" ,"GEO_HOF_uai"),
                                       Variable = c('Individualism vs collectivism', 'Indulgence vs restraint', 'Long term orientation vs short term normative orientation',
                                                    'Masculinity vs femininity', 'Power distance', 'Uncertainty avoidance'), stringsAsFactors = F), by = "variable") %>%
                  mutate(val = round(val, 2)) %>%
                  group_by(Variable) %>%
                  summarise(TOT = n(),
                            NAs = sum(is.na(val)),
                            MIN = min(val, na.rm = T),
                            MAX = max(val, na.rm = T),
                            MEAN = mean(val, na.rm = T),
                            STD = sd(val, na.rm = T),
                            VAR_COEFF = sd(val, na.rm = T) / abs(mean(val, na.rm = T))) %>%
                  ungroup() %>%
                  cbind(data.frame(Type = 'Hofstede', Frequency = 'Static', stringsAsFactors = F)) %>%
                  select(Type, Variable, Frequency, everything()) %>%
                  arrange(Variable)) %>%
      mutate_if(is.numeric, round, 2) %>%
      setNames(c('Type', 'Variable', 'Frequency', 'Total Observations', 'Missing Values', 'Min', 'Max', 'Mean', 'Standard Deviation', 'Variation Coefficient'))
    write.table(latex_stats, './Paper/Latex_Table_Figure/00_stats.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  }
  
  # scree plot for PCA
  df_type = 'Difference'
  recov_met = 'TENSOR_BF'
  res_PCA_list = readRDS('./Checkpoints/res_PCA_list.rds')
  max_variance = 35
  {
      n_row = length(names(res_PCA_list[[df_type]][[1]][[1]]))    # #_PCA_meth
      n_col = length(names(res_PCA_list[[df_type]][[1]])) - 1     # number of years
      
      row_list = list()
      for (yr in setdiff(names(res_PCA_list[[df_type]][[1]]), 'Avg')){
        i = 1
          for (pc_met in names(res_PCA_list[[df_type]][[1]][[1]])){
            row_list[[yr]][[i]] = ggplotGrob(
              res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]][['scree_plot']]+
                theme(axis.title.y=element_blank(),
                      axis.title.x=element_blank(),
                      plot.title = element_text(size=35),
                      axis.text = element_text(size = 25)) +
                ggtitle(paste0(pc_met, ' - ', yr)) +
                ylim(0, max_variance + 10)
            )
            i = i + 1
          } # pc_met
      } # yr
      
      col_list = list()
      for (i in c(1:n_col)){
        col_list[[i]] = do.call(rbind, c(row_list[[i]], size="last"))
      }
      g = do.call(cbind, c(col_list, size="last"))
      png(paste0('./Paper/Latex_Table_Figure/Scree_plot.png'), width = 35, height = 20, units = 'in', res=300)
      grid.draw(g)
      dev.off()
  }
  
  # PCA stats
  res_PCA_stats = readRDS('./Checkpoints/res_PCA_stats.rds')
  latex_PCA_stats = res_PCA_stats %>%
    filter(method == 'TENSOR_BF') %>%
    filter(data == 'Difference') %>%
    filter(year != 'Avg') %>%
    filter(PC <= 2) %>%
    rename(Method = PCA,
           `Number of PC` = PC) %>%
    group_by(Method, `Number of PC`) %>%
    summarise(`Mean Explained Variance` = paste0(round(mean(Explain_var_Loadings) * 100, 1), '±', round(sd(Explain_var_Loadings) * 100, 1), '%'),
              `Mean R^2` = paste0(round(mean(Explain_var) * 100, 1), '±', round(sd(Explain_var) * 100, 1), '%'),
              `Mean R^2 on 99th` = paste0(round(mean(Explain_var_99) * 100, 1), '±', round(sd(Explain_var_99) * 100, 1), '%'),
              `Mean R^2 on 95th` = paste0(round(mean(Explain_var_95) * 100, 1), '±', round(sd(Explain_var_95) * 100, 1), '%')
    )
  write.table(latex_PCA_stats, './Paper/Latex_Table_Figure/01_PCA_stats.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  # DFM stats
  VAR_alpha = 0.2
  kalm_Q_hat_mode = 'identity'
  res_DFM_stats = readRDS(paste0('./Checkpoints/res_DFM_stats_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  latex_DFM_stats = res_DFM_stats %>%
    filter(method == 'TENSOR_BF') %>%
    filter(data == 'Difference') %>%
    filter(DFM == 'DFM_multivar') %>%
    filter(Total_Factors <= 2) %>%
    rename(Method = DFM,
           `Number of Factors` = Total_Factors) %>%
    mutate(Method = 'DFM') %>%
    group_by(Method, `Number of Factors`) %>%
    summarise(`R^2` = paste0(round(mean(Explain_var) * 100, 1), '%'),
              `R^2 on 99th` = paste0(round(mean(Explain_var_99) * 100, 1), '%'),
              `R^2 on 95th` = paste0(round(mean(Explain_var_95) * 100, 1), '%')
    )
  write.table(latex_DFM_stats, './Paper/Latex_Table_Figure/02_DFM_stats.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  # PCA and DFM cluster full list 2 dimension only
  df_type = 'Difference'
  recov_met = 'TENSOR_BF'
  pc_met = 'RobPCA'
  dfm_met = 'DFM_multivar'
  year_to_plot_sample = '2014'
  {
    res_index_DFM = readRDS('./Checkpoints/res_index_DFM.rds')
    res_DFM_factors = readRDS(paste0('./Checkpoints/res_DFM_factors_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
    res_DFM_loadings = readRDS(paste0('./Checkpoints/res_DFM_loadings_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
    res_DFM_best_model = readRDS(paste0('./Checkpoints/res_DFM_best_model_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
    res_PCA_list = readRDS('./Checkpoints/res_PCA_list.rds')
    res_PCA_loadings = readRDS('./Checkpoints/res_PCA_loadings.rds')
    res_PCA_importance = readRDS('./Checkpoints/res_PCA_importance.rds')
    res_index_PCA = readRDS('./Checkpoints/res_index_PCA.rds')
    index_1_thresh = 0
    index_2_thresh = 0
    load_thresh = 0  # loadings with magnitude below threshold are set to 0 and scores are evaluated always as X * loadings
    leading_var = "FSI_Emb_Capital_to_assets"
    leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
    PC_to_keep = 2
    VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
    kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
    expl_var_to_show=95
    
    # PCA
    load_ref = sort_loading_PCA(res_PCA_loadings, leading_var, leading_sign) %>%
      filter(data == df_type) %>%
      filter(method == recov_met) %>%
      filter(PCA == pc_met) %>%
      filter(PC <= PC_to_keep)
    load_range = range(load_ref$loading)
    
    expl_var_out = c()
    scores_out_raw = scores_out_index = data.frame(COUNTRY = c('INDEX', rownames(res_PCA_list[[df_type]][[recov_met]][[1]][[pc_met]][['pca']][['x']])), stringsAsFactors = F)
    load_out = data.frame(VARIABLE = c('VARIABLE', rownames(res_PCA_list[[df_type]][[recov_met]][[1]][[pc_met]][['pca']][['rotation']])), stringsAsFactors = F)
    year_set = setdiff(sort(unique(load_ref$year)), 'Avg')
    fig_lab = 1
    for (year_split in list(year_set[1:ceiling(length(year_set)/2)], year_set[(ceiling(length(year_set)/2)+1):length(year_set)])){
      
      row_list = list()
      i = 1
      
      for (yr in year_split){
        
        load = load_ref %>%
          filter(year == yr) %>%
          mutate(PC = paste0('I_', PC)) %>%
          spread(PC, loading) %>%
          select(-PCA, -data, -year, -method)
        X = res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]][['pca']][['orig_data']] %>% scale()
        row_order = rownames(res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]][['pca']][['rotation']])
        load = data.frame(variable = row_order, stringsAsFactors = F) %>% left_join(load, by = "variable") %>% select(-variable) %>% as.matrix()
        load[abs(load) < load_thresh] = 0
        # raw scores
        scores_orig = X %*% load %>% round(digits = 2) %>% as.data.frame() %>% mutate(INDEX = rownames(X)) %>% select(INDEX, everything())
        scores = scores_out_raw %>% select(COUNTRY) %>% left_join(colnames(scores_orig) %>% rbind(scores_orig), by = c('COUNTRY' = 'INDEX'))
        colnames(scores) = c('COUNTRY', rep(yr, PC_to_keep))
        scores_out_raw = scores_out_raw %>% cbind(scores[, -1])
        # scores to index by thresholds
        scores_index = scores_orig %>% mutate(I_1 = ifelse(I_1 > index_1_thresh, 1, 0))
        if ('I_2' %in% colnames(scores_orig)){scores_index = scores_index %>% mutate(I_2 = ifelse(I_2 > index_2_thresh, 1, 0))}
        scores_index_mod = scores_out_raw %>% select(COUNTRY) %>% left_join(colnames(scores_index) %>% rbind(scores_index), by = c('COUNTRY' = 'INDEX'))
        colnames(scores_index_mod) = c('COUNTRY', rep(yr, PC_to_keep))
        scores_out_index = scores_out_index %>% cbind(scores_index_mod[, -1])
        
        load_out = load_out %>% cbind(colnames(load) %>% rbind(load) %>% as.data.frame() %>% setNames(rep(yr, PC_to_keep)))
        
        # score plot
        for (pc in 2:min(c(PC_to_keep, 2))){
          
          # variable importance plot
          p_load = plot_loadings(load_list = data.frame(Variable = row_order, Importance = load[, pc]),
                                 load_range = load_range,
                                 dim_to_show = pc)
          
          expl_var = res_PCA_importance %>%
            filter(PCA == pc_met & method == recov_met & data == df_type & year == yr & PC %in% paste0('PC', pc))
          
          expl_var_out = expl_var_out %>% bind_rows(
            c(year = yr, PC = pc, Expl_Var = expl_var$`Cumulative Proportion`)
          )
          
          # score plot
          p = plot_index(index_list = scores_orig,
                         index_1_thresh = index_1_thresh,
                         index_2_thresh = index_2_thresh,
                         dim_to_show = pc,
                         year = yr,
                         expl_var = expl_var$`Cumulative Proportion`)
          
          if (yr == year_to_plot_sample){
            p1 = p + ggtitle('')
            png(paste0('./Paper/Latex_Table_Figure/Index_plot_PCA.png'), width = 10, height = 6, units = 'in', res=200)
            plot(p1)
            dev.off()
          }
          
          row_list[['1']][[i]] = ggplotGrob(p)
          row_list[['2']][[i]] = ggplotGrob(p_load)
        } # pc
        i = i + 1
      } # yr
      col_list = list()
      for (j in 1:length(row_list)){
        col_list[[j]] =  do.call(rbind, c(row_list[[j]], size="last"))
      }
      if (length(col_list) == 2){
        g = cbind(col_list[[1]], col_list[[2]], size="last")
      } else if (length(col_list) > 2){
        g = g = do.call(cbind, c(col_list, size="last"))
      } else {
        g = col_list[[1]]
      }
      png(paste0('./Paper/Latex_Table_Figure/Cluster_full_PCA_', fig_lab, '.png'), width = 20, height = 6 * length(year_split), units = 'in', res=200)
      grid.draw(g)
      dev.off()
      fig_lab = fig_lab + 1
    } # year_split
    
    # DFM
    best_model = res_DFM_best_model %>%
      filter(data == df_type) %>%
      filter(method == recov_met) %>%
      filter(DFM == dfm_met) %>%
      filter(Best_model == 'YES')
    n_factor_best = best_model$Total_Factors
    expl_var_to_extract = ifelse(expl_var_to_show==0, 'Explain_var', paste0('Explain_var_', expl_var_to_show))
    expl_var = unname(unlist(best_model %>% select(expl_var_to_extract)))
    
    factors_ref = res_DFM_factors %>%
      filter(data == df_type) %>%
      filter(method == recov_met) %>%
      filter(DFM == dfm_met) %>%
      filter(Total_Factors == n_factor_best) %>%
      select(-method, -data, -Var_removed, -DFM)
    max_factor = max(factors_ref$Factor)
    
    # evaluate average loading over all countries
    load_ref = sort_loading_DFM(res_DFM_loadings, leading_var, leading_sign) %>%
      filter(data == df_type) %>%
      filter(method == recov_met) %>%
      filter(DFM == dfm_met) %>%
      filter(Total_Factors == n_factor_best) %>%
      filter(Factor <= 2) %>%
      group_by(variable, Factor) %>%
      summarize(load_min = min(loading),
                load_max = max(loading),
                loading = mean(loading)) %>%
      ungroup()
    load_range = c(min(load_ref$load_min), max(load_ref$load_max))
    
    scores_out_raw = scores_out_index = data.frame(COUNTRY = c('INDEX', sort(unique(factors_ref$country))), stringsAsFactors = F)
    load_out = data.frame(VARIABLE = c('VARIABLE', unique(load_ref$variable)), stringsAsFactors = F)
    
    year_set = sort(unique(factors_ref$year))
    fig_lab = 1
    for (year_split in list(year_set[1:ceiling(length(year_set)/2)], year_set[(ceiling(length(year_set)/2)+1):length(year_set)])){
      
      row_list = list()
      i = 1
      
      for (yr in year_split){
        
        factors = factors_ref %>%
          filter(year == yr) %>%
          select(-Total_Factors, - year) %>%
          spread(Factor, val) %>%
          setNames(c('INDEX', paste0('I_', c(1:max_factor))))
        
        load = load_ref %>%
          mutate(Factor = paste0('I_', Factor)) %>%
          select(-load_min, -load_max) %>%
          spread(Factor, loading)
        
        # raw scores
        scores = scores_out_raw %>% select(COUNTRY) %>% left_join(colnames(factors) %>% rbind(factors), by = c('COUNTRY' = 'INDEX'))
        colnames(scores) = c('COUNTRY', rep(yr, max_factor))
        scores_out_raw = scores_out_raw %>% cbind(scores[, -1])
        # scores to index by thresholds
        scores_index = factors %>% mutate(I_1 = ifelse(I_1 > index_1_thresh, 1, 0))
        if ('I_2' %in% colnames(factors)){scores_index = scores_index %>% mutate(I_2 = ifelse(I_2 > index_2_thresh, 1, 0))}
        scores_index_mod = scores_out_raw %>% select(COUNTRY) %>% left_join(colnames(scores_index) %>% rbind(scores_index), by = c('COUNTRY' = 'INDEX'))
        colnames(scores_index_mod) = c('COUNTRY', rep(yr, max_factor))
        scores_out_index = scores_out_index %>% cbind(scores_index_mod[, -1])
        
        load_out = load_out %>% cbind(setdiff(colnames(load), 'variable') %>% rbind(load %>% select(-variable)) %>% as.data.frame() %>% setNames(rep(yr, max_factor)))
        
        for (fact in 2:min(c(max_factor, 2))){
          
          # variable importance plot - is always the same
          p_load = plot_loadings(load_list = data.frame(Variable = load$variable,
                                                        Importance = unlist(load[, paste0('I_', fact)]),
                                                        load_ref %>% filter(Factor == fact) %>% select(load_min, load_max)),
                                 load_range = load_range,
                                 dim_to_show_lab = fact,
                                 add_title = ' - average over countries - same for all years',
                                 err_bar = T)
          
          # plot index
          p = plot_index(index_list = factors,
                         index_1_thresh = index_1_thresh,
                         index_2_thresh = index_2_thresh,
                         dim_to_show = fact,
                         year = yr,
                         expl_var = expl_var,
                         add_title_variance = ifelse(expl_var_to_show==0, '', paste0(' (', expl_var_to_show, 'th percentile)')))
          
          if (yr == year_to_plot_sample){
            p1 = p + ggtitle('')
            png(paste0('./Paper/Latex_Table_Figure/Index_plot_FA.png'), width = 10, height = 6, units = 'in', res=200)
            plot(p1)
            dev.off()
          }
          
          row_list[['1']][[i]] = ggplotGrob(p)
          row_list[['2']][[i]] = ggplotGrob(p_load)
        } # fact
        i = i + 1
      } # yr
      col_list = list()
      for (j in 1:length(row_list)){
        col_list[[j]] =  do.call(rbind, c(row_list[[j]], size="last"))
      }
      if (length(col_list) == 2){
        g = cbind(col_list[[1]], col_list[[2]], size="last")
      } else if (length(col_list) > 2){
        g = g = do.call(cbind, c(col_list, size="last"))
      } else {
        g = col_list[[1]]
      }
      png(paste0('./Paper/Latex_Table_Figure/Cluster_full_FA_', fig_lab, '.png'), width = 20, height = 6 * length(year_split), units = 'in', res=200)
      grid.draw(g)
      dev.off()
      fig_lab = fig_lab + 1
    } # year_split
  }
  
  # Index evolution over time and list of index (both binary and continous)
  df_type = 'Difference'
  recov_met = 'TENSOR_BF'
  algo_to_plot = c('RobPCA', 'DFM_multivar')
  res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
  {
    country_short = read.csv("./Data/country_short_names_mapping.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)
    
    comp_data = res_ALL_scores %>%
      left_join(country_short, by = "country") %>%
      mutate(country = ifelse(is.na(label), country, label)) %>%
      filter(data == df_type) %>%
      filter(method == recov_met)  %>%
      filter(algo %in% algo_to_plot) %>%
      mutate(factor = paste0('Index ', factor),
             color_id = paste0(family, '_', algo)) %>%
      mutate(algo = replace(algo, algo == 'DFM_multivar', 'DFM')) %>%
      select(-label)
    
    
    # create color palette for each family
    grad_bl = c('blue3', 'deepskyblue', 'deepskyblue3')  # bluescale
    grad_rd = c('brown', 'brown1', 'darkorange')  # redscale
    grad_gr = c('chartreuse4', 'chartreuse3', 'chartreuse')  # greenscale
    color_palette_list = comp_data %>%
      select(family, color_id) %>%
      unique() %>%
      arrange(color_id) %>%
      group_by(family) %>%
      summarise(rep = n())
    color_palette = c()
    for (pal in 1:nrow(color_palette_list)){
      if (pal == 1){color_palette = c(color_palette, grad_bl[1:color_palette_list$rep[pal]])}
      if (pal == 2){color_palette = c(color_palette, grad_rd[1:color_palette_list$rep[pal]])}
      if (pal == 3){color_palette = c(color_palette, grad_gr[1:color_palette_list$rep[pal]])}
    }
    
    png(paste0('./Paper/Latex_Table_Figure/Evolution_both.png'), width = 10, height = 6, units = 'in', res=200)
    d = comp_data %>% filter(country %in% c('Russia', 'U. K.')) %>% mutate(year = as.numeric(substr(year, 3, 4)))
    plot(ggplot(d,
                aes(fill=algo, y=index, x=year)) + 
           geom_line(aes(colour = algo, linetype = algo), size = 2) +
           scale_colour_manual("", values=color_palette) +
           scale_linetype_manual("", values=c('solid', 'twodash', 'dotted', 'dotdash')) +
           facet_wrap(country ~ factor, dir = 'v', scales = 'free', strip.position = 'top', ncol = 2) +
           scale_x_continuous(labels = paste0('\'', as.character(sort(unique(d$year)))), breaks = unique(d$year)) +
           theme(panel.background = element_rect(fill = "white", colour = "black"),
                 panel.grid = element_line(colour = "black", linetype = 'dashed', size = 0.4),
                 panel.grid.minor = element_blank(),
                 strip.text = element_text(size = 22),
                 axis.text.x = element_text(size = 20),
                 axis.text.y = element_text(size = 16),
                 legend.text=element_text(size=20),
                 legend.position = "bottom") +
           labs(x = '', y = ''))
    dev.off()
    
    
    # index couples for PCA only
    couple_list = list(c('Greece', 'Cyprus'),
                       c('Saudi Arabia', 'Russia'),
                       c('Argentina', 'Chile'),
                       c('Poland', 'Ukraine', 'Russia'),
                       c('India', 'China')
    )
    
    # https://en.wikipedia.org/wiki/List_of_economic_crises#21st_century
    counrty_crisis = data.frame(country = c('Greece', 'Cyprus', 'Saudi Arabia', 'Russia', 'Argentina', 'Chile',
                                            'Poland', 'Ukraine', 'India', 'China'),
                                start = c(10, 12, NA, 14, NA, NA, NA, 13, NA, NA),
                                end = c(17, 13, NA, 17, NA, NA, NA, 14, NA, NA), stringsAsFactors = F)
    
    for (coup in couple_list){
      png(paste0('./Paper/Latex_Table_Figure/Evolution_couple_', gsub(' ', '_', paste0(coup, collapse='_')), '_PCA_only.png'), width = 5*length(coup), height = 6, units = 'in', res=200) 
      d = comp_data %>% filter(country %in% coup) %>% filter(family == 'PCA') %>%
        mutate(year = as.numeric(substr(year, 3, 4))) %>% left_join(counrty_crisis, by = "country") %>%
        mutate(m = min(year),
               M = max(year)) %>%
        rowwise() %>%
        mutate(start = max(c(m, start)),
               end = min(c(M, end)))
      
      plot(ggplot(d,
                  aes(fill=algo, y=index, x=year)) + 
             geom_rect(
               aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
               fill = "red", alpha = 0.05) +
             geom_line(aes(colour = algo, linetype = algo), size = 2) +
             scale_colour_manual("", values=color_palette) +
             scale_linetype_manual("", values=c('solid', 'twodash', 'dotted', 'dotdash')) +
             facet_wrap(country ~ factor, dir = 'v', scales = 'free', strip.position = 'top', ncol = length(coup)) +
             scale_x_continuous(labels = paste0('\'', as.character(sort(unique(d$year)))), breaks = unique(d$year)) +
             theme(panel.background = element_rect(fill = "white", colour = "black"),
                   panel.grid = element_line(colour = "black", linetype = 'dashed', size = 0.4),
                   panel.grid.minor = element_blank(),
                   strip.text = element_text(size = 22),
                   axis.text.x = element_text(size = 20),
                   axis.text.y = element_text(size = 16),
                   legend.position = "none") +
             labs(x = '', y = ''))
      dev.off()
    }
    
    country_list = unique(comp_data$country)
    step = 30
    for (country_split in 1:4){
      country_range = ((country_split - 1) * step + 1):(country_split * step)
      
      png(paste0('./Paper/Latex_Table_Figure/Evolution_full_', country_split, '.png'), width = 15, height = 20, units = 'in', res=300) 
      plot(ggplot(comp_data %>% filter(country %in% country_list[country_range]),
                  aes(fill=algo, y=index, x=year)) + 
             geom_line(aes(colour = algo, linetype = algo), size = 1) +
             scale_colour_manual("", values=color_palette) +
             scale_linetype_manual("", values=c('solid', 'twodash', 'dotted', 'dotdash')) +
             facet_wrap(country ~ factor, dir = 'v', scales = 'free', strip.position = 'top', ncol = 6) +
             scale_x_continuous(labels = as.character(sort(unique(comp_data$year))), breaks = unique(comp_data$year)) +
             theme(panel.background = element_rect(fill = "white", colour = "black"),
                   panel.grid = element_line(colour = "black", linetype = 'dashed', size = 0.4),
                   panel.grid.minor = element_blank(),
                   strip.text = element_text(size = 15),
                   axis.text.x = element_text(size = 8),
                   legend.text=element_text(size=30),
                   legend.direction="horizontal",legend.position="top") +
             labs(x = '', y = ''))
      dev.off()
      
      png(paste0('./Paper/Latex_Table_Figure/Evolution_full_', country_split, '_PCA_only.png'), width = 15, height = 20, units = 'in', res=300) 
      plot(ggplot(comp_data %>% filter(country %in% country_list[country_range]) %>% filter(family == 'PCA'),
                  aes(fill=algo, y=index, x=year)) + 
             geom_line(aes(colour = algo, linetype = algo), size = 1) +
             scale_colour_manual("", values=color_palette) +
             scale_linetype_manual("", values=c('solid', 'twodash', 'dotted', 'dotdash')) +
             facet_wrap(country ~ factor, dir = 'v', scales = 'free', strip.position = 'top', ncol = 6) +
             scale_x_continuous(labels = as.character(sort(unique(comp_data$year))), breaks = unique(comp_data$year)) +
             theme(panel.background = element_rect(fill = "white", colour = "black"),
                   panel.grid = element_line(colour = "black", linetype = 'dashed', size = 0.4),
                   panel.grid.minor = element_blank(),
                   strip.text = element_text(size = 15),
                   axis.text.x = element_text(size = 8),
                   legend.position = "none") +
             labs(x = '', y = ''))
      dev.off()
    }
    
    # binary and continous for PCA only
    cont_list = comp_data %>%
      filter(algo == 'RobPCA') %>%
      mutate(year_factor = paste0(year, '_', factor),
             index = as.character(index)) %>%
      select(country, year_factor, index) %>%
      spread(year_factor, index)
    cont_list = data.frame(t(c('Country', sort(rep(unique(comp_data$year), each = 2)))), stringsAsFactors = F) %>%
      setNames(colnames(cont_list)) %>%
      bind_rows(
        data.frame(t(c('Country', rep(c('Ind 1', 'Ind 2'), uniqueN(comp_data$year)))), stringsAsFactors = F) %>%
          setNames(colnames(cont_list)),
        cont_list
      )
    write.table(cont_list, './Paper/Latex_Table_Figure/06_Cont_index_PCA_only.csv', sep = ";", col.names = F, row.names = F, append = F, dec = ".")
    
    bin_list = comp_data %>%
      filter(algo == 'RobPCA') %>%
      mutate(year_factor = paste0(year, '_', factor),
             index = as.character(ifelse(index > 0, 1, 0))) %>%
      select(country, year_factor, index) %>%
      spread(year_factor, index)
    bin_list = data.frame(t(c('Country', sort(rep(unique(comp_data$year), each = 2)))), stringsAsFactors = F) %>%
      setNames(colnames(bin_list)) %>%
      bind_rows(
        data.frame(t(c('Country', rep(c('Ind 1', 'Ind 2'), uniqueN(comp_data$year)))), stringsAsFactors = F) %>%
          setNames(colnames(bin_list)),
        bin_list
      )
    write.table(bin_list, './Paper/Latex_Table_Figure/07_Binary_index_PCA_only.csv', sep = ";", col.names = F, row.names = F, append = F, dec = ".")
  }
  
  # binary index validation plot
  df_type = 'Difference'
  recov_met = 'TENSOR_BF'
  res_outlier_stability = readRDS('./Checkpoints/res_outlier_stability.rds')
  res_thresh_sensitivity_best = readRDS('./Checkpoints/res_thresh_sensitivity_best.rds')
  {  
    res_outlier_stability_work = res_outlier_stability %>%
      mutate(Ind1 = ifelse(reg_type != 'index', NA, Ind1),
             Ind2 = ifelse(reg_type != 'index', NA, Ind2)) %>%
      mutate(Ind1 = as.numeric(Ind1),
             Ind2 = as.numeric(Ind2),
             maxAll = maxAll / tot_outlier,
             max1 = max1 / tot_outlier,
             max2 = max2 / tot_outlier,
             max3 = max3 / tot_outlier,
             max4 = max4 / tot_outlier) %>%
      mutate(outlier_score = (maxAll + 0.8 * max1 + 0.7 * max2 + 0.6 * max3 + 0.5 * max4) * 100) %>%
      mutate(outlier_score = ifelse(is.na(outlier_score), 0, outlier_score)) %>%
      mutate(mean_measure = ifelse(measure == 'abs_perc_err', mean_measure * 100, mean_measure),
             mean_measure_95 = ifelse(measure == 'abs_perc_err', mean_measure_95 * 100, mean_measure_95),
             mean_measure_99 = ifelse(measure == 'abs_perc_err', mean_measure_99 * 100, mean_measure_99),
             Model_Explain_var = Model_Explain_var * 100,
             Model_Explain_var_99 = Model_Explain_var_99 * 100,
             Model_Explain_var_95 = Model_Explain_var_95 * 100)
    z_lim = range(res_outlier_stability_work$outlier_score)
    z_lim_meas = res_outlier_stability_work %>%
      group_by(measure) %>%
      summarise(max = max(mean_measure)) %>%
      ungroup()
    z_lim_R2 = c(min(res_outlier_stability_work %>% select(starts_with('Model_Expl'))), 100)
    max_tot_outlier = max(res_outlier_stability_work$tot_outlier)
    upper_bin_max_tot_out = 150  # threshold to group all outlier above the value
    
    
    for(var_target in unique(res_outlier_stability_work$target_var)){
      
      set1 = res_outlier_stability %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        filter(target_var == var_target) %>%
        select(fit_method, algo_main) %>%
        unique()
      
      p_row = p_row_PCA = c()
      for (fit_met in unique(set1$fit_method)){
        for (algo_m in unique(set1$algo_main)){
          
          best_perf = res_thresh_sensitivity_best %>%
            filter(data == df_type) %>%
            filter(method == recov_met) %>%
            filter(target_var == var_target) %>%
            filter(fit_method == fit_met) %>%
            filter(algo_main == algo_m) %>%
            group_by(algo_main, regressor_type, PERF) %>%
            summarize(train_avg = round(mean(TRAIN_mean), 2),
                      train_std = round(sd(TRAIN_mean), 2),
                      test_avg = round(mean(TEST_mean), 2),
                      test_std = round(sd(TEST_mean), 2)) %>%
            ungroup() %>%
            mutate(train_lab = paste0(train_avg, ifelse(is.na(train_std), '', paste0('±', train_std))),
                   test_lab = paste0(test_avg, ifelse(is.na(test_std), '', paste0('±', test_std)))) %>%
            mutate(label = paste0(' -', ifelse(regressor_type == 'index', 'index (avg)', regressor_type), ' Train: ', train_lab, ' Test: ', test_lab)) %>%
            arrange(regressor_type)
          
          # row header
          image_lab = image_graph(res = 100, width = 570, height = 500, clip = F)
          plot(ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + xlim(0, 10) + ylim(2, 3.5) +
                 annotate("text", x = 0, y = 3.2, label = ifelse(fit_met == 'DFM_multivar', 'DFM', fit_met),
                          fontface = 2,
                          cex = 9,
                          hjust = 0, vjust = 0.7) +
                 annotate("text", x = 0, y = 3.17, label = paste0('\n\nModel: ', ifelse(algo_m == 'RandomForest', 'Random Forest', algo_m)),
                          fontface = 1,
                          cex = 8,
                          hjust = 0, vjust = 0.7) +
                 annotate("text", x = 0, y = 2.7, label = paste('\nPerformance:', toupper(unique(best_perf$PERF)),
                                                              '\nRegressor:',
                                                              paste0(c('', best_perf$label), collapse = '\n')),
                          cex = 7,
                          hjust = 0, vjust = 0.7) +
                 theme_bw() +
                 theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
                        plot.margin=unit(c(0,0.4,0,0.4),"cm"))
          )
          dev.off()
          
          # plot
          data = res_outlier_stability_work %>%
            filter(data == df_type) %>%
            filter(method == recov_met) %>%
            filter(fit_method == fit_met) %>%
            filter(target_var == var_target) %>%
            filter(algo_main == algo_m)
          
          p_hist = p_surf = c()
          meas = "abs_perc_err"
          
          # plot outlier distribution (histogram)
          data_t = data %>%
            filter(measure == meas)
          
          data_hist = data_t %>%
            filter(reg_type == 'index')
          
          z_mat = data_hist %>%
            select(Ind1, Ind2, outlier_score) %>%
            # mutate(outlier_score = ifelse(outlier_score <= 0, NA, outlier_score)) %>%
            spread(Ind2, outlier_score) %>%
            arrange(Ind1) %>%
            select(-Ind1) %>%
            as.matrix()
          
          z_color = data_hist %>%
            select(Ind1, Ind2, maxAll) %>%
            # mutate(maxAll = ifelse(maxAll <= 0, NA, maxAll)) %>%
            spread(Ind2, maxAll) %>%
            arrange(Ind1) %>%
            select(-Ind1) %>%
            as.matrix()
          
          floor_tot_outlier = data_hist %>%
            select(Ind1, Ind2, tot_outlier) %>%
            mutate(tot_outlier = ifelse(tot_outlier >= upper_bin_max_tot_out, upper_bin_max_tot_out * 2, tot_outlier)) %>%
            spread(Ind2, tot_outlier) %>%
            arrange(Ind1) %>%
            select(-Ind1) %>%
            as.matrix()
          
          hist = image_graph(res = 100, width = 800, height = 600, clip = F)
          par(mar=c(5,5,5,5))
          hist3D(z = z_mat, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                 border = "black",
                 xlab = "index 1", ylab = "index 2", zlab = "shared outlier (%)",
                 main = paste0('Outlier distribution for\n', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                 shade = 0,
                 phi = 20,  theta = -50,
                 ticktype = "detailed",
                 nticks = uniqueN(data_hist$Ind1),
                 space = 0.65,
                 alpha = 0.8,
                 bty = 'b2',
                 zlim = z_lim,
                 col = ramp.col (col = c("red", "blue3"), n = 100),
                 colvar = z_color,
                 cex.main = 2,
                 cex.lab = 1.5,
                 cex.axis = 1,
                 image = list(z = floor_tot_outlier, col = ramp.col (col = c("grey", "black"), n = 100, alpha = 0.5))
          )
          dev.off()
          hist = image_crop(hist, geometry_area(width = 500, height = 550, x_off = 120, y_off = 0))
          p_hist = c(p_hist, hist)
          
          
          # plot index stability according to meas (surface)
          z_mat = data_hist %>%
            select(Ind1, Ind2, mean_measure) %>%
            # mutate(mean_measure = ifelse(mean_measure <= 0, NA, mean_measure)) %>%
            spread(Ind2, mean_measure) %>%
            arrange(Ind1) %>%
            select(-Ind1) %>%
            as.matrix()
          
          z_mat_95 = data_hist %>%
            select(Ind1, Ind2, mean_measure_95) %>%
            # mutate(mean_measure_95 = ifelse(mean_measure_95 <= 0, NA, mean_measure_95)) %>%
            spread(Ind2, mean_measure_95) %>%
            arrange(Ind1) %>%
            select(-Ind1) %>%
            as.matrix()
          
          surf = image_graph(res = 100, width = 900, height = 600, clip = F)
          persp3D(z = z_mat, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                  facets = T, curtain = F,
                  xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error (%)', 'not found')),
                  main = paste0('Index stability for\n', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                  col = 'darkslategrey', border = 'black',
                  shade = 0,
                  phi = 20,  theta = -50,
                  ticktype = "detailed",
                  nticks = uniqueN(data_hist$Ind1),
                  zlim = c(0, unlist(z_lim_meas %>% filter(measure == meas) %>% select(max))),
                  cex.main = 2,
                  cex.lab = 1.5,
                  cex.axis = 1,
                  bty = 'b2'
          )
          # text3D(max(data_hist$Ind1), min(data_hist$Ind2), mean(z_mat, na.rm = T), paste0(' index ', round(mean(z_mat, na.rm = T), 2), '-', round(mean(z_mat_95, na.rm = T), 2)), colvar = NULL, add = T)
          # text3D(max(data_hist$Ind1), min(data_hist$Ind2), mean(z_mat, na.rm = T), ' index', colvar = NULL, add = T)
          persp3D(z = z_mat_95, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                  facets = T, curtain = F,
                  xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error (%)', 'not found')),
                  main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                  col = 'darkgrey', border = 'black',
                  shade = 0,
                  phi = 20,  theta = -50,
                  ticktype = "detailed",
                  nticks = uniqueN(data_hist$Ind1),
                  zlim = c(0, unlist(z_lim_meas %>% filter(measure == meas) %>% select(max))),
                  bty = 'b2',
                  add = T
          )
          # text3D(max(data_hist$Ind1), min(data_hist$Ind2), mean(z_mat_95, na.rm = T), ' 95% of index', colvar = NULL, add = T)
          label_list = data.frame(val = c(mean(z_mat, na.rm = T), mean(z_mat_95, na.rm = T)), label = c(' index', ' 95% of index'), stringsAsFactors = F)
          
          data_surf = data_t %>%
            filter(reg_type != 'index') %>%
            mutate(color = c('blue', 'coral', 'chartreuse'),
                   label = paste0(reg_type, ' ', round(mean_measure, 2), '-', round(mean_measure_95, 2)))
          for (rr in unique(data_surf$reg_type)){
            r_col = data_surf %>% filter(reg_type == rr) %>% select(color) %>% unlist() %>% setNames(NULL)
            r_val = data_surf %>% filter(reg_type == rr) %>% select(mean_measure) %>% unlist() %>% setNames(NULL)
            r_lab = data_surf %>% filter(reg_type == rr) %>% select(label) %>% unlist() %>% setNames(NULL)
            persp3D(z = z_mat_95, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                    facets = T, curtain = F,
                    xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error (%)', 'not found')),
                    main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                    col = 'darkgrey', border = 'black',
                    shade = 0,
                    phi = 20,  theta = -50,
                    ticktype = "detailed",
                    nticks = uniqueN(data_hist$Ind1),
                    zlim = c(0, unlist(z_lim_meas %>% filter(measure == meas) %>% select(max))),
                    bty = 'b2',
                    add = T,
                    image = list(side = r_val, facets = T, col = r_col, border = 'black', z = z_mat)
            )
            # text3D(max(data_hist$Ind1), min(data_hist$Ind2), r_val, r_lab, colvar = NULL, add = T)
            label_list = label_list %>% bind_rows(
              data.frame(val = r_val, label = r_lab, stringsAsFactors = F)
            )
          }
          label_list = space_label(label_list, 0.07 * unlist(z_lim_meas %>% filter(measure == meas) %>% select(max)))
          text3D(rep(max(data_hist$Ind1), nrow(label_list)), rep(min(data_hist$Ind2), nrow(label_list)), label_list$val, label_list$label, colvar = NULL, add = T, cex = 1.2)
          dev.off()
          surf = image_crop(surf, geometry_area(width = 650, height = 550, x_off = 170, y_off = 0))
          p_surf = c(p_surf, surf)
          
          
          
          
          # create legend
          bar_legend = image_graph(res = 100, width = 300, height = image_info(p_hist[[1]])$height / 2, clip = F)
          grid.draw(cowplot::get_legend(
            ggplot(data = data.frame(val = unique(res_outlier_stability_work$maxAll)) %>% mutate(x = 1:n()),
                   aes(x = x, y = val, fill = val)) +
              geom_point() +
              scale_fill_gradientn(colours=ramp.col (col = c("red", "blue3"), n = 100),
                                   breaks=c(0,1),labels=c("0% of combinations", "100% of combinations"),
                                   limits=c(0,1),
                                   name = 'Bar color:\nShare of outliers\n') +
              theme(legend.title=element_text(size=18), 
                    legend.text=element_text(size=16))
          ))
          dev.off()
          
          floor_legend= image_graph(res = 100, width = 300, height = image_info(p_hist[[1]])$height / 2, clip = F)
          grid.draw(cowplot::get_legend(
            ggplot(data = res_outlier_stability_work %>%
                     select(tot_outlier) %>%
                     mutate(tot_outlier = ifelse(tot_outlier >= upper_bin_max_tot_out, upper_bin_max_tot_out * 2, tot_outlier)) %>%
                     unique() %>%
                     mutate(x = 1:n()),
                   aes(x = x, y = tot_outlier, fill = tot_outlier)) +
              geom_point() +
              scale_fill_gradientn(colours=ramp.col (col = c("grey", "black"), n = 100, alpha = 0.5),
                                   breaks=c(0,1),labels=c(min(res_outlier_stability_work$tot_outlier), "300+"),
                                   limits=c(0,1),
                                   name = 'Floor color:\nNumber of total outliers\n') +
              theme(legend.title=element_text(size=18), 
                    legend.text=element_text(size=16))
          ))
          dev.off()
          
          regr_legend = image_graph(res = 100, width = 350, height = image_info(p_surf[[1]])$height, clip = F)
          grid.draw(cowplot::get_legend(
            ggplot(data = data.frame(reg = c('index', '95% of index', 'original (full-95%)', 'rank index (full-95%)', 'raw index (full-95%)')), aes(reg, fill = reg)) + 
              geom_bar() +
              scale_fill_manual(name = 'Regressor:', values = c('darkslategrey', 'darkgrey', 'blue', 'coral', 'chartreuse')) +
              theme(legend.title=element_text(size=23), 
                    legend.text=element_text(size=20))
          ))
          dev.off()
          
          # assemble row plot
          p_row = c(p_row, image_append(c(image_lab, p_hist, image_append(c(bar_legend, floor_legend), stack = T), p_surf, regr_legend)))
          if (fit_met == 'RobPCA'){
            p_row_PCA = c(p_row_PCA, image_append(c(image_lab, p_hist, image_append(c(bar_legend, floor_legend), stack = T), p_surf, regr_legend)))
          }
          
        } # algo_m
      } # fit_met
      
      eval(parse(text=paste0('final_plot = image_append(c(', paste0('p_row[[', 1:length(p_row), ']]', collapse = ','), '), stack = T)')))
      eval(parse(text=paste0('final_plot_PCA = image_append(c(', paste0('p_row_PCA[[', 1:length(p_row_PCA), ']]', collapse = ','), '), stack = T)')))
      
      # final_plot = image_crop(final_plot, geometry_area(width = 2500, height = 600, x_off = 0, y_off = 0))
      
      v_tag = substr(gsub('var_add_WB_WDI_', '', var_target), 1, 20)
      png(paste0('./Paper/Latex_Table_Figure/Index_Val_', v_tag, '.png'), width = 9, height = 2.1 * length(p_row), units = 'in', res=300)
      par(mar=c(0,0,0,0))
      par(oma=c(0,0,0,0))
      plot(final_plot)
      dev.off()
      
      png(paste0('./Paper/Latex_Table_Figure/Index_Val_', v_tag, '_PCA_only.png'), width = 9, height = 2.1 * length(p_row_PCA), units = 'in', res=300)
      par(mar=c(0,0,0,0))
      par(oma=c(0,0,0,0))
      plot(final_plot_PCA)
      dev.off()
    } # var_target
    
  }
}

