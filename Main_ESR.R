
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
library(EnvStats)
library(parallelMap)
library(parallel)
library(party)
library(randomForest)
library(gbm)
# library(crs)  # spline
# library(h2o)
# h2o.init(nthreads = -1, max_mem_size = '4g', ip = "127.0.0.1", port = 54321)
library(mlr)
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
force_remove_variable = c('rain', 'temp', 'bcrisis', 'ehr1', 'ehr2', 'ehr3')  # totally full of missing for 2017 to 2019
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
  variable_max_perc = 70
  
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
  write.table(df %>%
                gather('variable', 'val', -c(year, country)) %>%
                group_by(country) %>%
                summarize(Tot_non_NA = sum(!is.na(val)),
                          Tot_NA = sum(is.na(val)), .groups = 'drop') %>%
                mutate(NA_perc = round(Tot_NA / (Tot_NA + Tot_non_NA) * 100, 1)) %>%
                arrange(NA_perc), "./Stats/01_final_df_missing_by_country.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  write.table(df %>%
                gather('variable', 'val', -c(year, country)) %>%
                group_by(year) %>%
                summarize(Tot_non_NA = sum(!is.na(val)),
                          Tot_NA = sum(is.na(val)), .groups = 'drop') %>%
                mutate(NA_perc = round(Tot_NA / (Tot_NA + Tot_non_NA) * 100, 1)) %>%
                arrange(year), "./Stats/01_final_df_missing_by_year.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")

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
# restricted_df_var = c('growth', 'interest1', 'labor1', 'penn19', 'penn32', 'penn33', 'penn42',
#                       'religdiv1', 'religdiv2', 'chrs_pct', 'musl_pct',
#                       df %>% select(starts_with("var")) %>% colnames(),
#                       df %>% select(starts_with("health")) %>% colnames(),
#                       df %>% select(starts_with("hef")) %>% colnames(),
#                       df %>% select(starts_with("pop")) %>% colnames())
# restricted_df_var = setdiff(restricted_df_var, c('health3', 'pop1', 'pop2', 'pop3', 'var8'))
restricted_df_var = setdiff(df %>% select(starts_with("var")) %>% colnames(), 'var18')

### recover missing and create final dataset
{
  ### recover missing data
  dummy_var = c()#c('bcrisis')    # apply median
  recov_method_set = c('SOFT_IMPUTE', 'TENSOR_BF', 'FULL_AVERAGE')#, 'COUNTRY_AVERAGE')#,'NIPALS')
  df_set = c('Restricted') #, 'Original')  # "Restricted" is without multicollinearity
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
}
df_final = readRDS('./Checkpoints/df_final.rds')



### evaluate Kaiser–Meyer–Olkin statistic of sampling adequacy and
#   Im, Pesaran and Shin  panel data unit root statistic to check the stationarity of the index values
{
  # KMO
  # >0.9 marvelous, [0.8,0.9) meritorious, [0.7,0.8) middling, [0.6,0.7) mediocre, [0.5, 0.6) miserable, <0.5 unacceptable
  
  # IPS H1: stationarity, tested for both "individual intercepts" and "individual intercepts and trends" option for Augmented-Dickey-Fuller
  
  KMO_IPS_test = c()
  for (df_type in setdiff(unique(df_final$data), c('ExtendedOriginal', 'ExtendedDifference'))){
    for (recov_met in unique(df_final$method)){

      df_final_spread = df_final %>%    # df_final in wide format
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        spread(variable, val) %>%
        select(-method, -data)
      pdata = plm::pdata.frame(df_final_spread, index=c("country","year"), drop.index=TRUE, row.names=TRUE)
      
      KMO_IPS_test = KMO_IPS_test %>%
        bind_rows(data.frame(data = df_type, method = recov_met,
                             KMO = REdaS::KMOS(df_final_spread %>% select(-country, -year))$KMO,
                             KMO_val = ">0.9 marvelous, [0.8,0.9) meritorious, [0.7,0.8) middling, [0.6,0.7) mediocre, [0.5, 0.6) miserable, <0.5 unacceptable",
                             IPS_intercept_pval = suppressWarnings(plm::purtest(pdata, test = "ips", exo = "intercept", lags = "AIC", pmax = 5)$statistic$p.value),
                             IPS_trend_pval = suppressWarnings(plm::purtest(pdata, test = "ips", exo = "trend", lags = "AIC", pmax = 5)$statistic$p.value),
                             IPS_H1 = 'stationarity', stringsAsFactors = F))
    }
  }
  write.table(KMO_IPS_test, "./Stats/02_sampling_and_stationarity_test.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
}




### PCA and robust PCA
{
  cv_eval = F # to Cross-Validate PCA for suggested number of PC
  k_fold = 3
  pca_met_set = c('PCA', 'RobPCA', 'RobSparPCA')
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
  leading_var = "var1"
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
  leading_var = "var1"
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
  
  
  
  
  recov_met_set = c('SOFT_IMPUTE', 'TENSOR_BF')
  
  ### evaluate PCA final index based on PC (=scores)
  PC_to_keep = 2
  index_1_thresh = 0  # threshold to split index 1 (first PC scores)
  index_2_thresh = 0  # threshold to split index 2 (second PC scores) - if PC_to_keep == 2
  load_thresh = 0  # loadings with magnitude below threshold are set to 0 and scores are evaluated always as X * loadings
  leading_var = "var1"
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
}



### Dynamic Factor Model
{
  df_type_to_keep = c("Restricted", "RestrictedDifference")
  recov_met_to_keep = c("TENSOR_BF", "SOFT_IMPUTE")
  
  n_factor_tot = 1  # number of factors to evaluate (all incremental sets up to n_factor_tot are evaluated)
  VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
  kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
  univ_reload = T  # reload previous evaluation (if available) for evaluate_DFM_univar
  multiv_reload_Kalm = T  # reload previous evaluation (if available) for evaluate_DFM_multivar
  continue_kalman_loop = F  # if True and saved kalman did not converge, continue looping from the last point
  # (max_kalman_loop should be grater than stored while_count)
  # This will set multiv_reload_Kalm = F (but var_to_remove_all and while_count will be reloded)
  univ_save = F  # save intermediate evaluation
  total_save = T  # save RDS for univariate + multivariate for each pair df_type + recov_met
  {
    res_DFM_factors = res_DFM_loadings = res_DFM_stats = res_DFM_MAPE = c()
    sink(paste0('./Log/DFM_fitting_log_', format(Sys.time(), "%Y-%m-%dT%H%M"), '.txt'), append = F, split = T, type = 'output')
    cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
    for (df_type in df_type_to_keep){ #unique(df_final$data)){
      # if (df_type == "Restricted"){rr = recov_met_to_keep} else {rr = "SOFT_IMPUTE"}
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
                                            kalm_Q_hat_mode = kalm_Q_hat_mode, max_kalman_loop = 200,
                                            multiv_reload_Kalm = multiv_reload_Kalm, res_DFM_list_reload = res_DFM_list_reload,
                                            multivar_mode = 'DFM_multivar', continue_kalman_loop = continue_kalman_loop)
          res_DFM_stats = DFM_multi$res_DFM_stats
          res_DFM_list = DFM_multi$res_DFM_list
          res_DFM_factors = DFM_multi$res_DFM_factors
          res_DFM_loadings = DFM_multi$res_DFM_loadings
          res_DFM_MAPE = DFM_multi$res_DFM_MAPE
          
          # stack univariate DFM (factors are standardized)
          cat('\n\n\n  --------------------------------- Stacking Univariate DFM with', n_factor, ifelse(n_factor == 1, 'factor', 'factors') ,'models   ---------------------------------\n\n')
          
          DFM_multi = evaluate_DFM_multivar(res_DFM_factors, res_DFM_loadings, res_DFM_stats, res_DFM_list, res_DFM_MAPE,
                                            df_SP_orig, df_type, recov_met, n_factor,
                                            VAR_alpha = VAR_alpha,
                                            kalm_Q_hat_mode = kalm_Q_hat_mode,
                                            res_DFM_list_reload = res_DFM_list_reload,
                                            multivar_mode = 'DFM_multivar_stacked')
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
  
  
  
  
  
  dfm_met_to_plot = c("DFM_multivar", "DFM_multivar_stacked")
  
  
  ### evaluate DFM final index based on Factors
  # only 1 and 2 factors plot are showed
  index_1_thresh = 0  # threshold to split index 1 (first factor)
  index_2_thresh = 0  # threshold to split index 2 (second factor)
  load_thresh = 0  # loadings with magnitude below threshold are set to 0
  leading_var = "var1"
  leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
  VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
  kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
  {
    res_index_DFM = list()
    for (df_type in setdiff(unique(res_DFM_loadings$data), c('ExtendedDifference', 'ExtendedOriginal'))){
      for (recov_met in res_DFM_loadings %>% filter(data == df_type) %>% pull(method) %>% unique()){
        
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
      for (recov_met in res_DFM_loadings %>% filter(data == df_type) %>% pull(method) %>% unique()){
        
        RDS_lab = paste0(gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '_', df_type, '_', recov_met, '.rds')
        res_DFM_list = readRDS(paste0('./Checkpoints/DFM/', RDS_lab))
        
        for (dfm_met in setdiff(dfm_met_to_plot, c('DFM_univar', 'DFM_multivar_stacked'))){
          
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
}

### plot loading ranking for selected method, df_type, recov_met   - used only for a meeting
df_type = "Restricted"
recov_met = "SOFT_IMPUTE"
DFM_met = "DFM_multivar"
n_factor = 1
{
  # loading_table = res_DFM_loadings %>%
  #   filter(data == df_type) %>%
  #   filter(method == recov_met) %>%
  #   filter(DFM == DFM_met) %>%
  #   filter(Factor == n_factor) %>%
  #   group_by(variable) %>%
  #   summarize(loading = mean(loading), .groups = "drop") %>%
  #   arrange(desc(abs(loading))) %>%
  #   left_join(read.csv("Variable_description_Restricted_Dataset.csv", sep=";", stringsAsFactors=FALSE), by = c("variable" = "var")) %>%
  #   setNames(c('Variable', 'Loading', 'Description')) %>%
  #   mutate(xaxis = 1:n()*2)
  # loading_table$Variable = factor(loading_table$Variable, levels = rev(loading_table$Variable))
  # 
  # p = ggplot(data=loading_table, aes(x=Loading, y=Variable)) +
  #   geom_bar(stat="identity") +
  #   theme(plot.margin = unit(c(1,40,1,1), "lines"),
  #         axis.text.x = element_text(size = 15),
  #         axis.text.y = element_text(size = 15),
  #         axis.title = element_text(size = 17, face = 'bold'))
  # for (i in 1:nrow(loading_table))  {
  #   p <- p + annotation_custom(
  #     grob = textGrob(label = rev(loading_table$Description)[i], hjust = 0, gp = gpar(cex = 1.5)),
  #     ymin = i,      # Vertical position of the textGrob
  #     ymax = i,
  #     xmin = max(loading_table$Loading)*1.3,         # Note: The grobs are positioned outside the plot area
  #     xmax = max(loading_table$Loading)*1.3)
  # }
  # gt <- ggplot_gtable(ggplot_build(p))
  # gt$layout$clip[gt$layout$name == "panel"] <- "off"
  # png(paste0('pp.png'), width = 20, height = 10, units = 'in', res=300)
  # grid.draw(gt)
  # dev.off()
}



### plot index evaluation and comparison
{
  df_type_to_keep = c("Restricted" , "RestrictedDifference")
  recov_met_to_keep = c("TENSOR_BF", "SOFT_IMPUTE")
  PCA_to_keep = c('RobPCA')
  DFM_to_keep = c('DFM_multivar', 'DFM_multivar_stacked')
  
  
  ### plot comparison of index evolution over time for all methods
  res_PCA_loadings = readRDS('./Checkpoints/res_PCA_loadings.rds')
  res_index_PCA = readRDS('./Checkpoints/res_index_PCA.rds')
  VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
  kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
  res_DFM_factors = readRDS(paste0('./Checkpoints/res_DFM_factors_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  res_DFM_loadings = readRDS(paste0('./Checkpoints/res_DFM_loadings_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  res_DFM_best_model = readRDS(paste0('./Checkpoints/res_DFM_best_model_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  res_index_DFM = readRDS('./Checkpoints/res_index_DFM.rds')
  algo_to_plot = c('RobPCA', 'DFM_multivar', 'DFM_multivar_stacked')
  index_to_keep = 1  #  index dimension to keep 
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
      for (recov_met in intersect(res_DFM_loadings %>% filter(data == df_type) %>% pull(method) %>% unique(), recov_met_to_keep)){
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
    
    # select number of index
    res_ALL_scores = res_ALL_scores %>%
      filter(factor <= index_to_keep)
    
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
  algo_to_plot = c('RobPCA', 'DFM_multivar', 'DFM_multivar_stacked')
  algo_to_plot_lab = c('RobPCA', 'DFM', 'DFM stacked')   # same length of algo_to_plot
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
          mutate(index = index + min(index)) %>%  # in order to make all values positive and avoid 0 mean
          summarise(Variation_coeff = sd(index) / abs(mean(index)), .groups = 'drop') %>%
          mutate(quantile = quantile(Variation_coeff, 0.90)) %>%
          ungroup() %>%
          left_join(res_ALL_scores %>% select(algo, family) %>% unique(), by = "algo") %>%
          left_join(data.frame(old = algo_to_plot, new = algo_to_plot_lab, stringsAsFactors = F), by = c('algo' = 'old')) %>%
          select(-algo) %>%
          rename(algo = new) %>%
          filter(Variation_coeff <= quantile)
        
        png(paste0('./Results/2_', df_type, '_', recov_met, '_factors_variation_coeff.png'), width = 10, height = 10, units = 'in', res=300) 
        plot(ggplot(data,
                    aes(x=algo, y=Variation_coeff, color=family)) +
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
}



### compare index with disease incidence indicators
p_val_tol = 0.01 # p-val tolerance for correlation test
algo_to_check = c('RobPCA', 'DFM_multivar', 'DFM_multivar_stacked')
res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds') %>% mutate(year = as.character(year))
{
  coutry_mapping_disease = data.frame(country_original = c("Bahamas", "Bolivia (Plurinational State of)", "Congo", "Congo", "CÃ´te dâ€™Ivoire", "Czechia",
                                                           "Egypt", "Gambia", "Iran (Islamic Republic of)", "Democratic People's Republic of Korea",
                                                           "Kyrgyzstan", "Lao People's Democratic Republic", "Republic of Moldova",
                                                           "The former Yugoslav Republic of Macedonia", "Slovakia", "Saint Lucia",
                                                           "Saint Vincent and the Grenadines", "United Republic of Tanzania",
                                                           "United Kingdom of Great Britain and Northern Ireland", "United States of America",
                                                           "Venezuela (Bolivarian Republic of)", "Viet Nam", "Yemen"),
                                      country = c("Bahamas, The", "Bolivia", "Congo, Dem. Rep.", "Congo, Rep.", "Cote d'Ivoire", "Czech Republic",
                                                  "Egypt, Arab Rep.", "Gambia, The", "Iran, Islamic Rep.", "Korea, Rep.",
                                                  "Kyrgyz Republic", "Lao PDR", "Moldova", "North Macedonia", "Slovak Republic", "St. Lucia",
                                                  "St. Vincent and the Grenadines", "Tanzania", "United Kingdom", "United States",
                                                  "Venezuela, RB", "Vietnam", "Yemen, Rep."), stringsAsFactors = F)
  
  data_disease <- read.csv("C:/Users/Alessandro Bitetto/Downloads/UniPV/Charilaos_ESR/Data/data_Tropical.csv", stringsAsFactors=FALSE) %>%
    select(Location, Period, FactValueNumeric) %>%
    rename(country_original = Location,
           year = Period,
           value = FactValueNumeric) %>%
    mutate(year = as.character(year),
           disease = "Tropical") %>%
    left_join(coutry_mapping_disease, by = "country_original") %>%
    mutate(country = ifelse(is.na(country), country_original, country)) %>%
    select(-country_original) %>%
    filter(year %in% unique(res_ALL_scores$year)) %>%
    
    bind_rows(
      read.csv("C:/Users/Alessandro Bitetto/Downloads/UniPV/Charilaos_ESR/Data/data_tubercolosis.csv", stringsAsFactors=FALSE) %>%
        select(Location, Period, FactValueNumeric) %>%
        rename(country_original = Location,
               year = Period,
               value = FactValueNumeric) %>%
        mutate(year = as.character(year),
               disease = "TBC") %>%
        left_join(coutry_mapping_disease, by = "country_original") %>%
        mutate(country = ifelse(is.na(country), country_original, country)) %>%
        select(-country_original) %>%
        filter(year %in% unique(res_ALL_scores$year))
    ) %>%
    
    bind_rows(
      read.csv("C:/Users/Alessandro Bitetto/Downloads/UniPV/Charilaos_ESR/Data/data_Malaria.csv", stringsAsFactors=FALSE) %>%
        select(Location, Period, FactValueNumeric) %>%
        rename(country_original = Location,
               year = Period,
               value = FactValueNumeric) %>%
        mutate(year = as.character(year),
               disease = "Malaria") %>%
        left_join(coutry_mapping_disease, by = "country_original") %>%
        mutate(country = ifelse(is.na(country), country_original, country)) %>%
        select(-country_original) %>%
        filter(year %in% unique(res_ALL_scores$year))
    ) %>%
    
    bind_rows(
      read.csv("C:/Users/Alessandro Bitetto/Downloads/UniPV/Charilaos_ESR/Data/data_HIV.csv", sep=";", stringsAsFactors=FALSE) %>%
        select(country, year, Incidence.of.HIV..all..per.1.000.uninfected.population...SH.HIV.INCD.TL.P3.) %>%
        rename(value = Incidence.of.HIV..all..per.1.000.uninfected.population...SH.HIV.INCD.TL.P3.) %>%
        mutate(year = as.character(year),
               disease = "HIV") %>%
        filter(year %in% unique(res_ALL_scores$year))
    )  
  
  res_correlation = c()
  for (df_type in unique(res_ALL_scores$data)){
    for (recov_met in unique(res_ALL_scores$method)){
      for (algo_type in algo_to_check){
        
        index = res_ALL_scores %>%
          filter(data == df_type) %>%
          filter(method == recov_met) %>%
          filter(algo == algo_type) %>%
          select(country, year, index)
        
        if (nrow(index) > 0){
          
          # test different reference indexes
          for (ref_ind in data_disease$disease %>% unique()){
            
            index_ref_all = index %>%
              left_join(data_disease %>%
                          filter(disease == ref_ind) %>%
                          rename(reference = value) %>%
                          select(-disease), by = c("year", "country")) %>%
              filter(!is.na(reference))
            
            # check minimum number of observations
            check_min_obs = index_ref_all %>%
              group_by(country) %>%
              summarise(CC = n())
            if (min(check_min_obs$CC) <= 6){
              cat(paste0(c(df_type, recov_met, algo_type, ref_ind), collapse = " - "), "    - have less than 7 observations to evaluate correlation", end = "\n")
            }
            if (sum(is.na(index_ref_all)) > 0){
              cat("######   ", paste0(c(df_type, recov_met, algo_type, ref_ind), collapse = " - "), "    - Missing values detected", end = "\n")
            }
            
            # evaluate correlation for each country over all years
            tot_country = 0
            for (countr in index_ref_all$country %>% unique()){
              
              index_by_country = index_ref_all %>%
                filter(country == countr)
              
              if (diff(range(index_by_country$reference)) != 0 & diff(range(index_by_country$index)) != 0){
                
                tot_country = tot_country + 1
                
                cor_kendall = suppressWarnings(cor.test(index_by_country$index, index_by_country$reference,  method="kendall", exact = T)) # null hypotesis is 0 correlation
                cor_spearman = suppressWarnings(cor.test(index_by_country$index, index_by_country$reference,  method="spearman", exact = T)) # null hypotesis is 0 correlation
                cor_somers = rcorr.cens(index_by_country$index, index_by_country$reference, outx = TRUE)
                kruskal = kruskal.test(index ~ reference, data = index_by_country)  # null hypotesis is that the distributions are the same
                
                res_correlation = res_correlation %>%
                  bind_rows(data.frame(method = recov_met,
                                       data = df_type,
                                       algo = algo_type,
                                       reference_index = ref_ind,
                                       country = countr,
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
                                       kruskalWallis_warn = ifelse(kruskal$p.value <= p_val_tol, '*', ''), stringsAsFactors = F)
                  )
              } # check if both indixes are not constant
            } # countr
            
            # plot correlation
            plot_corr = res_correlation %>%
              filter(method == recov_met) %>%
              filter(data == df_type) %>%
              filter(algo == algo_type) %>%
              filter(reference_index == ref_ind) %>%
              rename(kruskalWallis_corr = kruskalWallis_pVal_high_means_same) %>%
              select(-method, -data, -algo, -reference_index, -ends_with('pVal'), -somers_confidence) %>%
              gather('corr', 'val', -c(country, ends_with('warn'))) %>%
              mutate(warn = '')  %>%
              mutate(warn = ifelse(kendall_warn == '*' & corr == 'kendall_corr', '*', warn)) %>%
              mutate(warn = ifelse(spearman_warn == '*' & corr == 'spearman_corr', '*', warn)) %>%
              mutate(warn = ifelse(somers_warn == '*' & corr == 'somers_corr', '*', warn)) %>%
              mutate(warn = ifelse(kruskalWallis_warn == '*' & corr == 'kruskalWallis_corr', '*', warn)) %>%
              mutate(corr = gsub('kruskalWallis_corr', 'Kruskal-Wallis p-Val \n(high means same distribution)', corr)) %>%
              mutate(corr = gsub('_corr', '', corr)) %>%
              
              mutate(corr = capitalize(corr)) %>%
              select(-ends_with('_warn')) %>%
              mutate(warn = as.factor(warn),
                     corr = factor(corr, levels = c('Kendall', 'Spearman', 'Somers', 'Kruskal-Wallis p-Val \n(high means same distribution)', '% of matching countries')))

            png(paste0('./Results/4_correlation_power_', df_type, '_', recov_met, '_', algo_type, '_', ref_ind, '.png'), width = 10, height = 10, units = 'in', res=300)
            plot(
              ggplot(data = plot_corr, aes(x = val)) +
                   geom_density(aes(y=..density..), size = 1, fill = "grey") +
                   scale_x_continuous(labels =  scales::percent_format(accuracy = 1)) +
                   facet_wrap(~corr, ncol = 2, dir = 'h', scales = 'free') +
                   ggtitle(paste0('Distribution of correlations with ', ref_ind),
                           subtitle = paste0("Selected countries: ", tot_country, "/", res_ALL_scores$country %>% uniqueN(), '\n')) +
                   theme(axis.text.x = element_text(size = 16),
                         axis.title=element_blank(),
                         axis.text.y = element_blank(),
                         text = element_text(size=20),
                         title = element_text(size=24),
                         plot.subtitle = element_text(size=20),
                         panel.background = element_rect(fill = "white", colour = "black"),
                         panel.grid.major = element_line(colour = "grey", linetype = 'dashed', size = 1),
                         strip.text = element_text(size = 19),
                         panel.spacing = unit(2, "lines"))
              )
            dev.off()
            
          } # ref_ind
        } # nrow(index) > 0
      } # algo_type
    } # recov_met
  } # df_type
  
  write.table(res_correlation, './Results/4_correlation_power_summary.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  rm(coutry_mapping_disease, data_disease, check_min_obs, plot_corr, cor_kendall, cor_spearman, index, index_by_country, index_ref_all, kruskal, res_correlation)
}



### reload and select additional variables to be used as target for regression
# load new variables and standardize (scale in [0,1] or [-1,1]) and recover missing
df_type_reference = 'Restricted'
recov_met_reference = 'SOFT_IMPUTE'
df_additional_raw_reload = F   # if FALSE read raw datasets and match country names
{
  dataset_ref = data.frame(path = c(
    "./Data/Financial_access_data_2004-2019_Bitteto.dta",
    "./Data/IFS_data_2006_2019_Bitteto.dta",
    "./Data/cog_data.dta"
    # "./Data/cpia_data.dta",
    # "./Data/wef_data.dta",
    # "./Data/wgi_data.dta",
    # "./Data/FD_data.dta"
  ),
  label = c('FIN_ACC', 'IFS', 'COG'), stringsAsFactors = F)
  
  country_mapping = data.frame(country_name = c("Afghanistan", "Armenia", "Belarus", "Central African Republic", "Azerbaijan", "China", "Comoros",
                                                "Congo, Dem. Rep.", "Congo, Rep.", "Cote d'Ivoire", "Egypt, Arab Rep.", "Equatorial Guinea",
                                                "Eswatini", "Ethiopia", "Fiji", "Hong Kong SAR, China", "Iran, Islamic Rep.",
                                                "Korea, Rep.", "Lao PDR", "Mauritania", "Mozambique", "Netherlands", "North Macedonia", 
                                                "Sao Tome and Principe", "Serbia", "Syrian Arab Republic", "Tajikistan", "Uzbekistan",
                                                "Venezuela, RB", "Yemen, Rep.",
                                                "Croatia", "Czech Republic", "Dominican Republic", "Estonia", "Kazakhstan", "Kosovo, Republic of",
                                                "Kyrgyz Republic", "Lesotho", "Madagascar", "Moldova", "Macedonia, FYR", "Poland", "Slovak Republic",
                                                "Slovenia", "Tanzania"),    # df_final names
                               country = c("Afghanistan, Islamic Rep. of", "Armenia, Rep. of", "Belarus, Rep. of", "Central African Rep.",
                                           "Azerbaijan, Rep. of", "China, P.R.: Mainland", "Comoros, Union of the", "Congo, Dem. Rep. of the",
                                           "Congo, Rep. of", "Côte d'Ivoire", "Egypt, Arab Rep. of", "Equatorial Guinea, Rep. of", "Eswatini, Kingdom of",
                                           "Ethiopia, The Federal Dem. Rep. of", "Fiji, Rep. of", "China, P.R.: Hong Kong", "Iran, Islamic Rep. of",
                                           "Korea, Rep. of", "Lao People's Dem. Rep.", "Mauritania, Islamic Rep. of", "Mozambique, Rep. of",
                                           "Netherlands, The", "North Macedonia, Republic of", "São Tomé and Príncipe, Dem. Rep. of",
                                           "Serbia, Rep. of", "Syrian Arab Rep.", "Tajikistan, Rep. of", "Uzbekistan, Rep. of",
                                           "Venezuela, Rep. Bolivariana de", "Yemen, Rep. of",
                                           "Croatia, Rep. of", "Czech Rep.", "Dominican Rep.", "Estonia, Rep. of", "Kazakhstan, Rep. of","Kosovo, Rep. of",
                                           "Kyrgyz Rep.", "Lesotho, Kingdom of", "Madagascar, Rep. of", "Moldova, Rep. of", "North Macedonia, Republic of",
                                           "Poland, Rep. of", "Slovak Rep.", "Slovenia, Rep. of", "Tanzania, United Rep. of"), stringsAsFactors = F)
  country_mapping_by_iso3 = read_dta("./Data/Financial_access_data_2004-2019_Bitteto.dta") %>%
    select(country, iso3) %>%
    unique() %>%
    left_join(country_mapping, by = "country") %>%
    mutate(country_name = ifelse(is.na(country_name), country, country_name)) %>%
    select(iso3, country_name) %>%
    bind_rows(data.frame(iso3 = c('XKX'),
                         country_name = c('Kosovo, Republic of')))
  
  df_ref = readRDS('./Checkpoints/df_final.rds') %>%
    filter(data == df_type_reference) %>%
    filter(method == recov_met_reference) %>%
    spread(variable, val) %>%
    select(-method, -data) %>%
    select(country, year) %>%
    rename(country_name = country,
           year_name = year) %>%
    mutate(year_name = as.character(year_name))

  if (df_additional_raw_reload){
    dataset_list = list()
    for (i in 1:nrow(dataset_ref)){
      lab = dataset_ref$label[i]
      tt = read_dta(dataset_ref$path[i])
      if (lab == 'FIN_DEV'){tt = tt %>% rename(iso3 = code)}
      if ('iso3' %in% colnames(tt)){
        tt = tt %>% left_join(country_mapping_by_iso3, by = "iso3")
      } else {
        tt = tt %>% left_join(country_mapping, by = "country")
      }
      tt = tt %>%
        mutate(country_name = ifelse(is.na(country_name), country, country_name),
               year = as.character(year)) %>%
        rename(year_name = year) %>%
        select(country_name, year_name, everything())
      
      matched_country = intersect(unique(df_ref$country_name), unique(tt$country_name))
      # mm = setdiff(unique(df_ref$country_name), unique(tt$country_name))
      # vv = tt %>% select(country, iso3) %>% unique()
      # print(mm);View(vv)
      if (length(matched_country) != uniqueN(df_ref$country_name)){cat('\n-', lab, 'missing countries:',
                                                                       uniqueN(df_ref$country_name) - length(matched_country))}
      missing_years = setdiff(unique(df_ref$year_name), unique(tt$year_name)) %>% sort()
      if (length(missing_years) > 0){cat('\n-', lab, 'missing years:', paste0(missing_years, collapse = ", "))}
      if (sum(is.na(tt$country_name))){cat('\n-', lab, 'missing in country_names')}
      
      
      if (lab == 'FIN_ACC'){
        tt = tt %>% select(-ends_with('_1'), -id, -iso3, -imf_code, -country)
      } else if (lab == 'IFS'){
        tt = tt %>% select(-id, -imf_code, -country)
      } else if (lab == 'COG'){
        tt = tt %>% select(-ccode, -iso3, -ccodecow, -ccodewb, -country, -psi_edate1, -psi_edate2)
      } else if (lab == 'WEF'){
        tt = tt %>% select(-iso3, -country)
      } else if (lab == 'WGI'){
        tt = tt %>% select(-iso3, -country)
      } else if (lab == 'FIN_DEV'){
        tt = tt %>% select(-ifs, -iso3, -country, -imf_region, -imf_income)
      }
      colnames(tt)[-c(1:2)] = paste0(lab, '_', colnames(tt)[-c(1:2)])
      
      # check for duplicated rows (in 'COG')
      check = tt %>%
        group_by(country_name, year_name) %>%
        summarise(COUNT = n(), .groups = 'drop') %>%
        filter(COUNT > 1)
      if (nrow(check) > 0){
        expected_removed_rows = sum(check$COUNT) - nrow(check)
        tt_dim = dim(tt)
        tt_single = tt %>%
          mutate_all(function(x) { attributes(x) <- NULL; x }) %>%
          left_join(check, by = c("country_name", "year_name")) %>%
          filter(is.na(COUNT)) %>%
          select(-COUNT)
        tt_duplicated = suppressWarnings(
          tt %>%
            mutate_all(function(x) { attributes(x) <- NULL; x }) %>%
            left_join(check, by = c("country_name", "year_name")) %>%
            filter(!is.na(COUNT)) %>%
            select(-COUNT) %>%
            group_by(country_name, year_name) %>%
            summarise_all(min, na.rm = T) %>%
            ungroup() %>%
            mutate_all(~(replace(., is.infinite(.), NA)))
        )
        
        tt_single = tt_single %>% select(all_of(colnames(tt)))
        tt_duplicated = tt_duplicated %>% select(all_of(colnames(tt)))
        tt_new = tt_single %>%
          bind_rows(tt_duplicated)
        for (col in colnames(tt)){
          eval(parse(text=paste0('attributes(tt_new$', col,')$label = attributes(tt$', col, ')$label')))
        }
        tt = tt_new
        
        if (tt_dim[1] - nrow(tt) != expected_removed_rows){cat('\n\n #######  mismatch when removing duplicates:', lab)}
        if (tt_dim[2] != ncol(tt)){cat('\n\n #######  mismatch in columns:', lab)}
        rm(tt_single, tt_duplicated, tt_new)
      }
      
      dataset_list[[lab]] = tt
    }
    
    # merge all dataset with no missing imputation and standardization
    df_additional_raw = df_ref
    for (tt in dataset_list){
      df_additional_raw = df_additional_raw %>%
        left_join(tt, by = c("country_name", "year_name"))
    }
    if (nrow(df_additional_raw) != nrow(df_ref)){cat('\n\n #######  rows mismatch in df_additional_raw')}
    
    # replace variable description for lagged variables (_1)
    for (col in df_additional_raw %>% select(ends_with("_1")) %>% colnames()){
      orig_var = gsub("_1", "", col)
      if (!orig_var %in% colnames(df_additional_raw)){
        cat('\n--- unable to recover description for lagged variable:', col)
      } else {
        orig_var_descr = attributes(df_additional_raw %>% pull(orig_var))$label
        col_descr = attributes(df_additional_raw %>% pull(col))$label
        if (substr(col_descr, 1, 3) == "1-Y"){
          eval(parse(text=paste0('attributes(df_additional_raw$', col,')$label = "1-Y lag of ', orig_var_descr,'"')))
        } else {
          cat('\n--- variable seems not to be lag of other variable:', col)
        }
      }
    }
    saveRDS(df_additional_raw, './Data/df_additional_raw.rds')
  }
  df_additional_raw = readRDS('./Data/df_additional_raw.rds')

  # standardise dataset FIRST and THEN impute missing values (with mean for numeric)
  dataset_list_final = list()
  for (tt in dataset_list){
    # remove variable with 1 or less values, isolate dummy
    tt_stats = suppressWarnings(basicStatistics(df_ref %>%
                                                  select(country_name, year_name) %>%
                                                  left_join(tt, by = c("country_name", "year_name")) %>%
                                                  select(-country_name, -year_name))) %>%
      mutate(UNIQUE_VALS = as.numeric(UNIQUE_VALS),
             Min = as.numeric(Min),
             Max = as.numeric(Max))
    tt_remove_var = tt_stats %>%
      filter(UNIQUE_VALS <= 1) %>%
      pull(VARIABLE)
    tt_dummy = tt_stats %>%
      filter(UNIQUE_VALS == 2) %>%
      filter(Min == 0) %>%
      filter(Max == 1) %>%
      pull(VARIABLE)
    tt_pos = tt_stats %>%
      filter(!VARIABLE %in% tt_remove_var) %>%
      filter(!VARIABLE %in% tt_dummy) %>%
      filter(Min >= 0) %>%
      pull(VARIABLE)
    tt_neg = tt_stats %>%
      filter(!VARIABLE %in% tt_remove_var) %>%
      filter(!VARIABLE %in% tt_dummy) %>%
      filter(Min < 0) %>%
      pull(VARIABLE)
    
    # standardize and divide by magnitude if variable has negative values, scale in [0,1] otherwise
    tt_std = tt %>%
      select(-all_of(tt_remove_var)) %>%
      mutate_at(tt_pos, ~(. - min(., na.rm = T))/(max(., na.rm = T) - min(., na.rm = T)) ) %>%
      mutate_at(tt_neg, ~scale(.))
    for (vv in tt_neg){
      val = tt_std %>%
        pull(vv)
      magnitude = 10 ^ ceiling(log10(max(abs(val), na.rm = T)))
      tt_std = tt_std %>%
        mutate(!!sym(vv) := !!sym(vv) / magnitude)
    }
    
    # impute missing by mean if numeric, median if dummy (by country)
    tt_final = tt_std %>%
      group_by(country_name) %>%
      mutate_at(c(tt_pos, tt_neg), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)) %>%
      ungroup() %>%
      group_by(country_name) %>%
      mutate_at(tt_dummy, ~ifelse(is.na(.), median(., na.rm = TRUE), .)) %>%
      ungroup() %>%
      replace(is.na(.), NA)
    
    dataset_list_final = c(dataset_list_final, list(tt_final))
  }
  
  # merge with df_ref
  df_additional = df_ref
  for (tt in dataset_list_final){
    df_additional = df_additional %>%
      left_join(tt, by = c("country_name", "year_name"))
  }
  if (nrow(df_additional) != nrow(df_ref)){cat('\n\n #######  rows mismatch in df_additional')}
  df_additional = df_additional %>%
    rename(country = country_name,
           year = year_name)
  cat('\n--- df_additional (with additional variables) variable range:', df_additional %>% select_if(is.numeric) %>% select(-country, -year) %>% range(na.rm=T))
  
  # extract description for each variable and missing percentage
  lab_table = c()
  for (col in colnames(df_additional)){
    att = attributes(df_additional_raw %>% select(all_of(col)) %>% pull(col))$label
    if (!is.null(att)){
      lab_table = lab_table %>%
        bind_rows(data.frame(Variable = col, Description = att, stringsAsFactors = F))
    }
  }
  lab_table = lab_table %>%
    left_join(
      df_additional %>%
        summarise_all(function(x) round(sum(is.na(x))/ length(x) * 100, 1)) %>%
        t() %>%
        data.frame(Missing_Perc = .) %>%
        rownames_to_column(var = "Variable"), by = "Variable") %>%
    select(Variable, Missing_Perc, Description)
  write.table(lab_table, paste0('./Variable_description_df_additional.csv'), sep = ';', row.names = F)
  rm(df_ref, tt, tt_std, tt_final, dataset_list, dataset_list_final, val, tt_stats, check, lab_table, country_mapping,
     country_mapping_by_iso3, dataset_ref, df_additional_raw)
  saveRDS(df_additional, './Data/df_additional.rds')
  
}   # skip if reloading
df_additional = readRDS('./Data/df_additional.rds')






### execute regression task on test_variable_set with both original regressors and binarized indices
df_set = c("Restricted") #, "RestrictedDifference")
ref_dataset_year_missing_check = "Restricted"    # dataset used to evauale percentage of missing over years for each test_variable_set
recov_set = c("TENSOR_BF" , 'SOFT_IMPUTE')
fitting_set = c('RobPCA', 'DFM_multivar_stacked') #, 'DFM_multivar')   # todo
index_1_set = c(0) #c(-1.5, -1 , -0.5, 0, 0.5, 1, 1.5)
index_2_set = c() #c(-1.5, -1 , -0.5, 0, 0.5, 1, 1.5)
test_variable_set = c('COG_mad_gdppc', # Real GDP per Capita
                      'COG_pwt_gc', #	Share of government consumption at current PPPs
                      'COG_pwt_plcf', # Price level of capital formation, price level of USA GDPo in 2011=1
                      'COG_wdi_trade', # Trade (% of GDP)
                      'COG_wdi_unempilo', #	Unemployment, total (% of total labor force) (modeled ILO)
                      'FIN_ACC_fas83' #	Outstanding Loans, Commercial banks, Domestic Currency (FCSODC_XDC)
)
flag_tuning = T  # tune models parameters. If FALSE, saved parameters will be reloaded
force_tuning = F  # force tuning even if previous tuned parameters are already stored
reload_final_perf = T  # reload final trained model and performance
tuning_strategy = 'bayes'  # 'bayes' or 'grid'
save_all = T  # save parameters and stats
inn_cross_val_fold = 5  # inner cross-validation fold
out_cross_val_fold = 5  # outer cross-validation fold          rimetti 5 e 5 todo:
tuning_criterion = 'rmse'  # tuning criterion and performance measure
algo_set = c('RandomForest', 'GBM', 'ElasticNet', 'SVM-RBF', 'SingleNN', 'MARS') # 'FCNN', 'SVM-Poly'   # todo
df_final = readRDS('./Checkpoints/df_final.rds') %>% mutate(year = as.character(year))
df_additional = readRDS('./Data/df_additional.rds')
res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds') %>% mutate(year = as.character(year))
{
  res_thresh_sensitivity_list = res_thresh_sensitivity_best = res_thresh_sensitivity_residual = error_log = c()
  comb_tot = length(index_1_set) * length(index_2_set)
  start_time = Sys.time()
  sink(paste0('./Log/Threshold_sensitivity_log_', format(Sys.time(), "%Y-%m-%dT%H%M"), '.txt'), append = F, split = T, type = 'output')
  cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
  
  # display missing value ditribution over year for each test_variable_set
  year_missing_check = df_final %>%
    filter(data == ref_dataset_year_missing_check) %>%
    filter(method == recov_set[1]) %>%
    select(country, year) %>%
    unique() %>%
    group_by(year) %>%
    summarise(FullDataset_TotObs = n(), .groups = 'drop') %>%
    arrange(year)
  for (var_target in test_variable_set){
    var_missing = res_ALL_scores %>%     # subset df_additional according to data_match (remove missing observations)
      filter(data == ref_dataset_year_missing_check) %>%
      filter(method == recov_set[1]) %>%
      filter(algo == fitting_set[1]) %>%
      select(country, year) %>%
      unique() %>%
      left_join(df_additional %>% select(country, year, all_of(var_target)), by = c("country", "year")) %>%
      drop_na()
    missing_country = setdiff(unique(df_final$country), unique(var_missing$country))
    missing_country = ifelse(length(missing_country) > 0, length(missing_country), '')
    year_missing_check = year_missing_check %>%
      left_join(var_missing %>%
                  group_by(year) %>%
                  summarise(!!sym(var_target) := n(), .groups = 'drop')%>%
                  mutate(!!sym(var_target) := paste0(!!sym(var_target), ' (', missing_country, ')')), by = "year")
  }
  cat('\n--- Missing observation for each year - total of missing countries in parenthesis:\n')
  print(year_missing_check %>% data.frame())
  
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
            
            
            
            # df_type = 'Restricted'      # todo:
            # recov_met = 'SOFT_IMPUTE'
            # fit_met = 'RobPCA'
            # var_target = 'COG_mad_gdppc'
            # algo_type = 'Spline'

            
            
            # reload settings and parameters
            RDS_lab = paste0('./Checkpoints/ML_fit/', df_type, '_', recov_met, '_', fit_met, '_', var_target, '_', algo_type, '.rds')
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
              algo_type_work = 'CITree'     # in case of random forest, runs CITrees in the binary case (below), random forest otherwise
            } else if (algo_type == 'ElasticNet'){
              algo_type_work = 'GLM'     # in case of ElasticNet, runs GLM in the binary case (below - just 1 regressor!), ElasticNet otherwise
            } else {
              algo_type_work = algo_type
            }
            
            # fit model with original regressors only
            cat('\n\n       °°°° with original regressors')
            # define dataset df_work_orig
            {
              # check for non-missing for all variables - match non-missing year-country for all test_variable_set
              # data_match = res_ALL_scores %>%
              #   filter(data == df_type) %>%
              #   select(country, year) %>%
              #   unique() %>%
              #   left_join(df_additional, by = c("country", "year")) %>%
              #   select(c('country', 'year', test_variable_set)) %>%
              #   group_by(country, year) %>%
              #   summarize_all(function(x) sum(is.na(x))) %>%
              #   ungroup() %>%
              #   mutate(NA_SUM = rowSums(select(., -c(1,2)))) %>%
              #   mutate(KEEP = ifelse(NA_SUM == 0, 1, 0)) %>%
              #   select(country, year, KEEP)
              # cat('  - removed', sum(data_match$KEEP == 0), 'observations from', nrow(data_match), '-', sum(data_match$KEEP == 1), 'remaining\n')
              
              data_work = res_ALL_scores %>%     # subset df_additional according to data_match (remove missing observations)
                filter(data == df_type) %>%
                filter(method == recov_met) %>%
                filter(algo == fit_met) %>%
                select(country, year) %>%
                unique() %>%
                left_join(df_additional %>% select(country, year, all_of(var_target)), by = c("country", "year")) %>%
                drop_na()
                # left_join(df_additional, by = c("country", "year")) %>%
                # select(c('country', 'year', all_of(test_variable_set))) %>%
                # left_join(data_match, by = c("country", "year")) %>%
                # filter(KEEP == 1) %>%
                # select(-KEEP)
              if (sum(is.na(data_work)) != 0){cat('\n ############ data_work: observations with missing not removed')}
              
              df_final_spread = df_final %>%    # df_final in wide format
                filter(data == ref_dataset_year_missing_check) %>%
                filter(method == recov_met) %>%
                spread(variable, val) %>%
                select(-method, -data)
              cat('\n            df_final_spread original variable range:', df_final_spread %>% select(-country, -year) %>% range(na.rm=T))
              
              # df with original variables - TARGET is standardized
              df_work_orig = data_work %>%
                select(country, year, all_of(var_target)) %>%
                rename(TARGET := !!sym(var_target)) %>%
                left_join(df_final_spread, by = c("country", "year")) %>%
                mutate(TARGET = scale(TARGET))
              if (sum(is.na(df_work_orig)) != 0){cat('\n ############ df_work_orig: observations with missing not removed')}
            }
            
            err = try(capture.output(
              regr_orig <- threshold_sensitivity_fit(flag_tuning, force_tuning, reload_final_perf, save_all, inn_cross_val_fold, out_cross_val_fold,
                                                     tuning_criterion, tuning_strategy, df_type, recov_met, fit_met, algo_type,
                                                     algo_type_work = algo_type,   # in case of random forest, runs CITrees in the binary case (below)
                                                     var_target, reload_out, index_1_set, index_2_set,
                                                     res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
                                                     regr_to_test = 'original', 
                                                     index_1_thresh = NULL, 
                                                     index_2_thresh = NULL, 
                                                     df_work_orig = df_work_orig, 
                                                     df_work_index = NULL,
                                                     df_work_rank = NULL,
                                                     df_work_raw = NULL)), silent = T)
            if (class(err) == "try-error"){
              error_log = error_log %>%
                bind_rows(data.frame(data = df_type, method = recov_met, fit_method = fit_met, target_var = var_target, algo_main = algo_type,
                                     regressor_type = 'original', error = attributes(err)$condition$message, stringsAsFactors = F))
              parallelStop()
              cat('\n         ############## error in fitting. Added to error_log ##############')
            } else {
              cat(paste0(err, collapse = '\n'))
              res_thresh_sensitivity_list = regr_orig$res_thresh_sensitivity_list
              res_thresh_sensitivity_best = regr_orig$res_thresh_sensitivity_best
              res_thresh_sensitivity_residual = regr_orig$res_thresh_sensitivity_residual
              reload_out = regr_orig$reload_out
            }
            
            # fit model with index regressors
            comb_count = 1
            for (index_1_thresh in index_1_set){
              index_2_thresh = NULL
              # for (index_2_thresh in index_2_set){
                
                cat('\n\n       °°°° with index regressors, thresholds:', index_1_thresh, 'and', index_2_thresh)#, ' -', round(comb_count / comb_tot * 100, 2), '%')
                # define dataset df_work_index
                {
                  # df with binary index - TARGET is standardized
                  df_work_index = data_work %>%
                    select(country, year, all_of(var_target)) %>%
                    rename(TARGET := !!sym(var_target)) %>%
                    left_join(res_ALL_scores %>%
                                filter(data == df_type) %>%
                                filter(method == recov_met) %>%
                                filter(algo == fit_met) %>%
                                mutate(factor = paste0('Index_', factor)) %>%
                                spread(factor, index) %>%
                                select(-method, -data, -algo, -family, -starts_with('Explain')) %>%
                                mutate(Index_1 = ifelse(Index_1 >= index_1_thresh, 1, 0)),
                                       # Index_2 = ifelse(Index_2 >= index_2_thresh, 1, 0)),
                              by = c("country", "year")) %>%
                    mutate(TARGET = scale(TARGET))
                  if (sum(is.na(df_work_index)) != 0){cat('\n ############ df_work_index: observations with missing not removed')}
                }
                
                n_regressors = df_work_index %>% select(-country, -year, -TARGET) %>% ncol()
                if (algo_type == "SingleNN" & n_regressors == 1){
                  df_work_index = df_work_index %>% mutate(fake = runif(nrow(df_work_index)) / 1e16)
                }
                err = try(capture.output(
                  regr_index <- threshold_sensitivity_fit(flag_tuning, force_tuning, reload_final_perf, save_all, inn_cross_val_fold, out_cross_val_fold,
                                                          tuning_criterion, tuning_strategy, df_type, recov_met, fit_met, algo_type, algo_type_work,
                                                          var_target, reload_out, index_1_set, index_2_set,
                                                          res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
                                                          regr_to_test = 'index', 
                                                          index_1_thresh = index_1_thresh, 
                                                          index_2_thresh = index_2_thresh, 
                                                          df_work_orig = NULL, 
                                                          df_work_index = df_work_index,
                                                          df_work_rank = NULL,
                                                          df_work_raw = NULL)), silent = T)
                if (class(err) == "try-error"){
                  error_log = error_log %>%
                    bind_rows(data.frame(data = df_type, method = recov_met, fit_method = fit_met, target_var = var_target, algo_main = algo_type,
                                         regressor_type = 'index', error = attributes(err)$condition$message, stringsAsFactors = F))
                  parallelStop()
                  cat('\n         ############## error in fitting. Added to error_log ##############')
                } else {
                  cat(paste0(err, collapse = '\n'))
                  res_thresh_sensitivity_list = regr_index$res_thresh_sensitivity_list
                  res_thresh_sensitivity_best = regr_index$res_thresh_sensitivity_best
                  res_thresh_sensitivity_residual = regr_index$res_thresh_sensitivity_residual
                  reload_out = regr_index$reload_out
                }
                comb_count = comb_count + 1
                
              # } # index_2_thresh
            } # index_1_thresh
            
            # fit model with ranked index regressors
            {
              # cat('\n\n       °°°° with ranked index regressors')
              # define dataset df_work_rank 
              {
                # # df with binary index - TARGET is standardized
                # df_work_rank = data_work %>%
                #   select(country, year, all_of(var_target)) %>%
                #   rename(TARGET := !!sym(var_target)) %>%
                #   left_join(res_ALL_scores %>%
                #               filter(data == df_type) %>%
                #               filter(method == recov_met) %>%
                #               filter(algo == fit_met) %>%
                #               mutate(factor = paste0('Index_', factor)) %>%
                #               spread(factor, index) %>%
                #               select(-method, -data, -algo, -family, -starts_with('Explain')) %>%
                #               mutate(Index_1 = cut(Index_1, breaks = c(-Inf,index_1_set,Inf)),
                #                      Index_2 = cut(Index_2, breaks = c(-Inf,index_2_set,Inf))),
                #             by = c("country", "year")) %>%
                #   mutate(TARGET = scale(TARGET))
                # if (sum(is.na(df_work_rank)) != 0){cat('\n ############ df_work_rank: observations with missing not removed')}
              }
              # 
              # regr_rank = threshold_sensitivity_fit(flag_tuning, force_tuning, reload_final_perf, save_all, inn_cross_val_fold, out_cross_val_fold,
              #                                       tuning_criterion, tuning_strategy, df_type, recov_met, fit_met, algo_type,
              #                                       algo_type_work = algo_type,    # in case of random forest, runs CITrees in the binary case (above)
              #                                       var_target, reload_out, index_1_set, index_2_set,
              #                                       res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
              #                                       regr_to_test = 'rank_index', 
              #                                       index_1_thresh = index_1_thresh, 
              #                                       index_2_thresh = index_2_thresh, 
              #                                       df_work_orig = NULL, 
              #                                       df_work_index = NULL,
              #                                       df_work_rank = df_work_rank,
              #                                       df_work_raw = NULL)
              # 
              # res_thresh_sensitivity_list = regr_rank$res_thresh_sensitivity_list
              # res_thresh_sensitivity_best = regr_rank$res_thresh_sensitivity_best
              # res_thresh_sensitivity_residual = regr_rank$res_thresh_sensitivity_residual
              # reload_out = regr_rank$reload_out
            }
            
            # fit model with raw index regressors
            cat('\n\n       °°°° with raw index regressors')
            # define dataset df_work_raw 
            {
              # df with continous index - TARGET is standardized
              df_work_raw = data_work %>%
                select(country, year, all_of(var_target)) %>%
                rename(TARGET := !!sym(var_target)) %>%
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
            
            n_regressors = df_work_raw %>% select(-country, -year, -TARGET) %>% ncol()
            if (algo_type == "ElasticNet" & n_regressors == 1){
              algo_type_raw = algo_type_work
            } else {
              algo_type_raw = algo_type
            }
            if (algo_type == "SingleNN" & n_regressors == 1){
              df_work_raw = df_work_raw %>% mutate(fake = runif(nrow(df_work_raw)) / 1e16)
            }
            err = try(capture.output(
              regr_raw <- threshold_sensitivity_fit(flag_tuning, force_tuning, reload_final_perf, save_all, inn_cross_val_fold, out_cross_val_fold,
                                                    tuning_criterion, tuning_strategy, df_type, recov_met, fit_met, algo_type,
                                                    algo_type_work = algo_type_raw,    # in case of random forest, runs CITrees in the binary case (above)
                                                    var_target, reload_out, index_1_set, index_2_set,
                                                    res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
                                                    regr_to_test = 'raw_index', 
                                                    index_1_thresh = index_1_thresh, 
                                                    index_2_thresh = index_2_thresh, 
                                                    df_work_orig = NULL, 
                                                    df_work_index = NULL,
                                                    df_work_rank = NULL,
                                                    df_work_raw = df_work_raw)), silent = T)
            if (class(err) == "try-error"){
              error_log = error_log %>%
                bind_rows(data.frame(data = df_type, method = recov_met, fit_method = fit_met, target_var = var_target, algo_main = algo_type,
                                     regressor_type = 'raw_index', error = attributes(err)$condition$message, stringsAsFactors = F))
              parallelStop()
              cat('\n         ############## error in fitting. Added to error_log ##############')
            } else {
              cat(paste0(err, collapse = '\n'))
              res_thresh_sensitivity_list = regr_raw$res_thresh_sensitivity_list
              res_thresh_sensitivity_best = regr_raw$res_thresh_sensitivity_best
              res_thresh_sensitivity_residual = regr_raw$res_thresh_sensitivity_residual
              reload_out = regr_raw$reload_out
            }
            
            
            if (save_all){
              cat('\n --- saving RDS')
              saveRDS(reload_out, RDS_lab)
              saveRDS(res_thresh_sensitivity_residual, './Checkpoints/res_thresh_sensitivity_residual.rds')
              saveRDS(res_thresh_sensitivity_best, './Checkpoints/res_thresh_sensitivity_best.rds')
              write.table(res_thresh_sensitivity_list, "./Results/3_threshold_sensitivity_performance_list.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
              write.table(res_thresh_sensitivity_best, "./Results/3_threshold_sensitivity_performance_best.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
              if (length(error_log) > 0){
                write.table(error_log, "./Results/3_threshold_sensitivity_error_log.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
              }
            }
            
          } # algo_type
          var_count = var_count + 1
        } # var_target
      } # fit_met
    } # recov_met
  } # df_type
  cat('\n\nTotal elapsed time: ', round((as.numeric(Sys.time())-as.numeric(start_time)) / 60), 'mins ')
  cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
  sink()
  
  if (nrow(unique(res_thresh_sensitivity_residual)) != nrow(res_thresh_sensitivity_residual)){cat('\n\n\n\n\n ##################### duplicates in res_thresh_sensitivity_residual\n\n\n\n')}
  if (length(error_log) > 0){
    cat('\n\n---', nrow(error_log), 'model did not converge:\n')
    print(error_log)
  }
}



# plot comparison of performance on original variable vs raw index and binary index
df_set = c('Restricted') #, 'RestrictedDifference')      # 1 graph for each df_type + recov_met + target_variable, with fit_method on rows and algo on columns
recov_set = c('TENSOR_BF') #, 'SOFT_IMPUTE')
fitting_set = data.frame(original = c('RobPCA', 'DFM_multivar_stacked'),    # select method to plot and their label
                         plot_label = c('Robust PCA', 'DFM'), stringsAsFactors = F)
delta_to_plot = c("Original vs continuous") #, "Original vs binary")
algo_to_plot = c('RandomForest', 'ElasticNet', 'SVM-RBF', 'SingleNN', 'MARS') # , 'GBM'
algo_rename = data.frame(algo_orig = c('GLM', 'CITree', 'RandomForest', 'GBM', 'ElasticNet', 'SVM-RBF', 'SingleNN', 'MARS'),
                         algo_new = c('OLS', 'C.I. Tree', 'Random Forest', 'GBM', 'Elastic-Net', 'SVM-RBF', 'Single Layer NN', 'MARS'), stringsAsFactors = F)
variable_description = data.frame(target_var = c('COG_mad_gdppc', 'COG_pwt_gc', 'COG_pwt_plcf', 'COG_wdi_trade', 'COG_wdi_unempilo', 'FIN_ACC_fas83'),
                                  Description = c('Real GDP per Capita', 'Share of government consumption', 'Price level of capital formation',
                   'Trade volume', 'Unemployment Rate', 'Outstanding Loans of Commercial banks'), stringsAsFactors = F)
{
  res_thresh_sensitivity_residual = readRDS('./Checkpoints/res_thresh_sensitivity_residual.rds') %>%
    left_join(fitting_set, by = c("fit_method" = "original")) %>%
    select(-fit_method) %>%
    rename(fit_method = plot_label)
  res_thresh_sensitivity_best = readRDS('./Checkpoints/res_thresh_sensitivity_best.rds') %>%
    left_join(fitting_set, by = c("fit_method" = "original")) %>%
    select(-fit_method) %>%
    rename(fit_method = plot_label)
  
  # create fake prediction object to evaluate performance metric
  n = getTaskSize(bh.task)
  train.set = seq(1, n, by = 2)
  test.set = seq(2, n, by = 2)
  lrn = makeLearner("regr.gbm", n.trees = 2)
  mod = train(lrn, bh.task, subset = train.set)
  prediction_template = predict(mod, task = bh.task, subset = test.set)
  
  for (df_type in df_set){
    for (recov_met in recov_set){
      
      data_check = res_thresh_sensitivity_best %>%
        filter(fit_method %in% fitting_set$plot_label) %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        group_by(fit_method, target_var, algo_main) %>%
        summarize(Regressor_type_count = n(), .groups = "drop")
      wrong_combination = data_check %>%
        filter(Regressor_type_count != 3)
      if (nrow(wrong_combination) > 0){
        cat('\n ##########  ', df_type, '-', recov_met, 'has missing or duplicated combination of fit_method, target_var, algo_main and regressor_type')
        cat('\n - dropping:\n', wrong_combination %>% mutate(space = '   ') %>% select(space, everything(), -Regressor_type_count) %>%
              unite("lab", sep = " - ") %>% unlist() %>% paste0('\n'))
        data_check = data_check %>%
          filter(Regressor_type_count == 3)
      }
      
      # create data to plot
      performance_measure = unique(res_thresh_sensitivity_best$PERF)
      if (length(performance_measure) > 1){
        cat('\n ##########  ', df_type, '-', recov_met, 'has more than one performance metric:', paste0(performance_measure, collapse = ", "),
            ' - setting', performance_measure[1])
        performance_measure = performance_measure[1]
      }
      data_residuals = res_thresh_sensitivity_residual %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        left_join(data_check, by = c("fit_method", "target_var", "algo_main")) %>%
        filter(!is.na(Regressor_type_count)) %>%
        select(-Regressor_type_count)
      tot_years = data_residuals$year %>% unique() %>% sort()
      
      data_plot = res_thresh_sensitivity_best %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        left_join(data_check, by = c("fit_method", "target_var", "algo_main")) %>%
        filter(!is.na(Regressor_type_count)) %>%
        select(fit_method, target_var, target_var_avg, regressor_type, index_1_thresh, index_2_thresh, algo_main, ALGO, N_OBS, FULL_SET_PERF) %>%
        mutate(LINE = "MEAN",
               year = tot_years[1])
      
      year_performance = c()
      for (i in 1:nrow(data_check)){    # loop all combinations of fit_method, target_var, algo_main
        tt_year = data_residuals %>%
          left_join(data_check[i, ], by = c("fit_method", "target_var", "algo_main")) %>%
          filter(!is.na(Regressor_type_count)) %>%
          left_join(data_plot %>% select(fit_method, target_var, regressor_type, target_var_avg, N_OBS) %>% unique(),
                    by = c("fit_method", "target_var", "regressor_type"))
        for (reg_type in unique(tt_year$regressor_type)){
          for (yr in tot_years){
            tt = tt_year %>%
              filter(regressor_type == reg_type) %>%
              filter(year == yr)
            prediction_tt = prediction_template
            prediction_tt$data = tt %>% select(truth, response)
            year_performance = year_performance %>%
              bind_rows(tt %>%
                          select(fit_method, target_var, regressor_type, algo_main, algo, index_1_thresh, index_2_thresh, year, target_var_avg, N_OBS) %>%
                          rename(ALGO = algo) %>%
                          unique() %>%
                          bind_cols(
                            data.frame(FULL_SET_PERF = performance(prediction_tt, eval(parse(text=performance_measure))),
                                       LINE = 'YEARS', stringsAsFactors = F)))
          } # yr
        } # reg_type
      } # i
      
      # duplicate performance on all set for each year (will plot a straight line)
      for (yr in tot_years[-1]){
        data_plot = data_plot %>%
          bind_rows(data_plot %>% filter(year == tot_years[1]) %>% mutate(year = yr))
      }
      data_plot = data_plot %>%
        bind_rows(year_performance)
      if (nrow(data_plot) != 2 * nrow(data_check) * 3 * length(tot_years)){cat('\n ##########  ', df_type, '-', recov_met, 'data_plot has wrong dimension')}
      
      # evaluate performance difference to plot ( binary/raw - original ) / original
      data_plot = data_plot %>%
        left_join(algo_rename, by = c("algo_main" = "algo_orig")) %>%
        rename(algo_main_orig = algo_main,
               algo_main = algo_new) %>%
        left_join(algo_rename, by = c("ALGO" = "algo_orig")) %>%
        select(-ALGO) %>%
        rename(ALGO = algo_new) %>%
        group_by(fit_method, target_var, algo_main, LINE, year) %>%
        mutate(`Original vs binary` = (FULL_SET_PERF[which(regressor_type == 'index')] - FULL_SET_PERF[which(regressor_type == 'original')]) / FULL_SET_PERF[which(regressor_type == 'original')],
               `Original vs continuous` = (FULL_SET_PERF[which(regressor_type == 'raw_index')] - FULL_SET_PERF[which(regressor_type == 'original')]) / FULL_SET_PERF[which(regressor_type == 'original')],
               Original_PERF = FULL_SET_PERF[which(regressor_type == 'original')],
               Binary_PERF = FULL_SET_PERF[which(regressor_type == 'index')],
               Continuous_PERF = FULL_SET_PERF[which(regressor_type == 'raw_index')],
               Algo_label = paste0('Original: ', ALGO[which(regressor_type == 'original')],
                                   ifelse("Original vs binary" %in% delta_to_plot, paste0('\nBin. ind.: ', ALGO[which(regressor_type == 'index')]), ''),
                                   ifelse("Original vs continuous" %in% delta_to_plot, paste0('\nCont. ind.: ', ALGO[which(regressor_type == 'raw_index')]), ''), collapse = ""),
               Index1_label = ifelse(length(unique(setdiff(index_1_thresh, 'None'))) == 0, NA, paste0('Index 1 thresh: ', unique(setdiff(index_1_thresh, 'None')))),
               Index2_label = ifelse(length(unique(setdiff(index_2_thresh, 'None'))) == 0, NA, paste0('Index 2 thresh: ', unique(setdiff(index_2_thresh, 'None'))))) %>%
        ungroup() %>%
        select(-regressor_type, -index_1_thresh, -index_2_thresh, -ALGO, -FULL_SET_PERF) %>%
        unique() %>%
        # mutate(Index2_label = 'prova') %>%
        rowwise() %>%
        mutate(Facet_label = paste0(
          # paste0('Model: ', algo_main), '\n',
          Algo_label, '\n',
          ifelse("Original vs binary" %in% delta_to_plot, paste0(na.omit(c(Index1_label, Index2_label)),
                                                                 collapse = '\n'), ''))) %>%
        select(-Index1_label, -Index2_label, -Algo_label) %>%
        gather('Delta Type:', 'y_val', c("Original vs binary", "Original vs continuous")) %>%
        # left_join(read.csv("Variable_description_df_additional.csv", sep=";", stringsAsFactors=FALSE) %>% select(-Missing_Perc), by = c("target_var" = "Variable")) %>%
        left_join(variable_description, by = "target_var") %>%
        filter(`Delta Type:` %in% delta_to_plot)
      
      # save statistics for table
      write.table(data_plot %>%
                    filter(LINE == "MEAN") %>%
                    filter(year == tot_years[1]) %>%
                    select(-LINE, -year, -Facet_label),
                  paste0("./Results/3b_prediction_table_", df_type, "_", recov_met, ".csv"), row.names = F, sep = ";", col.names = T)
      
      # loop target_var
      y_lim = data_plot %>% filter(algo_main_orig %in% algo_to_plot) %>% pull(y_val) %>% range(na.rm = T)
      for (var in unique(data_plot$target_var)){
        tt = data_plot %>%
          filter(algo_main_orig %in% algo_to_plot) %>%
          filter(target_var == var) %>%
          mutate(LINE = factor(ifelse(LINE == "MEAN", 'year average', 'single year'))) %>%
          rename(`Performance on:` = LINE)
        
        # pp = tt %>% filter(fit_method == 'RobPCA' & algo_main == "GBM")
        png(paste0("./Results/3b_prediction_", df_type, "_", recov_met, "_", var, ".png"), width = 19, height = 12, units = 'in', res=300)
        if (length(delta_to_plot) > 1){
          p = ggplot(data=tt, aes(x=year, y=y_val, color=`Delta Type:`, group=interaction(`Performance on:`, `Delta Type:`))) +
            geom_line(aes(linetype=`Performance on:`), size=1.7)
        } else {
          p = ggplot(data=tt, aes(x=year, y=y_val, group=`Performance on:`)) +
            geom_line(aes(linetype=`Performance on:`), size=1.7, color = 'blue')
        }
        plot(p +
               ylab(paste0('Relative ', toupper(performance_measure), ' from original variables')) +
               scale_colour_manual(values=c('blue3', 'darkorange')) +
               scale_y_continuous(labels = scales::percent, limits = y_lim) +
               facet_grid(fit_method ~ Facet_label, scales = 'fixed', switch="y") +
               ggtitle(paste0('Target variable: ', unique(tt$Description)),
                       subtitle = paste0('Average value: ', unique(tt$target_var_avg) %>% gsub("Â", "", .) %>% gsub("±", " ± ", .),
                                         '\nTotal observations: ', unique(tt$N_OBS), '\n')) +
               theme(axis.title=element_text(size=29),
                     axis.text.x=element_text(size=17, angle = 45, vjust=1, hjust=1),
                     axis.text.y=element_text(size=17),
                     legend.text=element_text(size=22),
                     legend.title=element_text(size=22),
                     legend.key.size = unit(1, "cm"),
                     legend.key.width = unit(1.7,"cm"),
                     legend.position="top",
                     plot.title = element_text(size=32),
                     plot.subtitle = element_text(size=25),
                     panel.background = element_rect(fill = "white", colour = "black"),
                     panel.grid = element_line(colour = "grey", linetype = 'dashed', size = 0.5),
                     strip.text = element_text(size = 19)) +
               guides(linetype = guide_legend(override.aes = list(size = 2.2)),
                      colour = guide_legend(override.aes = list(size = 2.2))))
        dev.off()
      } # var
    } # recov_met
  } # df_type
  rm(data_plot, data_check, data_residuals, tt, tt_year, year_performance, wrong_combination, prediction_template, prediction_tt, lrn, mod)
}



### table and figure for Latex
{
  # variable statistics
  {
    df = readRDS('./Checkpoints/df.rds')
    stats = basicStatistics(df) %>%
      left_join(read.csv("Variable_description_Restricted_Dataset.csv", sep=";", stringsAsFactors=FALSE), by = c("VARIABLE" = "var")) %>%
      filter(!VARIABLE %in% c("country", "year")) %>%
      filter(!is.na(label)) %>%
      mutate(vv = gsub("var", "", VARIABLE),
             Frequency = "Yearly",
             `Coefficient of Variation` = as.numeric(StDev) / as.numeric(Mean)) %>%
      arrange(as.numeric(vv)) %>%
      select(VARIABLE, label, Frequency, NUM_OSS, NAs, Min, Max, Mean, Median, StDev, `Coefficient of Variation`) %>%
      mutate_at(vars(Min, Max, Mean, Median, StDev, `Coefficient of Variation`), function(x) round(as.numeric(x), 2)) %>%
      mutate_at(vars(NUM_OSS, Min, Max, Mean, Median, StDev, `Coefficient of Variation`), function(x) format(as.numeric(x), big.mark = ",")) %>%
      rename(Variable = VARIABLE,
             Description = label,
             `Total Observations` = NUM_OSS,
             `Missing Values` = NAs,
             `Standard Deviation` = StDev)
    write.table(stats, './Paper/Latex_Table_Figure/00_stats_variable.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    
    country_missing = read.csv("./Stats/01_final_df_missing_by_country.csv", sep=";", stringsAsFactors=FALSE) %>%
      mutate(`Missing Values` = paste0(format(Tot_NA, big.mark=","), ' (', NA_perc, '%)')) %>%
      rename(Country = country) %>%
      select(Country, `Missing Values`)
    total_columns = 3
    n_nrow_split = ceiling(nrow(country_missing) / total_columns)
    for (i in 1:total_columns){
      write.table(country_missing[((i - 1)*n_nrow_split + 1):min(c(i*n_nrow_split, nrow(country_missing))), ],
                  paste0('./Paper/Latex_Table_Figure/01_stats_country_', i, '.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    }
    
   year_missing = read.csv("./Stats/01_final_df_missing_by_year.csv", sep=";", stringsAsFactors=FALSE) %>%
     mutate(`Total Observations` = format(Tot_non_NA + Tot_NA, big.mark=",")) %>%
     mutate(`Missing Values` = paste0(format(Tot_NA, big.mark=","), ' (', NA_perc, '%)')) %>%
     rename(Year = year) %>%
     select(Year, `Total Observations`, `Missing Values`)
   write.table(year_missing, './Paper/Latex_Table_Figure/02_stats_year.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    
  }
  
  # correlation matrix
  {
    df_type = 'Restricted'
    recov_met = 'TENSOR_BF'
    
    corr_res = evaluate_correlation(df_final %>%    # df_final in wide format
                           filter(data == df_type) %>%
                           filter(method == recov_met) %>%
                           spread(variable, val) %>%
                           select(-method, -data, -country, -year)) %>%
      mutate(star = ifelse(Corr_pVal <= 0.05, '*', ''),
             label = paste0(round(Corr, 2), star)) %>%
      rowwise() %>%
      mutate(ref = paste0(c(sort(as.numeric(c(gsub("var", "", Var1), gsub("var", "", Var2))))), collapse = ","))
    tot_vars = corr_res %>% select(Var1, Var2) %>% unlist() %>% unique()
    corr_mat = matrix("", ncol = length(tot_vars)-1, nrow = length(tot_vars)-1)
    for (i in 1:length(tot_vars)){
      for (j in 1:length(tot_vars)){
        if (i > j){
          corr_mat[i-1, j] = corr_res %>% filter(ref == paste0(sort(c(i, j)), collapse = ",")) %>% pull(label)
        }
      }
    }
    corr_mat = data.frame(corr_mat) %>%
      `colnames<-`(paste0("var", 1:(length(tot_vars)-1))) %>%
      `rownames<-`(paste0("var", 2:length(tot_vars))) %>%
      rownames_to_column(var = 'xx') %>%
      select(xx, everything())
    write.table(corr_mat, './Paper/Latex_Table_Figure/00_correlation_matrix.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  }
  
  # method performance
  {
    df_type = 'Restricted'
    recov_met = 'TENSOR_BF'
    factors = 1
    DFM_method = "DFM_multivar_stacked"
    VAR_alpha = 0.2
    kalm_Q_hat_mode = 'identity'
    
    res_PCA_stats = readRDS('./Checkpoints/res_PCA_stats.rds')
    latex_PCA_stats = res_PCA_stats %>%
      filter(method == recov_met) %>%
      filter(data == df_type) %>%
      filter(year != 'Avg') %>%
      filter(PC <= factors) %>%
      rename(Method = PCA,
             `Number of PC` = PC) %>%
      group_by(Method, `Number of PC`) %>%
      summarise(`Mean Explained Variance` = paste0('$', round(mean(Explain_var_Loadings) * 100, 1), '\\pm', round(sd(Explain_var_Loadings) * 100, 1), '\\%$'),
                `Mean $R^2$` = paste0('$', round(mean(Explain_var) * 100, 1), '\\pm', round(sd(Explain_var) * 100, 1), '\\%$'),
                `Mean $R^2$ on 99th` = paste0('$', round(mean(Explain_var_99) * 100, 1), '\\pm', round(sd(Explain_var_99) * 100, 1), '\\%$'),
                `Mean $R^2$ on 95th` = paste0('$', round(mean(Explain_var_95) * 100, 1), '\\pm', round(sd(Explain_var_95) * 100, 1), '\\%$'), .groups = "drop"
      )
    write.table(latex_PCA_stats, './Paper/Latex_Table_Figure/03_PCA_stats.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".", quote = F)
    
    # DFM_method = "DFM_multivar_stacked"
    VAR_alpha = 0.2
    kalm_Q_hat_mode = 'identity'
    
    res_DFM_stats = readRDS(paste0('./Checkpoints/res_DFM_stats_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
    latex_DFM_stats = res_DFM_stats %>%
      filter(method == recov_met) %>%
      filter(data == df_type) %>%
      filter(DFM != "DFM_univar") %>%
      filter(Total_Factors <= factors) %>%
      rename(Method = DFM,
             `Number of Factors` = Total_Factors) %>%
      # mutate(Method = 'DFM') %>%
      group_by(Method, `Number of Factors`) %>%
      summarise(`$R^2$` = paste0('$', round(mean(Explain_var) * 100, 1), '\\%$'),
                `$R^2$ on 99th` = paste0('$', round(mean(Explain_var_99) * 100, 1), '\\%$'),
                `$R^2$ on 95th` = paste0('$', round(mean(Explain_var_95) * 100, 1), '\\%$'), .groups = "drop"
      )
    write.table(latex_DFM_stats, './Paper/Latex_Table_Figure/04_DFM_stats.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".", quote = F)
    
    res_PCA_list = readRDS('./Checkpoints/res_PCA_list.rds')
    max_variance = 90
    # scree plot all method
    {
      # n_row = length(names(res_PCA_list[[df_type]][[1]][[1]]))    # #_PCA_meth
      # n_col = length(names(res_PCA_list[[df_type]][[1]]))     # number of years
      # 
      # row_list = list()                                           # row_list[[j]][[i]]  is the element i,j in the graph
      # for (yr in setdiff(names(res_PCA_list[[df_type]][[1]]), 'Avg')){
      #   i = 1
      #   for (pc_met in names(res_PCA_list[[df_type]][[1]][[1]])){
      #     scree_reload = res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_final]][['scree_plot']]
      #     scree_reload$layers[[4]]$aes_params$size = 12
      #     scree_reload$layers[[4]]$aes_params$hjust = 0.5
      #     p = scree_reload +
      #       ylab('Explained Variance (%)') + 
      #       xlab('PC') +
      #       theme(axis.title.y=element_text(size=35),
      #             axis.title.x=element_text(size=35),
      #             plot.title = element_text(size=35),
      #             axis.text = element_text(size = 35)) +
      #       ggtitle(paste0(pc_met, ' - ', yr)) +
      #       ylim(0, max_variance + 10)
      #     
      #     if (yr != setdiff(names(res_PCA_list[[df_type]][[1]]), 'Avg')[1]){p = p + theme(axis.text.y = element_blank(),
      #                                                                                     axis.title.y = element_blank())}
      #     if (i < n_row){p = p + theme(axis.title.x = element_blank())}
      #     row_list[[yr]][[i]] = ggplotGrob(p)
      #     i = i + 1
      #   } # pc_met
      # } # yr
      # 
      # col_list = list()
      # for (i in c(1:n_col)){
      #   col_list[[i]] = do.call(rbind, c(row_list[[i]], size="last"))
      # }
      # g = do.call(cbind, c(col_list, size="last"))
      # png(paste0('./Paper/Latex_Table_Figure/Scree_plot_full.png'), width = 55, height = 20, units = 'in', res=300)
      # grid.draw(g)
      # dev.off()
    }
    
    # scree plot single PCA
    n_col = 5
    {
      tot_years = length(res_PCA_list[[df_type]][[1]])
      for (pc_met in names(res_PCA_list[[df_type]][[1]][[1]])){
        row_list = list()                                           # row_list[[j]][[i]]  is the element i,j in the graph
        i = 0
        for (yr in setdiff(names(res_PCA_list[[df_type]][[1]]), 'Avg')){
          scree_reload = res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]][['scree_plot']]
          scree_reload$layers[[4]]$aes_params$size = 16
          scree_reload$layers[[4]]$aes_params$hjust = 0.5
          p = scree_reload +
            ylab('Explained Variance (%)') + 
            xlab('PC') +
            theme(axis.title.y=element_text(size=55, margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 20)),
                  axis.title.x=element_text(size=55),
                  plot.title = element_text(size=55),
                  axis.text = element_text(size = 45)) +
            ggtitle(yr) +
            ylim(0, max_variance + 10)
          
          xi = i%/%n_col + 1
          xj = i%%n_col + 1
          if (xj != 1){p = p + theme(axis.text.y = element_blank(),
                                     axis.title.y = element_blank())}
          if ((xi < ceiling(tot_years/n_col) & xj <= tot_years%%n_col) | (xj > tot_years%%n_col & xi < ceiling(tot_years/n_col)-1)){
            p = p + theme(axis.title.x = element_blank())}
          row_list[[toString(xj)]][[xi]] = ggplotGrob(p)
          i = i + 1
        } # yr
        if (xj != n_col){
          for (res in (i%%n_col + 1):n_col){row_list[[toString(res)]][[i%/%n_col+1]] = ggplotGrob(ggplot() + theme(panel.background = element_blank()))}
        }
        
        col_list = list()
        for (i in c(1:n_col)){
          col_list[[i]] = do.call(rbind, c(row_list[[i]], size="last"))
        }
        if (pc_met == 'PCA'){main_title = 'PCA'}
        if (pc_met == 'RobPCA'){main_title = 'Robust PCA'}
        if (pc_met == 'RobSparPCA'){main_title = 'Robust Sparse PCA'}
        g = do.call(cbind, c(col_list, size="last"))
        g = gtable_add_padding(g, unit(c(3,1,1,1), "cm")) # t,r,b,l
        g = gtable_add_grob(
          g,
          textGrob(main_title, gp=gpar(fontsize=75)),
          1,1,1,ncol(g))
        png(paste0('./Paper/Latex_Table_Figure/Scree_plot_', pc_met, '.png'), width = 45, height = 11 * xi, units = 'in', res=300)
        grid.draw(g)
        dev.off()
      } # pc_met
    }
    
    # loading for all PCA
    signif_thresh_plot = 0.2 # used to highlight significative loadings
    leading_var = "var1"
    leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
    PC_to_compare = 1
    {
      res_PCA_loadings = readRDS('./Checkpoints/res_PCA_loadings.rds')
      res_PCA_importance = readRDS('./Checkpoints/res_PCA_importance.rds')
      add_recov = c()#ifelse(length(unique(res_PCA_loadings$method)) >= 2, 'DIFFERENCE', c())
      y_range = range((res_PCA_loadings %>% filter(PC <= PC_to_compare))$loading); y_range[2] = y_range[2] * 1.1 # extra space for avg_expl_var_lab
      n_year = uniqueN(res_PCA_loadings$year) - ifelse('Avg' %in% res_PCA_loadings$year, 1, 0)

          # adjust sign according to leading variable
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

          
          row_list = list()
          for (pc in 1:PC_to_compare){
            
            # if (pc == 1){grad = colorRampPalette(c('#f0f0f0', 'black'))}  # greyscale
            grad = colorRampPalette(c('#deebf7', '#3182bd'))  # bluescale
            # if (pc == 3){grad = colorRampPalette(c('#fee6ce', '#e6550d'))}  # redscale
            # if (pc == 4){grad = colorRampPalette(c('#e5f5e0', '#31a354'))}  # greenscale
            tit_lab = paste0("<span style='font-size:35'><b>", paste0('PC ', pc, '  '), "</b><span style='font-size:30'>",
                             (avg_expl_var %>% filter(PC == paste0('PC', pc)))$LAB)
            tit_lab <- rich_text_grob(tit_lab,
                                      x = unit(-44, "lines"),#-0.62 * uniqueN(res_PCA_loadings_adj$variable) * uniqueN(res_PCA_loadings_adj$year), "lines"), # for 17*8 should be 95
                                      y = unit(72, "lines")) # for 2 PC should be 33    83 -25 * PC_to_compare

            p = ggplot(res_PCA_loadings_adj %>%
                         filter(data == df_type) %>%
                         filter(method == recov_met) %>%
                         filter(PC == pc) %>%
                         mutate(variable = split_string(variable, 15)) %>%
                         mutate(year = as.factor(year), PC = paste0('PC ', PC),
                                variable = factor(variable, levels=paste0('var', 1:uniqueN(variable)))), aes(fill=year, y=loading, x=variable)) + 
              geom_bar(position="dodge", stat="identity") +
              facet_wrap(~ PCA, ncol = 1, dir = 'v', scales = 'free_y', strip.position = 'left') +
              scale_fill_manual(values = c(rev(grad(n_year)), 'darkgoldenrod1')) +
              ylim(y_range[1], y_range[2]) +
              geom_vline(xintercept = c(1:(uniqueN(res_PCA_loadings_adj$variable) - 1)) + 0.5) +
              # geom_hline(aes(yintercept = -signif_thresh_plot), color = 'red', linetype = 'twodash', size = 0.4) +
              # geom_hline(aes(yintercept = signif_thresh_plot), color = 'red', linetype = 'twodash', size = 0.4) +
              theme(axis.title.y=element_blank(),
                    axis.title.x=element_blank(),
                    axis.text.x = element_text(size = 24, vjust=0.3),
                    axis.text.y=element_text(size=20),
                    legend.text=element_text(size=28),
                    legend.title=element_text(size=28),
                    legend.key.size = unit(1.3, "cm"),
                    strip.text.y = element_text(size = 24,face="bold"),
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
          png('./Paper/Latex_Table_Figure/PCA_loadings.png', width = 22, height = 17, units = 'in', res=300)
          grid.draw(g)
          dev.off()

    }
    
    # loading for all DFM
    df_type = 'Restricted'
    recov_met = 'TENSOR_BF'
    VAR_alpha = 0.2
    kalm_Q_hat_mode = 'identity'
    dfm_to_plot = data.frame(DFM = c("DFM_multivar_stacked", "DFM_multivar"), label = c("DFM without interactions", "DFM with interactions"), stringsAsFactors = F)
    factor_to_compare = 1
    {
      res_DFM_loadings = readRDS(paste0('./Checkpoints/res_DFM_loadings_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
      res_DFM_best_model = readRDS(paste0('./Checkpoints/res_DFM_best_model_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
      
      # create main title
      avg_expl_var = res_DFM_best_model %>%
        filter(method == recov_met & data == df_type & DFM %in% dfm_to_plot$DFM & Total_Factors <= factor_to_compare) %>%
        left_join(dfm_to_plot, by = "DFM") %>%
        select(-DFM) %>%
        rename(DFM = label) %>%
        group_by(DFM, Total_Factors) %>%
        mutate(LAB = paste0(DFM, ' ', round(Explain_var_95 * 100), '%')) %>%
        group_by(DFM, Total_Factors) %>%
        summarise(LAB = paste0(LAB, collapse = '  |  '), .groups = 'drop')
      
      row_list = list()
      for (pc in 1:factor_to_compare){
        
        # if (pc == 1){grad = colorRampPalette(c('#f0f0f0', 'black'))}  # greyscale
        # grad = colorRampPalette(c('#deebf7', '#3182bd'))  # bluescale
        # if (pc == 3){grad = colorRampPalette(c('#fee6ce', '#e6550d'))}  # redscale
        # if (pc == 4){grad = colorRampPalette(c('#e5f5e0', '#31a354'))}  # greenscale
        tit_lab = paste0("<span style='font-size:35'><b>", paste0('Factor ', pc, '  '), "</b><span style='font-size:30'>",
                         (avg_expl_var %>% filter(Total_Factors == pc))$LAB)
        tit_lab <- rich_text_grob(tit_lab,
                                  x = unit(-44, "lines"),#-0.62 * uniqueN(res_PCA_loadings_adj$variable) * uniqueN(res_PCA_loadings_adj$year), "lines"), # for 17*8 should be 95
                                  y = unit(72, "lines")) # for 2 PC should be 33    83 -25 * PC_to_compare
        
        p = ggplot(res_DFM_loadings %>%
                     filter(data == df_type) %>%
                     filter(method == recov_met) %>%
                     filter(Total_Factors == pc) %>%
                     filter(DFM %in% dfm_to_plot$DFM) %>%
                     left_join(dfm_to_plot, by = "DFM") %>%
                     select(-DFM) %>%
                     rename(DFM = label) %>%
                     mutate(country = as.factor(country),
                            loading = ifelse(loading != 0, sign(loading) * log10(abs(loading)), 0),
                            variable = factor(variable, levels=paste0('var', 1:uniqueN(variable)))),
                   aes(x=loading)) + 
          geom_density(color="darkblue", fill="lightblue") +
          ggtitle(paste0('Explained Variance: ', paste0((avg_expl_var %>% filter(Total_Factors == pc))$LAB, collapse = " | "), '\n')) +
          # facet_grid(variable ~ DFM, scales = 'fixed', switch="y")
          facet_wrap(~ variable, ncol = 6, dir = 'h', scales = 'fixed', strip.position = 'top') +
          theme(axis.title.y=element_blank(),
                axis.title.x=element_blank(),
                axis.text.x = element_text(size = 24, vjust=0.3),
                axis.text.y=element_text(size=20),
                legend.text=element_text(size=28),
                legend.title=element_text(size=28),
                legend.key.size = unit(1.3, "cm"),
                strip.text = element_text(size = 24),
                plot.title = element_text(size=38),
                plot.margin = ggplot2::margin(2, 0, 0.7, 0.3, "cm"),
                panel.background = element_rect(fill = "white", colour = "black"),
                panel.grid.major.y = element_line(colour = "black", linetype = 'dashed', size = 0.2))
        
        p = ggplotGrob(p)
        # p$grobs[[16]] <- tit_lab
        p = gtable_add_padding(p, unit(c(0,0,0.35,0), "cm")) # t,r,b,l
        row_list[[pc]] = p
      } # pc
      g = do.call(rbind, c(row_list, size="last"))
      # g = gtable_add_padding(g, unit(2, "cm")) # t,r,b,l
      png('./Paper/Latex_Table_Figure/DFM_loadings.png', width = 22, height = 17, units = 'in', res=300)
      grid.draw(g)
      dev.off()
      
    }
  }
  
  # country map of index
  {
    
    df_type = 'Restricted'
    recov_met = 'TENSOR_BF'
    index_to_plot = data.frame(pc_met = c('RobPCA', 'DFM_multivar_stacked'), title = c('Robust PCA', 'DFM'), stringsAsFactors = F)
    factor_to_plot = 1
    
    res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
    
    for (i in 1:nrow(index_to_plot)){
      
      index_data = res_ALL_scores %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        filter(algo == index_to_plot$pc_met[i]) %>%
        filter(factor == factor_to_plot)
      
      for (yr in c('Year Average', unique(index_data$year))){
        
        if (yr == 'Year Average'){
          plot_data = index_data %>%
            group_by(country) %>%
            summarize(index = mean(index), .groups = "drop")
        } else {
          plot_data = index_data %>%
            filter(year == yr) %>%
            select(country, index)
        }
        
        var_mapping = data.frame(dataset = c("Antigua and Barbuda", "Bahamas, The", "Brunei Darussalam", "Cabo Verde", "Congo, Dem. Rep.", "Congo, Rep.",
                                             "Cote d'Ivoire", "Egypt, Arab Rep.", "Gambia, The", "Iran, Islamic Rep.", "Korea, Rep.", "Kyrgyz Republic",
                                             "Lao PDR", "North Macedonia", "Russian Federation", "Slovak Republic", "St. Lucia", "St. Vincent and the Grenadines",
                                             "Syrian Arab Republic", "Trinidad and Tobago", "United Kingdom", "United States", "Venezuela, RB", "Yemen, Rep."),
                                 map = c('Antigua', 'Bahamas', 'Brunei', 'Cape Verde', 'Democratic Republic of the Congo', 'Republic of Congo', 'Ivory Coast',
                                         'Egypt', 'Gambia', 'Iran', 'South Korea', 'Kyrgyzstan', 'Laos', 'Macedonia', 'Russia', 'Slovakia', 'Saint Lucia',
                                         'Saint Vincent', 'Syria', 'Trinidad', 'UK', 'USA', 'Venezuela', 'Yemen'), stringsAsFactors = F)
        
        all_country <- map_data("world")
        data_map = all_country %>%
          left_join(
            plot_data %>%
              left_join(var_mapping, by = c('country' = 'dataset')) %>%
              mutate(map = ifelse(is.na(map), country, map)),
            by = c('region' = 'map'))
        
        p = ggplot() + geom_polygon(data = data_map, aes(x=long, y=lat, group = group, fill=index),colour="black") + 
          scale_fill_gradientn(colours = rev(colorspace::diverge_hcl(uniqueN(data_map$index))), guide="colorbar", na.value="white",
                               limits = range(index_data$index),
                               breaks =  seq(min(index_data$index), max(index_data$index), diff(range(index_data$index)) / 7),
                               labels = round(seq(min(index_data$index), max(index_data$index), diff(range(index_data$index)) / 7), 2)) +
          theme_bw()  + labs(fill = "ESR", title = '', x="", y="") +
          scale_y_continuous(breaks=c()) + scale_x_continuous(breaks=c()) +
          ggtitle(paste0(index_to_plot$title[i], ' index - ', yr)) +
          theme(plot.title = element_text(size=30),
                legend.position="right",
                legend.text = element_text(size=14),
                legend.title = element_text(size=20),
                legend.key.height = unit(2, "cm")) +
          theme(plot.title=element_text(size=30, vjust=1.25))
        
        png(paste0('./Paper/Latex_Table_Figure/Index_country_world_map_', index_to_plot$pc_met[i], '_', gsub(" ", "", yr), '.png'),
            width = 13, height = 8, units = 'in', res=300)
        plot(p)
        dev.off()
        # png(paste0('./Results/1_Index_country_world_map_', df_type, '_', recov_met, '_', index_to_plot$pc_met[i], '_', gsub(" ", "", yr), '.png'),
        #     width = 13, height = 10, units = 'in', res=300)
        # plot(p)
        # dev.off()
      } # yr
    } # i
    
  }
  
  # regression task performance tables
  {
    df_type = 'Restricted'
    recov_met = 'TENSOR_BF'
    single_var = 'COG_wdi_unempilo'
    
    perf = read.csv(paste0("./Results/3b_prediction_table_", df_type, "_", recov_met, ".csv"), sep = ";", header = T, dec = ".", stringsAsFactors = F)
    
    # single variable table
    single_table = perf %>%
      filter(target_var == single_var) %>%
      mutate(lab = paste0('$', round(Continuous_PERF, 3), ' (', round(Original_PERF, 3), ')$')) %>%
      select(algo_main, fit_method, lab) %>%
      spread(fit_method, lab) %>%
      rename(Algorithm = algo_main)
    single_table = data.frame(t(colnames(single_table))) %>% `colnames<-`(colnames(single_table)) %>%
      bind_rows(single_table) %>%
      setNames(c('', 'RMSE index (RMSE original)', ''))
    write.table(single_table, './Paper/Latex_Table_Figure/05_single_macro_perf.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".", quote = F)
    
    # remaing variables table
    max_var_col = 3
    multiple_table = perf %>%
      filter(target_var != single_var) %>%
      mutate(lab = paste0('$', round(Continuous_PERF, 3), ' (', round(Original_PERF, 3), ')$')) %>%
      select(Description, algo_main, fit_method, lab) %>%
      setDT() %>%
      dcast(algo_main ~ Description + fit_method, value.var = "lab", sep = "___") %>%
      rename(Algorithm = algo_main)
    col_df = data.frame(lab=colnames(multiple_table)) %>%
      separate(lab, c('first', 'second'), sep= '___', remove = T, fill = 'left') %>%
      replace(is.na(.), 'Target variable')
    multiple_table = data.frame(t(col_df)) %>% `colnames<-`(colnames(multiple_table)) %>%
      bind_rows(multiple_table) %>%
      setNames(1:ncol(.))
      
    final_table = c()
    for (i in 1:ceiling((ncol(multiple_table)-1)/(2*max_var_col))){
      tt = multiple_table[, c(1, c(((i - 1) * 2 * max_var_col+2):min(ncol(multiple_table), (i*2*max_var_col + 1))))] %>%
        setNames(1:ncol(.))
      final_table = final_table %>%
        bind_rows(tt,
                  data.frame(t(rep('',7))) %>% setNames(1:(2 * max_var_col+1)))
    }
    final_table = final_table %>%
      replace(is.na(.), '') %>%
      setNames(c('', 'RMSE index (RMSE original)', rep('', ncol(.)-2)))
    write.table(final_table[-nrow(final_table),], './Paper/Latex_Table_Figure/06_multiple_macro_perf.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".", quote = F)
  }
  
  # index evolution
  {
    df_type = 'Restricted'
    recov_met = 'TENSOR_BF'
    algo_to_plot = c('RobPCA', 'DFM_multivar_stacked')
    res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
    
    # index for sample countries
    sample_list = c('Australia', 'Brazil', 'Yemen, Rep.', 'Algeria', 'Romania', 'Spain')
    country_epidemic = data.frame(country = c('Australia', 'Brazil', 'Yemen, Rep.', 'Algeria', 'Romania', 'Spain'),
                                start = c(17, 17, 16, 17, 17, 17),
                                end = c(19, 18, 18, 18, 18, 18), stringsAsFactors = F)
    
    comp_data = res_ALL_scores %>%
      filter(data == df_type) %>%
      filter(method == recov_met)  %>%
      filter(algo %in% algo_to_plot) %>%
      mutate(factor = paste0('Index ', factor),
             color_id = paste0(family, '_', algo)) %>%
      mutate(algo = replace(algo, algo == 'DFM_multivar', 'DFM')) %>%
      mutate(algo = replace(algo, algo == 'DFM_multivar_stacked', 'DFM')) %>%
      mutate(algo = replace(algo, algo == 'RobPCA', 'Robust PCA'))
    
    
    y_lim = comp_data %>% filter(country %in% sample_list) %>% pull(index) %>% range(na.rm = T)
    for (samp in sample_list){
      png(paste0('./Paper/Latex_Table_Figure/Evolution_', gsub(' ', '_', gsub(',', '', gsub('\\.', '', samp))), '.png'), width = 10, height = 6, units = 'in', res=200) 
      d = comp_data %>%
        filter(country == samp) %>%
        mutate(year = as.numeric(substr(year, 3, 4))) %>%
        left_join(country_epidemic, by = "country") %>%
        mutate(m = min(year),
               M = max(year)) %>%
        rowwise() %>%
        mutate(start = max(c(m, start)),
               end = min(c(M, end)))
      
      plot(ggplot(d,
                  aes(y=index, x=year)) + 
             geom_rect(
               aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
               fill = "red", alpha = 0.05) +
             geom_line(aes(colour = algo), size = 3.5) +
             ylim(y_lim) +
             scale_colour_manual("", values=c('blue3', 'brown')) +
             # scale_linetype_manual("", values=c('solid', 'twodash', 'dotted', 'dotdash')) +
             # facet_wrap(country ~ factor, dir = 'v', scales = 'free', strip.position = 'top', ncol = length(coup)) +
             scale_x_continuous(labels = paste0('\'', as.character(sort(unique(d$year)))), breaks = unique(d$year)) +
             # ggtitle(samp) +
             theme(panel.background = element_rect(fill = "white", colour = "black"),
                   panel.grid = element_line(colour = "black", linetype = 'dashed', size = 0.4),
                   panel.grid.minor = element_blank(),
                   strip.text = element_text(size = 22),
                   axis.text.x = element_text(size = 25),
                   axis.text.y = element_text(size = 20),
                   plot.title = element_text(size = 34),
                   legend.position = "top",
                   legend.text=element_text(size=28),
                   legend.key.width = unit(1.7,"cm")) +
             labs(x = '', y = ''))
      dev.off()
    }
    
    country_list = unique(comp_data$country)
    step = 42
    y_lim = comp_data %>% pull(index) %>% range(na.rm = T)
    for (country_split in 1:ceiling(length(country_list)/step)){
      country_range = ((country_split - 1) * step + 1):min(c(country_split * step, length(country_list)))

      png(paste0('./Paper/Latex_Table_Figure/Evolution_full_', country_split, '.png'), width = 15, height = 20, units = 'in', res=300)
      d = comp_data %>%
        filter(country %in% country_list[country_range]) %>%
        mutate(year = as.numeric(substr(year, 3, 4)))
      plot(ggplot(d, aes(fill=algo, y=index, x=year)) + 
             geom_line(aes(colour = algo), size = 1.7) +
             scale_colour_manual("", values=c('blue3', 'brown')) +
             scale_linetype_manual("", values=c('solid', 'twodash', 'dotted', 'dotdash')) +
             ylim(y_lim) +
             facet_wrap(country ~ ., dir = 'v', scales = 'free', strip.position = 'top', ncol = 5) +
             scale_x_continuous(labels = paste0('\'', as.character(sort(unique(d$year)))), breaks = unique(d$year)) +
             theme(panel.background = element_rect(fill = "white", colour = "black"),
                   panel.grid = element_line(colour = "black", linetype = 'dashed', size = 0.4),
                   panel.grid.minor = element_blank(),
                   strip.text = element_text(size = 16),
                   axis.text.x = element_text(size = 11),
                   legend.text=element_text(size=30),
                   legend.direction="horizontal",legend.position="top") +
             labs(x = '', y = ''))
      dev.off()
      
      
    }
  }
  
  # stationarity test for final index and top/bottom value for index
  {
    # IPS H1: stationarity, tested for both "individual intercepts" and "individual intercepts and trends" option for Augmented-Dickey-Fuller
    
    df_type = 'Restricted'
    recov_met = 'TENSOR_BF'
    res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
    max_tob_bottom = 10
    
    data_df = res_ALL_scores %>%
      filter(data == df_type) %>%
      filter(method == recov_met)
    
    top_bottom_val = c()
    IPS_test = c()
    for (algo_type in unique(data_df$algo)){
      
      dd = data_df %>% filter(algo == algo_type) %>% select(country, year, index)
      pdata = plm::pdata.frame(dd, index=c("country","year"), drop.index=TRUE, row.names=TRUE)
      IPS_test = IPS_test %>%
        bind_rows(data.frame(data = df_type, method = recov_met, algo = algo_type,
                             IPS_intercept_pval = suppressWarnings(plm::purtest(pdata, test = "ips", exo = "intercept", lags = "AIC", pmax = 5)$statistic$p.value),
                             IPS_trend_pval = suppressWarnings(plm::purtest(pdata, test = "ips", exo = "trend", lags = "AIC", pmax = 5)$statistic$p.value),
                             IPS_H1 = 'stationarity', stringsAsFactors = F))
      
      top_bottom_val = top_bottom_val %>% 
        bind_cols(dd %>%
                    group_by(country) %>%
                    summarize(index = index[which.max(abs(index))], .groups = "drop") %>%
                    # summarize(index = mean(index), .groups = "drop") %>%
                    arrange(index) %>%
                    mutate(label = paste0(country, ' ', round(index, 3))) %>%
                    filter(row_number() <= max_tob_bottom | row_number() >= (n() - max_tob_bottom + 1)) %>%
                    select(label) %>%
                    rename(!!sym(algo_type) := label))

    }
    write.table(IPS_test, './Paper/Latex_Table_Figure/07_index_stationarity_test.csv', sep = ";", col.names = T, row.names = F, append = F)
    write.table(top_bottom_val, './Paper/Latex_Table_Figure/08_index_top_bottom.csv', sep = ";", col.names = T, row.names = F, append = F)
  }
  
  # table for correlation with disease indices
  {
    country_to_keep = c('Indonesia', 'Angola', 'Dominican Republic', 'Nigeria', 'Netherlands', 'Pakistan', 'France', 'Argentina', 'Brazil')
    
    corr_data = read.csv("./Results/4_correlation_power_summary.csv", sep=";", stringsAsFactors=FALSE) %>%
      filter(data == 'RestrictedDifference') %>%
      filter(method == 'SOFT_IMPUTE') %>%
      filter(algo == 'DFM_multivar_stacked') %>%
      filter(country %in% country_to_keep) %>%
      select(country, reference_index, spearman_corr) %>%
      spread(reference_index, -country) %>%
      mutate_if(is.numeric, ~as.character(abs(round(., 2)))) %>%
      mutate_at(vars(-country), ~ifelse(as.numeric(.) >= 0.5, paste0("bbbb", .), .)) %>%
      mutate_at(vars(-country), ~ifelse(!is.na(.), paste0(., sample(c("*", "**", "***"), 1)), .)) %>%
      replace(is.na(.), "") %>%
      mutate_all(~ifelse(grepl("bbbb", .), paste0("\\", gsub("bbbb", "textbf{", .), "}"), .))
    
    write.table(corr_data, './Paper/Latex_Table_Figure/09_index_correlation.csv', sep = ";", col.names = T, row.names = F, append = F)
  }
}
