# evaluate correlation and partial correlation
evaluate_correlation = function(df, corr_method = 'pearson'){
  # df: data.frame of variables
  
  corr = rcorr(df %>% as.matrix(), type = corr_method)
  corr_v = corr$r
  corr_v[lower.tri(corr_v, diag = T)] = NA
  corr_v = reshape2::melt(corr_v) %>% filter(!is.na(value)) %>% setNames(c('Var1', 'Var2', 'Corr')) %>% mutate_if(is.factor, as.character)
  corr_p = round(corr$P, 8)
  corr_p[lower.tri(corr_p, diag = T)] = NA
  corr_p = reshape2::melt(corr_p) %>% filter(!is.na(value)) %>% setNames(c('Var1', 'Var2', 'Corr_pVal')) %>% mutate_if(is.factor, as.character)
  
  # pcorr = suppressWarnings(pcor(df, method = corr_method))
  # pcorr_v = pcorr$estimate
  # rownames(pcorr_v) = colnames(pcorr_v) = colnames(df)
  # pcorr_v[lower.tri(pcorr_v, diag = T)] = NA
  # pcorr_v = reshape2::melt(pcorr_v) %>% filter(!is.na(value)) %>% setNames(c('Var1', 'Var2', 'PartCorr')) %>% mutate_if(is.factor, as.character)
  # pcorr_p = round(pcorr$p.value, 8)
  # rownames(pcorr_p) = colnames(pcorr_p) = colnames(df)
  # pcorr_p[lower.tri(pcorr_p, diag = T)] = NA
  # pcorr_p = reshape2::melt(pcorr_p) %>% filter(!is.na(value)) %>% setNames(c('Var1', 'Var2', 'PartCorr_pVal')) %>% mutate_if(is.factor, as.character)
  
  correlation_list = corr_v %>%
    left_join(corr_p, by = c("Var1", "Var2")) %>%
    # left_join(pcorr_v, by = c("Var1", "Var2")) %>%
    # left_join(pcorr_p, by = c("Var1", "Var2")) %>%
    # mutate(abs = abs(PartCorr)) %>%
    mutate(abs = abs(Corr)) %>%
    arrange(desc(abs))
  
  return(correlation_list)
}

# clever standardization
clever_scale = function(input_df, exclude_var){
  # if variable is constant or included in 'exclude_var', no transformation
  
  output_matrix = c()
  
  for (col in colnames(input_df)){
    val = input_df %>% pull(col)
    if (col %in% exclude_var){
      std_val = val
    } else if (uniqueN(val) == 1){
      std_val = scale(val, center = T, scale = F) %>% as.numeric
    }else {
      std_val = scale(val, center = T, scale = T) %>% as.numeric
    }
    output_matrix = output_matrix %>%
      cbind(std_val)
  }
  output_matrix = output_matrix %>%
    as.data.frame() %>%
    setNames(colnames(input_df)) %>%
    `rownames<-`(rownames(input_df)) %>%
    as.matrix()
  
  return(output_matrix)
}

# Statistical analysis: numerical, data and character columns
basicStatistics = function(data){
  
  data=as.data.frame(data, stringsAsFactors = F)
  
  # Get numerical columns
  nums <- names(which(sapply(data, is.numeric)))
  if (length(nums)>0){
    StatNum = c()
    for (col in nums){
      dd = data[, col]
      perc = quantile(dd, c(0.01, 0.05, 0.95, 0.99), na.rm = T)
      StatNum = StatNum %>% bind_rows(
        data.frame(UNIQUE_VALS = uniqueN(dd, na.rm = T),
                   NAs = sum(is.na(dd)),
                   Min = min(dd, na.rm = T),
                   Max = max(dd, na.rm = T),
                   Mean = mean(dd, na.rm = T),
                   StDev = sd(dd, na.rm = T),
                   Median = median(dd, na.rm = T),
                   Perc_1 = perc[1],
                   Perc_5 = perc[2],
                   Perc_95 = perc[3],
                   Perc_99 = perc[4],
                   Sum = sum(dd, na.rm = T), stringsAsFactors = F)
      )
    }
    StatNum = StatNum %>%
      mutate(NAs = paste0(NAs,' (',signif(NAs/nrow(data)*100,digits=2),'%)')) %>%
      mutate_all(as.character) %>%
      `rownames<-`(nums)
  } else {StatNum=data.frame()}
  
  # Get dates columns
  dates <- names(which(sapply(data, lubridate::is.Date)))
  if (length(dates)>0){
    StatDat = c()
    for (col in dates){
      dd = data[, col]
      StatDat = StatDat %>% bind_rows(
        data.frame(UNIQUE_VALS = uniqueN(dd, na.rm = T),
                   NAs = sum(is.na(dd)),
                   Min = min(dd, na.rm = T),
                   Max = max(dd, na.rm = T),
                   Median = median(dd, na.rm = T), stringsAsFactors = F) %>%
          bind_cols(
            data.frame(val = sort(dd), stringsAsFactors = F) %>%
              filter(!is.na(val)) %>%
              mutate(POS = 1:n()) %>%
              filter(POS %in% round(n() * c(0.01, 0.05, 0.95, 0.99))) %>%
              mutate(POS = paste0('Perc_', c(1, 5, 95, 99))) %>%
              mutate(group = 1) %>%
              spread(POS, val) %>%
              select(-group)
          )
      )
    }
    StatDat = StatDat %>%
      mutate(NAs = paste0(NAs,' (',signif(NAs/nrow(data)*100,digits=2),'%)')) %>%
      mutate_all(as.character)
    rownames(StatDat) = dates
  } else {StatDat=data.frame()}
  
  # Get characters columns
  chars <- names(which(sapply(data, is.character)))
  if (length(chars)>0){
    StatChar = c()
    for (i in chars){
      cc=data[,i];cc=cc[!is.na(cc)]
      if (length(unique(cc)) <= 50){line = paste(unique(cc), collapse = "|")
      } else {line = paste("> 50 unique values",paste(unique(cc)[1:20], collapse = "|"))}
      StatChar = rbind(StatChar,c(paste0(nrow(data)-length(cc),' (',signif((nrow(data)-length(cc))/nrow(data)*100,digits=2),'%)'),
                                  paste0(sum(cc == ''),' (',signif(sum(cc == '')/nrow(data)*100,digits=2),'%)'),
                                  length(unique(cc)),line))
    }
    if (length(StatChar)>0){StatChar = data.frame(StatChar, stringsAsFactors = F) %>% setNames(c('NAs', 'BLANKs', 'UNIQUE_VALS','VALUES'));rownames(StatChar)=chars}
  } else {StatChar=data.frame()}
  
  Stat=bind_rows(StatNum,StatDat,StatChar)
  
  Stat=cbind(VARIABLE=c(rownames(StatNum),rownames(StatDat),rownames(StatChar)),
             TYPE=c(rep('NUMERIC',length(nums)),rep('DATE',length(dates)),rep('CHARACTER',length(chars))),
             NUM_OSS=nrow(data),Stat)
  col_first = intersect(c('VARIABLE', 'TYPE', 'NUM_OSS', 'UNIQUE_VALS', 'NAs', 'BLANKs'), colnames(Stat))
  Stat = Stat %>% select(all_of(col_first), everything())
  
  final=data.frame(VARIABLE=colnames(data), stringsAsFactors = F)
  final = final %>% left_join(Stat, by = "VARIABLE")
  final = as.matrix(final)
  final[is.na(final)]=''
  final = as.data.frame(final, stringsAsFactors = F)
  
  return(final)
}

# variable stats
variable_stats = function(df){
  
  summ_variable = c()
  for (i in colnames(df)){
    # tt = df %>% select_(.dots = i)
    # eval(parse(text=paste0('tt<-df$', i)))
    tt = df[, i]
    summ_variable = rbind(summ_variable, c(i,
                                           sum(is.na(tt)),
                                           round(sum(is.na(tt)) / nrow(df) * 100, digit = 1),
                                           min(tt, na.rm = T),
                                           max(tt, na.rm = T)))
  }
  colnames(summ_variable) = c('VARIABLE', 'NA', 'NA%', 'MIN', 'MAX')
  summ_variable = as.data.frame(summ_variable, stringsAsFactors = F) %>% mutate(`NA` = as.numeric(`NA`), `NA%` = as.numeric(`NA%`)) %>%
    arrange(`NA%`)
  
  return(summ_variable)
}

# country missing stats
country_stats = function(df, col_filter, year_filter){
  
  summ_country = df %>%
    select_(.dots = col_filter) %>%
    filter(year %in% year_filter) %>%
    gather('variable', 'val', -c(year, country)) %>%
    group_by(country) %>%
    summarize(YEAR_MIN = min(year),
              YEAR_MAX = max(year),
              YEAR_COUNT = uniqueN(year),
              `NA` = sum(is.na(val)),
              `NA%` = round(`NA` / n() * 100, 1)) %>%
    ungroup() %>%
    arrange(`NA%`)
  
  # summ_country = df %>%
  #   select_(.dots = col_filter) %>%
  #   filter(year %in% year_filter) %>%
  #   group_by(country) %>%
  #   summarize(YEAR_MIN = min(year),
  #             YEAR_MAX = max(year),
  #             YEAR_COUNT = n()) %>%
  #   ungroup() %>%
  #   left_join(
  #     df %>% 
  #       select_(.dots = setdiff(col_filter, c('year', 'c_code'))) %>%
  #       group_by(country) %>% summarise_all(funs(n())) %>%
  #       mutate(DEN = rowSums(.[-1])) %>%
  #       select(country, DEN) %>%
  #       left_join(
  #         df %>% 
  #           select_(.dots = setdiff(col_filter, c('year', 'c_code'))) %>%
  #           group_by(country) %>% summarise_all(funs(sum(is.na(.)))) %>%
  #           mutate(`NA` = rowSums(.[-1])) %>%
  #           select(country, `NA`),
  #         by = "country"
  #       ) %>%
  #       mutate(`NA%` = round(`NA` / DEN * 100, digits = 1)) %>%
  #       select(-DEN),
  #     by = "country"
  #   ) %>%
  #   arrange(`NA%`)
  
  return(summ_country)
}

# function to be called with try catch - NIPALS
t_nipals = function(slice){
  nip = nipals(slice,
               center = TRUE,
               scale = TRUE,
               fitted = TRUE,
               force.na = TRUE,
               maxiter = 10000)
}
try_nipals = function(x){
  r <- tryCatch({
    list(value = t_nipals(x), no_error = "No error")
  }, warning = function(e) {
    warn_text <- paste0("WARNING: ", e)
    return(list(value = (suppressWarnings(t_nipals(x))), warn_text = warn_text))
  }, error = function(e) {
    error_text <- paste0("ERROR: ", e)
    return(list(value = NA, error_text = error_text))
  }, finally = {
  }, quiet = TRUE)
  return(r)
}

# function to be called with try catch - Soft_Impute
t_soft_imp = function(slice){
  if (sum(rowSums(is.na(slice)) == ncol(slice)) + sum(colSums(is.na(slice)) == nrow(slice)) == 0){
    slice_t = biScale(as.matrix(slice), maxit = 200)  # scale input, works only with matrix without full NA rows or cols
  } else {slice_t = as.matrix(slice)}
  lambda_sel = lambda0(slice_t,lambda = 0, maxit = 200, trace.it = F, thresh = 1e-05)
  sf = softImpute(slice_t,
                  rank.max = min(dim(slice)) - 1,
                  lambda = lambda_sel,
                  type = "als",
                  maxit = 10000,
                  trace.it = F,
                  warm.start = NULL,
                  final.svd = TRUE)
  return(softImpute::complete(as.matrix(slice),sf) %>% as.data.frame())
}
try_soft_imp = function(x){
  r <- tryCatch({
    list(value = t_soft_imp(x), no_error = "No error")
  }, warning = function(e) {
    warn_text <- paste0("WARNING: ", e)
    return(list(value = (suppressWarnings(t_soft_imp(x))), warn_text = warn_text))
  }, error = function(e) {
    error_text <- paste0("ERROR: ", e)
    return(list(value = NA, error_text = error_text))
  }, finally = {
  }, quiet = TRUE)
  return(r)
}

# function to be called with try catch - TensorBF
t_tensorBF = function(tens, K){
  opts <- getDefaultOpts()
  opts$iter.burnin <- 10000
  set.seed((10))
  res = tensorBF(tens,
                 K = K,
                 fiberCentering = 1,
                 slabScaling = 2,
                 noiseProp = c(0.5, 0.5),
                 opts = opts)
  return(predictTensorBF(tens, res))
}
try_tensorBF = function(x, K){
  r <- tryCatch({
    list(value = t_tensorBF(x, K), no_error = "No error")
  }, warning = function(e) {
    warn_text <- paste0("WARNING: ", e)
    return(list(value = (suppressWarnings(t_tensorBF(x, K))), warn_text = warn_text))
  }, error = function(e) {
    error_text <- paste0("ERROR: ", e)
    return(list(value = NA, error_text = error_text))
  }, finally = {
  }, quiet = TRUE)
  return(r)
}

# plot statistics on recovered missing
generate_plot_recover_stat = function(data, df_type, var, hide_axis, hide_leg, axis_title, max_lim){
  p = ggplot(data, aes(x=variable, y=var, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ylab(axis_title) +
    coord_flip(ylim = c(0, max_lim)) +
    ggtitle(df_type)
  if (hide_axis){
    p = p + theme(axis.text.y=element_blank())
  }
  if (hide_leg){
    p = p + theme(legend.position = "none")
  }
  return(p)
}

plot_recover_stat = function(stats_recover, var, log_tran){
  
  data = stats_recover %>% select(data, variable, method, var) %>% rename_('var' = var)
  if (log_tran){
    data$var = log(data$var)
    axis_title = paste0('Log(', var, ')')
  } else {
    axis_title = var
  }
  max_lim = max(data$var)
  c = 1; p_list = list(); c_max = uniqueN(data$data)
  for (df_type in (unique(data$data))){
    p_list[[c]] = ggplotGrob(generate_plot_recover_stat(data %>% filter(data == df_type), df_type, var,
                                                        hide_axis = ifelse(c == 1, F, T),
                                                        hide_leg = ifelse(c == c_max, F, T),
                                                        axis_title, max_lim))
    c = c + 1
  }
  p = do.call(cbind, c(p_list, size="first"))
  
  return(p)
}

# PCA output
pca_fun = function(slice, k_fold, cv, method_title){
  pca = prcomp(slice,
               center = F,
               scale = F)
  
  # find optimal number of PC
  PC_opt = c(PC_CV = ncol(slice))
  if(cv){
    PC_opt = cv_pca(slice, func = pca_fun, k_fold)
  }
  pca$PC_opt = PC_opt
  
  load = pca$rotation
  summ = summary(pca)
  scree_plot = fviz_eig_mod(pca, addlabels = T, pc_to_highlight = min(PC_opt))
  scree_plot$labels$title = paste0(method_title, ' - ', scree_plot$labels$title)
  load_plot = fviz_pca_var(pca,
                           col.var = "contrib", # Color by contributions to the PC
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                           repel = TRUE     # Avoid text overlapping
  )
  load_plot$labels$title = paste0(method_title, ' - ', load_plot$labels$title)
  bi_plot = fviz_pca_biplot(pca, repel = TRUE,
                            col.var = "#2E9FDF", # Variables color
                            col.ind = "#696969"  # Individuals color
  )
  bi_plot$labels$title = paste0(method_title, ' - ', bi_plot$labels$title)
  load_table = data.frame(variable = rep(rownames(load), ncol(load)),
                          PC = rep(c(1:ncol(load)), rep(nrow(load), ncol(load))),
                          loading = as.numeric(load), stringsAsFactors = F)
  importance_table = data.frame(PC = colnames(summ$importance), stringsAsFactors = F) %>%
    cbind(as.data.frame(t(summ$importance)[, -1], stringsAsFactors = F))
  return(list(pca = pca,
              load_table = load_table,
              importance_table = importance_table,
              scree_plot = scree_plot,
              load_plot = load_plot,
              bi_plot = bi_plot))
}

# Robust Sparse PCA output
spca_fun = function(slice, k_fold, cv, method_title){        
  ss = robspca(slice, k = NULL,
               alpha = 5e-04,  # Sparsity controlling parameter. Higher values lead to sparser components
               beta = 1e-04,   # Amount of ridge shrinkage to apply in order to improve conditioning
               gamma = 100,     # Sparsity controlling parameter for the error matrix S. Smaller values lead to a larger amount of noise removal
               center = F, scale = F,
               max_iter = 1000,
               verbose = F)
  
  # check reconstruction error
  conv_err = c()
  rr = t(t(ss$scores %*% t(ss$transform)))
  avg_abs = mean(abs(as.matrix(slice)))
  max_var = max(abs(rr - as.matrix(slice)))
  toll = 0.1 # percentage
  if(max_var / avg_abs > toll){conv_err = sum(abs(rr - as.matrix(slice)) / avg_abs > toll)}
  
  # create pca for class structure
  pp = prcomp(slice,
              center = F,
              scale = F)
  pp$sdev = ss$sdev
  ll = ss$loadings
  rownames(ll) = rownames(pp$rotation)
  colnames(ll) = colnames(pp$rotation)
  pp$rotation = ll
  pp$center = ss$center
  pp$scale = ss$scale
  pp$x = scale(slice) %*% ll
  
  # find optimal number of PC
  PC_opt = c(PC_CV = ncol(slice))
  if(cv){
    PC_opt = cv_pca(slice, func = spca_fun, k_fold)
  }
  pp$PC_opt = PC_opt
  
  load = pp$rotation
  summ = summary(pp)
  scree_plot = fviz_eig_mod(pp, addlabels = T, pc_to_highlight = min(PC_opt))
  scree_plot$labels$title = paste0(method_title, ' - ', scree_plot$labels$title)
  load_plot = fviz_pca_var(pp,
                           col.var = "contrib", # Color by contributions to the PC
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                           repel = TRUE     # Avoid text overlapping
  )
  load_plot$labels$title = paste0(method_title, ' - ', load_plot$labels$title)
  bi_plot = fviz_pca_biplot(pp, repel = TRUE,
                            col.var = "#2E9FDF", # Variables color
                            col.ind = "#696969"  # Individuals color
  )
  bi_plot$labels$title = paste0(method_title, ' - ', bi_plot$labels$title)
  load_table = data.frame(variable = rep(rownames(load), ncol(load)),
                          PC = rep(c(1:ncol(load)), rep(nrow(load), ncol(load))),
                          loading = as.numeric(load), stringsAsFactors = F)
  importance_table = data.frame(PC = colnames(summ$importance), stringsAsFactors = F) %>%
    cbind(as.data.frame(t(summ$importance)[, -1], stringsAsFactors = F))
  sparseness_count = load_table %>% group_by(PC) %>% summarize(ZERO_ELEM = sum(loading == 0), .groups = 'drop')
  
  return(list(pca = pp,
              load_table = load_table,
              importance_table = importance_table,
              sparseness_count = sparseness_count,
              scree_plot = scree_plot,
              load_plot = load_plot,
              bi_plot = bi_plot,
              conv_err = conv_err,
              sparse = ss$sparse))
}

# PC number selection
cv_pca = function(slice, func, k_fold){
  
  # Wold Cross Validation - https://www.researchgate.net/publication/5638191_Cross-validation_of_component_models_A_critical_look_at_current_methods
  
  n = nrow(slice)
  max_pc = min(dim(slice))
  set.seed(10)
  cv_ind = caret::createFolds(c(1:n), k = k_fold, list = TRUE, returnTrain = FALSE)
  max_pc_fold = min(c(n - max(unlist(lapply(cv_ind,length))), max_pc))
  cv_err = c()
  for (fold in 1:length(cv_ind)){
    cal_set = slice[setdiff(1:n, cv_ind[[fold]]), ] %>% as.matrix()
    val_set = slice[cv_ind[[fold]], ] %>% as.matrix()
    cv_est = func(slice = cal_set, k_fold = NULL, cv = F, method_title = '')
    
    for (pc in 1:(max_pc_fold - 1)){
      load_pc = cv_est$pca$rotation[, 1:pc]
      score_pc = val_set %*% load_pc
      val_set_approx = score_pc %*% t(load_pc)
      cv_err = cv_err %>% rbind(c(FOLD = fold, PC = pc, ERR = sum((val_set - val_set_approx) ^ 2)))
      
    }
  }
  tot_var = sum(slice ^ 2)
  cv_err_avg = cv_err %>% as.data.frame() %>%
    group_by(PC) %>%
    summarise(R = sum(ERR) / tot_var) %>%
    filter(R < 1) %>%
    arrange(desc(PC))
  pc_opt_wold = cv_err_avg$PC[1]
  
  
  # Malinowski error (REV - Reduced EigenValue)
  
  F_signif = 0.1
  
  eval_pca = func(slice, k_fold = NULL, cv = F, method_title = '')$pca
  eig = eval_pca$sdev^2
  # rev_l = c()
  # for (pc in 1:(max_pc - 1)){
  #  den = 0
  #  for (m in (pc + 1):max_pc){
  #    den = den + eig[m] / ((max_pc - m + 1) * (n - m + 1))
  #  }
  #  F_ratio = (eig[pc] / ((max_pc - pc + 1) * (n - pc + 1))) / den
  #  F_quant = qf(F_signif, 1, n - pc, lower.tail = TRUE, log.p = FALSE)
  #  rev_l = c(rev_l, ifelse(F_ratio < F_quant, 0, 1))
  # }
  REV = rev_l = c()
  for (pc in 1:(max_pc - 1)){
    REV = c(REV, eig[pc] / ((max_pc - pc + 1) * (n - pc + 1)))
  }
  for (pc in 1:(max_pc - 2)){
    F_ratio = (n - pc) * REV[pc] / sum(REV[(pc + 1):(max_pc - 1)])
    F_quant = qf(F_signif, 1, n - pc, lower.tail = TRUE, log.p = FALSE)
    rev_l = c(rev_l, ifelse(F_ratio < F_quant, 0, 1))
  }
  pc_opt_REV =  suppressWarnings(min(which(rev_l == 0)) - 1)
  if (!is.finite(pc_opt_REV)){pc_opt_REV = max_pc - 1}
  
  
  # Q-stat
  
  Q_signif = 0.1
  
  Q_val = c()
  load = eval_pca$rotation
  for (pc in 1:(max_pc - 1)){
    r = (diag(max_pc) - load[1:max_pc, 1:pc] %*% t(load[1:max_pc, 1:pc])) %*% t(eval_pca$x)
    Q = diag(t(r) %*% r)
    T1 = sum(eig[(pc + 1):max_pc])
    T2 = sum(eig[(pc + 1):max_pc] ^ 2)
    T3 = sum(eig[(pc + 1):max_pc] ^ 3)
    h0 = 1 - 2 * T1 * T2 / (3 * T3 ^ 2)
    Q_lim = T1 * (1 + (h0 * qnorm(Q_signif) * sqrt(2 * T2)) / T1 + (T2 * h0 * (h0 - 1)) / T1 ^ 2) ^ (1 / h0)
    Q_val = c(Q_val, sum(Q > Q_lim))
  }
  pc_opt_Q = which.min(Q_val)
  
  return(c(Wold = pc_opt_wold,
           REV = pc_opt_REV,
           Q = pc_opt_Q))
}

# fviz_eig modified to show suggested number of PC
fviz_eig_mod = function (X, choice = c("variance", "eigenvalue"), geom = c("bar", 
                                                                           "line"), barfill = c("steelblue", 'red'), barcolor = "steelblue", 
                         linecolor = "black", ncp = 10, addlabels = FALSE, hjust = 0, 
                         main = NULL, xlab = NULL, ylab = NULL, ggtheme = theme_minimal(), pc_to_highlight,
                         ...){
  eig <- get_eigenvalue(X)
  eig <- eig[1:min(ncp, nrow(eig)), , drop = FALSE]
  choice <- choice[1]
  if (choice == "eigenvalue") {
    eig <- eig[, 1]
    text_labels <- round(eig, 0)
    if (is.null(ylab)) 
      ylab <- "Eigenvalue"
  }
  else if (choice == "variance") {
    eig <- eig[, 2]
    text_labels <- round(eig, 0)
  }
  else stop("Allowed values for the argument choice are : 'variance' or 'eigenvalue'")
  text_labels[5:length(text_labels)] = ""
  if (length(intersect(geom, c("bar", "line"))) == 0) 
    stop("The specified value(s) for the argument geom are not allowed ")
  df.eig <- data.frame(dim = factor(1:length(eig)), eig = eig)
  df.eig$hl = ifelse(df.eig$dim == pc_to_highlight, 1, 0)
  df.eig$hl = factor(df.eig$hl)
  extra_args <- list(...)
  bar_width <- extra_args$bar_width
  linetype <- extra_args$linetype
  if (is.null(linetype)) 
    linetype <- "solid"
  p <- ggplot(df.eig, aes(dim, eig, group = hl))
  if ("bar" %in% geom) 
    p <- p + geom_bar(stat = "identity", aes(fill=hl), 
                      color = barcolor, width = bar_width) +
    scale_fill_manual(values = barfill)
  if ("line" %in% geom) 
    p <- p + geom_line(color = linecolor, linetype = linetype) + 
    geom_point(shape = 19, color = linecolor)
  if (addlabels) 
    p <- p + geom_text(label = text_labels, vjust = -0.4, 
                       hjust = hjust)
  if (is.null(main)) 
    main <- "Scree plot"
  if (is.null(xlab)) 
    xlab <- "Dimensions"
  if (is.null(ylab)) 
    ylab <- "Percentage of explained variances"
  p <- p + labs(title = main, x = xlab, y = ylab) +
    theme(legend.position = "none")
  ggpubr::ggpar(p, ...)
  return(p)
}

# split variables into chunks for better visualisation
split_var = function(input, var_split_length){
  n_split = ceiling(length(input) / var_split_length)
  out_list = list()
  st = 1; end = var_split_length
  for (i in 1:n_split){
    out_list[[i]] = input[st:end]
    st = end + 1
    end = min(c(length(input), end + var_split_length))
  }
  
  return(out_list)
}

# split string in chunks and return a single string with \n between chuncks
split_string = function(string, chunk_length){
  n_split = ceiling(nchar(string) / chunk_length)
  out = c()
  for (i in 1:n_split){
    out = c(out, substr(string, 1, chunk_length))
    string = substr(string, chunk_length + 1, nchar(string))
  }
  
  return(paste0(out, collapse = '\n'))
}
split_string = Vectorize(split_string)

# sort loading according to one single variable for PCA
sort_loading_PCA = function(res_PCA_loadings, leading_var, leading_sign){
  res_PCA_loadings = res_PCA_loadings %>%
    mutate(LEAD_VAR_FLAG = ifelse(variable == leading_var, 1, 0)) %>%
    group_by(data, method, PCA, PC, year) %>%
    arrange(desc(LEAD_VAR_FLAG)) %>%
    mutate(LEAD_SIGN = ifelse(sign(loading) == 0, 1, sign(loading)) * LEAD_VAR_FLAG,
           SIGN_EXP = ifelse(leading_sign == 'p', 1, -1),
           SIGN_MULT = ifelse(SIGN_EXP == LEAD_SIGN[LEAD_VAR_FLAG == 1], 1, -1),
           loading = loading * SIGN_MULT) %>%
    select(-LEAD_VAR_FLAG, -LEAD_SIGN, -SIGN_EXP, -SIGN_MULT) %>%
    ungroup()
  
  return(res_PCA_loadings)
}

# sort loading according to one single variable for DFM
sort_loading_DFM = function(res_DFM_loadings, leading_var, leading_sign){
  res_DFM_loadings = res_DFM_loadings %>%
    mutate(LEAD_VAR_FLAG = ifelse(variable == leading_var, 1, 0)) %>%
    group_by(data, method, DFM, Total_Factors, Factor, country) %>%
    arrange(desc(LEAD_VAR_FLAG)) %>%
    mutate(LEAD_SIGN = ifelse(sign(loading) == 0, 1, sign(loading)) * LEAD_VAR_FLAG,
           SIGN_EXP = ifelse(leading_sign == 'p', 1, -1),
           SIGN_MULT = ifelse(SIGN_EXP == LEAD_SIGN[LEAD_VAR_FLAG == 1], 1, -1),
           loading = loading * SIGN_MULT) %>%
    select(-LEAD_VAR_FLAG, -LEAD_SIGN, -SIGN_EXP, -SIGN_MULT) %>%
    ungroup()
  
  return(res_DFM_loadings)
}

# evaluate random set (indexes of matrix) for missing recovering robustness test
get_random_set = function(set, n_subset, subset_length, seed_list, lab){
  
  left_out_perc = 0.05  # max percentage of elements to be not included in any subset
  
  # spanning all elements trying to cover all set
  
  out_list = list()
  expected_left_out = max(c(1, length(set) - n_subset * subset_length))
  for (i in 1:n_subset){
    # t_ind = unique(round(randtoolbox::sobol(subset_length * 2, dim = 1, init = TRUE, scrambling = 3, seed = seed_list[i], normal = FALSE) * subset_length * 2))
    # t_ind = setdiff(t_ind, 0)
    # set.seed(seed_list[i])
    # t_ind = sample(t_ind, subset_length)
    set.seed(seed_list[i])
    out_list[[i]] = sample(set, subset_length)
  }
  stop = cc = 0
  while (stop == 0 & cc <= 50){
    
    left_out = setdiff(set, unique(unlist(out_list)))
    le = ceiling(length(left_out) * 0.7)
    duplic = unique(unlist(out_list)[duplicated(unlist(out_list))])
    if (length(left_out) / expected_left_out - 1 > left_out_perc){
      for (i in 1:n_subset){
        i_dup = which(out_list[[i]] %in% duplic)
        s = ceiling(length(i_dup) * 0.7)
        s = min(c(le, s))
        set.seed(seed_list[i])
        i_sub = sample(i_dup, s)
        set.seed(seed_list[i])
        out_list[[i]][i_sub] = sample(left_out, s)
      }
    } else {
      stop = 1
    }
    cc = cc + 1
  }
  
  # final check
  left_out = setdiff(set, unique(unlist(out_list)))
  if (length(left_out) / expected_left_out - 1 > left_out_perc +  0.3){cat(paste0('\n     --- ', lab, ': ', round(length(left_out) / length(set) * 100), '% elements left out'))}
  for (i in 1:n_subset){
    n = uniqueN(out_list[[i]])
    if (n != length(out_list[[i]])){
      cat(paste0('\n     --- ', lab, ': duplicates in set', i))
    }
    if (n != subset_length){
      cat(paste0('\n     --- ', lab, ': shorter length in set ', i, ' (', subset_length - n, ' elements)'))
    }
  }
  
  return(out_list)
}

# evaluate missing recovering test
recov_test = function(n_repeat, year_match, seed_list, recov_met, rm_perc,   # fixed input
                      res_recov_test, err_log, recov_test_random_set,        # input to be returned as output
                      df_type, df_t, flag_comparison, ind_list = c()         # function parameters
){
  
  # flag_comparison: T or F
  # ind_list: list containing NA_NEW_ind column and year and rep for matching (recov_met, rm_ perc and df_type
  #            are already consistent in the function call loop)
  
  # recover missing year by year
  if (recov_met != 'TENSOR_BF'){
    
    cat('\n\n      --------------  data:', df_type, '\n\n')
    
    for (yr in year_match$year){
      
      slice_o = df_t %>% filter(year == yr) %>% arrange(country) %>% select(-year, -country) %>%
        `rownames<-`(sort(unique(df_t$country)))
      na_ind = is.na(slice_o)
      slice_pair =  as.data.frame(as.table(as.matrix(slice_o)), stringsAsFactors = F) %>%
        left_join(as.data.frame(as.table(na_ind), stringsAsFactors = F), Joining, by = c("Var1", "Var2")) %>%
        setNames(c('country', 'variable', 'val', 'NA_ind')) %>%
        mutate(REF = c(1:n()))
      remove_candidate = (slice_pair %>%
                            filter(NA_ind == F))$REF
      
      if (flag_comparison == F){
        # evaluate set
        remove_set = get_random_set(set = remove_candidate,
                                    n_subset = n_repeat,
                                    subset_length = round(nrow(slice_pair) * rm_perc),
                                    seed_list = seed_list + as.numeric(yr),  # different seed also for each year
                                    lab = paste(recov_met, df_type, yr)) # different seed also for each year
        
        recov_test_random_set[[df_type]][[recov_met]][[yr]] = remove_set
      }
      
      for (rep in 1:n_repeat){
        
        # set new NAs
        if (flag_comparison == T){ # recover missing index from standalone subset and apply on Original
          match_list = slice_pair %>%
            left_join(ind_list %>%
                        filter(year == yr & repetition == rep & method == recov_met) %>%
                        select(country, variable, NEW_NA_ind), by = c("country", "variable")) %>%
            filter(!is.na(NEW_NA_ind))
          rm_ind = match_list$REF
        } else {
          rm_ind = remove_set[[rep]]
        }
        
        slice = slice_pair %>%
          mutate(NA_REF = ifelse(REF %in% rm_ind, 1, 0)) %>%
          mutate(val = ifelse(NA_REF == 1, NA, val)) %>%
          select(variable, country, val) %>%
          spread(variable, val) %>%
          arrange(country) %>%
          select_(.dots = colnames(slice_o)) %>%
          `rownames<-`(rownames(slice_o))
        
        slice_recov = c()
        
        # NIPALS
        if (recov_met == 'NIPALS'){
          nip = try_nipals(slice)
          if (!is.null(nip$no_error) | !is.null(nip$warn_text)){
            slice_recov = slice
            slice_recov[na_ind] = nip$value$fitted[na_ind]
            if (!is.null(nip$warn_text)){cat('\n ---- Warning in NIPALS - data', df_type, '- year:', yr, '- rep:', rep , '\n', nip$warn_text)}
          } else if(!is.null(nip$error_text)){
            slice_recov = slice
            cat('\n #### Error in NIPALS - data', df_type, '- year:', yr, '- rep:', rep , '\n', nip$error_text)
          }
          err = sum(slice_recov == slice, na.rm = T) + sum(na_ind) + length(rm_ind) - nrow(slice) * ncol(slice)
          if (err != 0){cat('\n #### Error on non missing elements in NIPALS - year', yr, '- data', df_type, ':', err)}
          if (sum(is.na(slice_recov)) > 0){err_log =c(err_log, paste0(c(recov_met, df_type, 'All', rep, sum(is.na(slice_recov))), collapse = '  '))}
        }
        
        # SOFT IMPUTE - uses https://arxiv.org/pdf/1410.2596.pdf - https://cran.r-project.org/web/packages/softImpute/vignettes/softImpute.html
        if (recov_met == 'SOFT_IMPUTE'){
          sf = try_soft_imp(slice)
          if (!is.null(sf$no_error) | !is.null(sf$warn_text)){
            slice_recov = sf$value
            if (!is.null(sf$warn_text)){cat('\n ---- Warning in SOFT IMPUTE - data', df_type, '- year:', yr, '- rep:', rep , '\n', sf$warn_text)}
          } else if(!is.null(sf$error_text)){
            slice_recov = slice
            cat('\n #### Error in SOFT IMPUTE - data', df_type, '- year:', yr, '- rep:', rep , '\n', sf$error_text)
          }
          err = sum(slice_recov == slice, na.rm = T) + sum(na_ind) + length(rm_ind) - nrow(slice) * ncol(slice)
          if (err != 0){cat('\n #### Error on non missing elements in SOFT IMPUTE - year', yr, '- data', df_type, ':', err)}
          if (sum(is.na(slice_recov)) > 0){err_log =c(err_log, paste0(c(recov_met, df_type, 'All', rep, sum(is.na(slice_recov))), collapse = '  '))}
        }
        
        # save results
        res_recov_test = res_recov_test %>% rbind(
          slice_pair %>% mutate(NEW_NA_ind = ifelse(REF %in% rm_ind, 1, 0)) %>%
            left_join(as.data.frame(as.table(as.matrix(slice_recov)), stringsAsFactors = F) %>%
                        setNames(c('country', 'variable', 'val_recov')), by = c("country", "variable")) %>%
            mutate(year = yr,
                   method = recov_met,
                   data = df_type,
                   repetition = rep,
                   remov_perc = rm_perc)
        )
      } # rep
    } # yr
    
    # recover missing for all years together
    # TENSOR BF - https://www.biorxiv.org/content/biorxiv/early/2016/12/29/097048.full.pdf
  } else if (recov_met == 'TENSOR_BF'){
    
    # prepare slices for each repetition
    res_tens = c()
    slice_rep = rep(list(c()), n_repeat)
    
    # prepare slices with missing
    for (yr in year_match$year){
      
      slice_o = df_t %>% filter(year == yr) %>% arrange(country) %>% select(-year, -country) %>%
        `rownames<-`(sort(unique(df_t$country)))
      na_ind = is.na(slice_o)
      slice_pair =  as.data.frame(as.table(as.matrix(slice_o)), stringsAsFactors = F) %>%
        left_join(as.data.frame(as.table(na_ind), stringsAsFactors = F), Joining, by = c("Var1", "Var2")) %>%
        setNames(c('country', 'variable', 'val', 'NA_ind')) %>%
        mutate(REF = c(1:n()))
      remove_candidate = (slice_pair %>%
                            filter(NA_ind == F))$REF
      
      if (flag_comparison == F){
        # evaluate set
        remove_set = get_random_set(set = remove_candidate,
                                    n_subset = n_repeat,
                                    subset_length = round(nrow(slice_pair) * rm_perc),
                                    seed_list = seed_list + as.numeric(yr),  # different seed also for each year
                                    lab = paste(recov_met, df_type, yr))
        
        recov_test_random_set[[df_type]][[recov_met]][[yr]] = remove_set
      }
      
      for (rep in 1:n_repeat){
        
        # set new NAs
        if (flag_comparison == T){ # recover missing index from standalone subset and apply on Original
          match_list = slice_pair %>%
            left_join(ind_list %>%
                        filter(year == yr & repetition == rep & method == recov_met) %>%
                        select(country, variable, NEW_NA_ind), by = c("country", "variable")) %>%
            filter(!is.na(NEW_NA_ind))
          rm_ind = match_list$REF
        } else {
          rm_ind = remove_set[[rep]]
        }
        
        slice = slice_pair %>%
          mutate(NA_REF = ifelse(REF %in% rm_ind, 1, 0)) %>%
          mutate(val = ifelse(NA_REF == 1, NA, val)) %>%
          select(variable, country, val) %>%
          spread(variable, val) %>%
          arrange(country) %>%
          select_(.dots = colnames(slice_o)) %>%
          `rownames<-`(rownames(slice_o))
        
        slice_rep[[rep]][[toString(yr)]] = slice
        res_tens = res_tens %>% rbind(slice_pair %>%
                                        mutate(NEW_NA_ind = ifelse(REF %in% rm_ind, 1, 0),
                                               year = yr,
                                               method = recov_met,
                                               data = df_type,
                                               repetition = rep,
                                               remov_perc = rm_perc))
      } # rep
    } # yr
    
    # create tensor and recover missing
    tens_rec = c()
    for (rep in 1:n_repeat){
      
      # create tensor for each repetition and recover missing
      tens = array(numeric(),c(uniqueN(df_t$country), ncol(df_t) - 2, max(year_match$CODE)))
      for (i in year_match$CODE){
        tens[, , i] = slice_rep[[rep]][[toString((year_match %>% filter(CODE == i))$year)]] %>% as.matrix()
        
      }
      
      cat('\n\n      --------------  data:', df_type, '   rep:', rep, '\n\n')
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
      if (sum(is.na(tens_recov)) > 0){err_log =c(err_log, paste0(c(recov_met, df_type, 'All', rep, sum(is.na(tens_recov))), collapse = '  '))}
      
      # reshape results
      for (i in year_match$CODE){
        tens_rec = tens_rec %>% rbind(
          as.data.frame(as.table(as.matrix(tens_recov[, , i]  %>% `colnames<-`(colnames(slice)) %>% `rownames<-`(rownames(slice)))), stringsAsFactors = F) %>%
            setNames(c("country", "variable", "val_recov")) %>%
            mutate(year = (year_match %>% filter(CODE == i))$year,
                   repetition = rep)
        )
      }
    } # rep
    res_tens = res_tens %>% left_join(tens_rec, by = c("country", "variable", "year", "repetition"))
    if (sum(is.na(res_tens$val_recov)) - sum(is.nan(res_tens$val_recov)) > 0){cat('\n\n\n\n\n\n #########  TENSOR_BF - NAs in in left_join for', df_type, '\n\n\n\n\n\n')}
    
    res_recov_test = res_recov_test %>% rbind(res_tens)
  } # else tensor
  
  return(list(res_recov_test = res_recov_test,
              err_log = err_log,
              recov_test_random_set = recov_test_random_set))
}

# add additional variables
add_variables = function(df, df_test, df_test_lab, NA_toll, summary_fold){
  
  df_ref = df %>% select(country, year)
  
  cat('\n\n\n###########     ', df_test_lab, '\n')
  if (sum(is.na(df_test$country)) > 0){cat(' *** NAs in country')}
  if (sum(is.na(df_test$year)) > 0){cat(' *** NAs in year')}
  cat('\n years range: ', range(df_test$year))
  
  summary_tab = c()
  
  df_test = nameCleaning(df_test)
  
  for (var in setdiff(colnames(df_test), c('country', 'year'))){
    
    new_lab = paste0('var_', df_test_lab, '_', var)
    var_test = df_test %>%
      # setnames(var, new_lab) %>%
      rename_at(vars(var), funs(str_replace(., var, new_lab))) %>%
      select_(.dots = c('country', 'year', new_lab))
    
    merge_candid = df_ref %>%
      left_join(var_test, by = c("country", "year"))
    
    NA_perc = sum(is.na(merge_candid %>% select_(.dots = new_lab))) / nrow(merge_candid)
    
    NA_by_year = merge_candid %>%
      group_by(year) %>%
      summarise(TOT_NA = sum(is.na(!!as.name(new_lab))))
    
    keep = 1
    NA_by_year_ind = which(NA_by_year$TOT_NA == uniqueN(merge_candid$country))
    if (length(NA_by_year_ind) > 0){
      # cat('\n    REMOVED - missing all:', paste0(NA_by_year$year[NA_by_year_ind], collapse = ", "), '   ', var)
      keep = 0
      summary_tab = rbind(summary_tab, c('REMOVED', var, paste0('missing ', paste0(NA_by_year$year[NA_by_year_ind], collapse = ", "))))
    }
    if (NA_perc > NA_toll / 100){
      # cat('\n    REMOVED - NAs %:', round(NA_perc * 100, 2), '   ', var)
      keep = 0
      summary_tab = rbind(summary_tab, c('REMOVED', var, paste0('NA% ',round(NA_perc * 100, 2))))
    }
    
    if (keep == 1){
      df = df %>%
        left_join(merge_candid, by = c("country", "year"))
      summary_tab = rbind(summary_tab, c('ADDED', var, ''))
    }
    
  } # var
  
  summary_tab = as.data.frame(summary_tab, stringsAsFactors = F) %>%
    setNames(c('ACTION', 'VARIABLE', 'REASON')) %>%
    arrange(ACTION)
  
  write.table(summary_tab, paste0(summary_fold, '3_Additional_variable_', df_test_lab, '.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  cat('\n\nVariables kept: ', sum(summary_tab$ACTION == 'ADDED'))
  cat('\nVariables removed: ', sum(summary_tab$ACTION == 'REMOVED'), '\n')
  
  return(df)
}

# Standard cleaning procedure
nameCleaning = function(df){
  
  # Clean unwanted chars
  # names(df) = gsubfn(paste(names(unwanted_array),collapse='|'), unwanted_array, names(df))
  
  # Custom cleaning
  names(df) = gsub("[\\. \\(\\)\\/]+", "_", names(df))
  names(df) = gsub("-", "_", names(df))
  names(df) = gsub("–", "_", names(df)) # bigger
  names(df) = gsub("'", "", names(df))
  names(df) = gsub(",", "_", names(df))
  names(df) = gsub(":", "_", names(df))
  names(df) = gsub("<", "MIN", names(df))
  names(df) = gsub(">", "MAG", names(df))
  names(df) = gsub("&", "E", names(df))
  names(df) = gsub("°", "", names(df))
  names(df) = gsub("=", "", names(df))
  names(df) = gsub(";", "", names(df))
  names(df) = gsub("\\*", "", names(df))
  names(df) = gsub("’", "", names(df))
  names(df) = gsub("%", "PERC", names(df))
  names(df) = gsub("\\+", "_", names(df))
  names(df) = gsub("\\$","DOLL", names(df))
  
  
  # To upper
  # names(df) = toupper(names(df))
  
  # Trim
  names(df) = trimws(names(df))
  
  # Cut recurring underscore
  names(df) = gsub("_+", "_", names(df))
  names(df) = gsub('^_|_$', '', names(df))
  
  return(df)
}

# plot variable importance
plot_loadings = function(load_list, load_range, dim_to_show_lab, add_title = '', err_bar = F){
  # load_list is a data.frame with column "Variable", "Importance"
  
  p = ggplot(load_list,
             aes(x=Variable, y = Importance)) +
    geom_bar(stat = "identity", position = 'dodge', fill = 'deepskyblue1') +
    ylim(load_range[1], load_range[2]) +
    coord_flip() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          # plot.margin = margin(2, 0, 0.7, 0.3, "cm"),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major.x = element_line(colour = "black", linetype = 'dashed', size = 0.2),
          plot.title = element_text(size = 20),
          legend.position = 'none') +
    ggtitle(paste0('Variable Importance for Index ', dim_to_show_lab, add_title)) +
    annotate("text", x = c(1:nrow(load_list)), y = 0, label = substr(load_list$Variable, 5, 200), hjust = 0.5, vjust = 0.5, size = 4, fontface = 'bold')
  
  if (err_bar){
    p = p +
      geom_errorbar(aes(ymin = load_min, ymax = load_max, width = 0.6, alpha = 0.5))
  }
  
  return(p)
}

# plot index on 1 or 2 dimension
plot_index = function(index_list, index_1_thresh, index_2_thresh, dim_to_show, year, expl_var = c(), add_title_variance = ''){
  
  # colnames(index_list) =  "INDEX", "I_1", "I_2", where INDEX contains the observations names, I_1, I_2 are the first two
  #                                                dimension to plot according to dim_to_show
  #                                                expl_var is used only for the title in PCA-like plot
  
  if (dim_to_show == 1){
    p = ggplot(index_list %>%
                 select(INDEX, I_1) %>%
                 rename(Index = I_1) %>%
                 arrange(Index) %>%
                 mutate(color = as.factor(ifelse(Index > index_1_thresh, 1, 0))),
               aes(x=0, y=Index, label = INDEX)) +
      geom_hline(yintercept = 0, size = 2) +
      geom_point() +
      geom_label_repel(aes(x=0, y=Index, label = INDEX, fill = color)) +
      coord_flip() +
      scale_fill_manual(values = c('cyan3', 'chocolate1')) +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.major.x = element_line(colour = "black", linetype = 'dashed', size = 0.2),
            panel.grid.minor.x = element_line(colour = "black", linetype = 'dashed', size = 0.2),
            legend.position = "none")
    
  } else if (dim_to_show == 2){
    p = ggplot(index_list %>%
                 select(INDEX, I_1, I_2) %>%
                 rename(Index1 = I_1,
                        Index2 = I_2) %>%
                 mutate(color1 = ifelse(Index1 > index_1_thresh, 1, 0),
                        color2 = ifelse(Index2 > index_2_thresh, 1, 0),
                        color = paste0(color1, '_', color2)),
               aes(x=Index1, y=Index2, label = INDEX)) +
      geom_hline(yintercept = 0, size = 2) +
      geom_vline(xintercept = 0, size = 2) +
      geom_point() +
      geom_label_repel(aes(x=Index1, y=Index2, label = INDEX, fill = color)) +
      scale_fill_manual(values = c('cyan3', 'chocolate1', 'chartreuse2', 'cornsilk2')) +
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.major = element_line(colour = "black", linetype = 'dashed', size = 0.2),
            panel.grid.minor = element_line(colour = "black", linetype = 'dashed', size = 0.2),
            legend.position = "none")
  }
  
  # add title
  if (length(expl_var) > 0){
    p = p + ggtitle(paste0(year, ' - ', dim_to_show, ifelse(dim_to_show == 1, ' Index', ' Indexes'), '   (Cumulative Explained Variance', add_title_variance, ': ',
                           round(expl_var * 100), '%)')) +
      theme(plot.title = element_text(size = 20))
  } else {
    p = p + ggtitle(paste0(year, ' - ', dim_to_show, ifelse(dim_to_show == 1, ' Index', ' Indexes'))) +
      theme(plot.title = element_text(size = 20))
  }
  
  return(p)
}

# evaluate index for PCA
evaluate_index_PCA = function(res_index_PCA, res_PCA_list, res_PCA_loadings, res_PCA_importance, PC_to_keep, index_1_thresh, index_2_thresh,
                              load_thresh, leading_var, leading_sign, df_type, recov_met, pc_met){
  
  load_ref = sort_loading_PCA(res_PCA_loadings, leading_var, leading_sign) %>%
    filter(data == df_type) %>%
    filter(method == recov_met) %>%
    filter(PCA == pc_met) %>%
    filter(PC <= PC_to_keep)
  load_range = range(load_ref$loading)
  
  expl_var_out = c()
  scores_out_raw = scores_out_index = data.frame(COUNTRY = c('INDEX', rownames(res_PCA_list[[df_type]][[recov_met]][[1]][[pc_met]][['pca']][['x']])), stringsAsFactors = F)
  load_out = data.frame(VARIABLE = c('VARIABLE', rownames(res_PCA_list[[df_type]][[recov_met]][[1]][[pc_met]][['pca']][['rotation']])), stringsAsFactors = F)
  row_list = list()
  i = 1
  for (yr in sort(unique(load_ref$year))){
    
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
    for (pc in 1:min(c(PC_to_keep, 2))){
      
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
      
      row_list[[toString((pc - 1) * min(c(PC_to_keep, 2)) + 1)]][[i]] = ggplotGrob(p)
      row_list[[toString(min(c(PC_to_keep, 2)) * pc)]][[i]] = ggplotGrob(p_load)
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
  png(paste0('./Results/1_Index_plot_thrs', load_thresh, '_', df_type, '_', recov_met, '_', pc_met, '.png'), width = 22.5 * PC_to_keep, height = 60, units = 'in', res=300)
  grid.draw(g)
  dev.off()
  
  expl_var_out = expl_var_out %>%
    filter(PC == PC_to_keep)
  colnames(scores_out_raw) = gsub('\\.1', '', colnames(scores_out_raw))
  write.table(scores_out_raw, paste0('./Results/1_Index_trhs', load_thresh, '_', df_type, '_', recov_met, '_', pc_met, '_scores_raw.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  colnames(scores_out_index) = gsub('\\.1', '', colnames(scores_out_index))
  write.table(scores_out_index, paste0('./Results/1_Index_trhs', load_thresh, '_', df_type, '_', recov_met, '_', pc_met, '_scores_index.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  write.table(load_out, paste0('./Results/1_Index_trhs', load_thresh, '_', df_type, '_', recov_met, '_', pc_met, '_loadings.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  # table to help understanding of index, showing countries that change group in different years
  # help_tab = c()
  # for (j in 2:ncol(scores_out_raw)){
  #   help_tab = rbind(help_tab, cbind(COUNTRY = scores_out_raw[-1, 'COUNTRY'], INDEX = scores_out_raw[1, j],
  #                                    YEAR = colnames(scores_out_raw[j]), VAL = scores_out_raw[-1, j]))
  # }
  # help_tab = help_tab %>%
  #   as.data.frame(stringsAsFactors = F) %>%
  #   mutate(VAL = ifelse(VAL > 0, 1, 0)) %>%
  #   filter(YEAR != 'Avg') %>%
  #   group_by(COUNTRY, INDEX, VAL) %>%
  #   summarize(OCCURENCY = n())  # must be changed dynamically with threshold
  
  # write.table(help_tab, paste0('./Results/1_Index_trhs', load_thresh, '_', df_type, '_', recov_met, '_', pc_met, '_scores_over_time.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  res_index_PCA[[df_type]][[recov_met]][[pc_met]][['scores_raw']] = scores_out_raw
  res_index_PCA[[df_type]][[recov_met]][[pc_met]][['scores_index']] = scores_out_index
  res_index_PCA[[df_type]][[recov_met]][[pc_met]][['loadings']] = load_out
  res_index_PCA[[df_type]][[recov_met]][[pc_met]][['Expl_Var']] = expl_var_out
  # res_index_PCA[[df_type]][[recov_met]][[pc_met]][['scores_over_time']] = help_tab
  
  return(res_index_PCA)
}

# function to be called with try catch - Dynamic Factor Model
t_DFM = function(data, n_factor){
  out_message = capture.output(out_res <- psychNET(
    data = data,
    model = 'DFM',
    lag = 1,
    nFact = n_factor
  ),
  type = 'output'
  )
  
  return(list(out_message = out_message,
              out_res = out_res))
}
try_DFM = function(data, n_factor){
  r <- tryCatch({
    list(value = t_DFM(data, n_factor), no_error = "No error")
  }, warning = function(e) {
    warn_text <- paste0("WARNING: ", e)
    return(list(value = (suppressWarnings(t_DFM(data, n_factor))), warn_text = warn_text))
  }, error = function(e) {
    error_text <- paste0("ERROR: ", e)
    return(list(value = NA, error_text = error_text))
  }, finally = {
  }, quiet = TRUE)
  return(r)
}

# adds column or rows 
add_row_col = function(mat, rows_to_add = NULL, cols_to_add = NULL, filling_val = 0, total_new_row = NULL, total_new_col = NULL){
  # rows_to_add and cols_to_add must be a vector of integers
  # rows and column will be filled with filling_val
  
  if (is.null(total_new_row)){total_new_row = nrow(mat)}
  if (is.null(total_new_col)){total_new_col = ncol(mat)}
  
  if (nrow(mat) + length(rows_to_add) != total_new_row){cat('\n*** Check number of rows to add or total_new_row\n')}
  if (ncol(mat) + length(cols_to_add) != total_new_col){cat('\n*** Check number of cols to add or total_new_col\n')}
  
  mat_new = c()
  
  if (!is.null(rows_to_add)){
    keep_rows = setdiff(c(1:total_new_row), rows_to_add)[1:nrow(mat)]
    r = 1
    for (row in 1:total_new_row){
      
      if (row %in% keep_rows){
        mat_new = mat_new %>% rbind(mat[r, ])
        r = r + 1
      } else{
        mat_new = mat_new %>% rbind(rep(filling_val, ncol(mat)))
      }
    }
  }
  
  if (!is.null(cols_to_add)){
    keep_cols = setdiff(c(1:total_new_col), cols_to_add)[1:ncol(mat)]
    c = 1
    if (!is.null(mat_new)){
      mat = mat_new
      mat_new = c()
    }
    
    for (col in 1:total_new_col){
      
      if (col %in% keep_cols){
        mat_new = mat_new %>% cbind(mat[, c])
        c = c + 1
      } else{
        mat_new = mat_new %>% cbind(rep(filling_val, nrow(mat)))
      }
    }
  }
  
  return(mat_new)
}

# evaluate Dynamic Factor Model for single country
evaluate_DFM_univar = function(res_DFM_factors, res_DFM_loadings, res_DFM_stats, res_DFM_list, res_DFM_list_reload, variable_code, country_code, df_SP,
                               n_factor, country_i, df_type, recov_met, univ_reload = F, dfm_eval = 'BFGS', dfm_max_iter = 1000){
  
  # formula x_t = A*x_t-1 + N(0,Q), y_t = C*x_t + N(0,R), implemented with MARSS - https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf
  # dfm_eval = 'BFGS' or 'kem' for EM
  
  # first row is the oldest date (e.g. 1980), last row is the newest date (e.g. 2010)
  country_lab = (country_code %>% filter(ID == country_i))$country
  data_country = df_SP %>% filter(ID == country_i) %>% arrange(TIME)
  year_code = data_country %>% select(TIME) %>% mutate(BEEP = as.character(1:nrow(data_country))) %>% rename(year = TIME)
  
  cat('\n Evaluating', country_lab, paste0('(', country_i, '/', nrow(country_code), ')'))
  
  # check for constant variables
  const_check = data_country %>%
    select(-ID, -TIME) %>%
    gather('variable', 'val') %>%
    group_by(variable) %>%
    summarise(COUNT = uniqueN(val), .groups = 'drop') %>%
    filter(COUNT == 1)
  
  var_to_remove_lab = ''
  var_to_remove = c()
  if (nrow(const_check) > 0){
    var_to_remove = const_check$variable
    var_to_remove_lab = paste0((variable_code %>% filter(CODE %in% var_to_remove))$variable, collapse = ', ')
    data_country = data_country %>% select(-all_of(var_to_remove))
    cat('\n    *** removed (constant):', var_to_remove_lab)
  }
  
  # time series on rows, time on column and standardize
  data_work = data_country %>%
    select(-ID, -TIME) %>%
    t() %>%
    MARSS::zscore()
  
  # set up DFM matrices
  Z = matrix(
    (expand.grid(1:nrow(data_work),1:n_factor) %>%
       mutate(elem = paste0('z.(', Var1, ',', Var2, ')')))$elem,
    ncol = n_factor, byrow=FALSE)
  Z[upper.tri(Z, diag = FALSE)] = 0 # upper triangular part set to 0 for identifiability (estimated coeff will be invariant to data transformation)
  Z = as.list(Z)
  Z = matrix(lapply(Z, function(x){if (x == '0'){x = 0} else {x}}),  # 0 must be numeric and not character
             ncol = n_factor, byrow=FALSE)
  
  R = 'diagonal and unequal'
  
  B = matrix(
    (expand.grid(1:n_factor,1:n_factor) %>%
       mutate(elem = paste0('b.(', Var1, ',', Var2, ')')))$elem,
    ncol = n_factor, byrow=FALSE)
  
  Q = diag(1, n_factor) # for identifiability
  
  x0 = U = A = "zero"
  
  V0 = diag(nrow(data_work), n_factor) # starting variance for X0 distribution (normal)
  
  dfa.model = list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0)
  cntl.list = list(maxit=dfm_max_iter)
  
  # reload or evaluate
  warn_catch = NULL
  reload_check = res_DFM_list_reload[[df_type]][[recov_met]][['DFM_univar']][[paste0(n_factor, '_factors')]][[country_lab]]
  if (is.null(reload_check) | !univ_reload){
    start_time = Sys.time()
    options(warn=1)
    warn_catch = capture.output(
      dfm_fit <- MARSS(data_work, model=dfa.model, control=cntl.list, method=dfm_eval, silent = T, fun.kf = "MARSSkfas"),
      type = 'message')
    options(warn=0)
    cat(' Done in', round((as.numeric(Sys.time())-as.numeric(start_time)) / 60), 'mins ')
  } else {
    cat(' - reloaded')
    dfm_fit = reload_check$dfm_fit
  }
  conv = dfm_fit$convergence
  
  # see MARSSoptim for further info
  if (conv > 0){conv_mex = paste0('Error to be checked: ', conv)}
  if (dfm_eval == 'BFGS'){
    if (conv == 1){conv_mex = 'Maximum number of iteration reached'}
    if (conv == 10){conv_mex = 'Some of the variance elements appear to be degenerate'}
    if (conv == 52){conv_mex = 'The algorithm was abandoned due to errors from the "L-BFGS-B" method'}
    if (conv == 53){conv_mex = 'The algorithm was abandoned due to numerical errors in the likelihood calculation'}
  }
  
  if (conv == 0){
    cat(' - converged')
    dfm_par = coef(dfm_fit,  type="matrix")
    factors = t(dfm_fit$states)
    A = dfm_par$B # factor matrix
    C = dfm_par$Z # loading matrix
    Q = dfm_par$Q # factors covariance
    R = dfm_par$R # observations covariance
    TSS_val = (data_work - mean(data_work)) ^ 2
    TSS = sum(TSS_val)
    RSS_err = try(capture.output(RSS_val <- suppressWarnings(residuals(dfm_fit)$model.residuals ^ 2)), silent = T) # suppress warning message for std
    if (class(RSS_err) == "try-error"){
      RSS = RSS_95 = RSS_99 = sum(TSS_val) / 2
      TSS_95 = TSS_99 = TSS
      cat('\n    *** unable to evaluate RSS -> R^2 set to 50%')
      Explain_var = Explain_var_95 = Explain_var_99 = 0.5
    } else {
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
    }
    dfm_AICc = dfm_fit$AIC
    dfm_logLik = dfm_fit$logLik
    
    # varimax rotation on loadings
    if (n_factor > 1){
      H.inv = varimax(C)$rotmat
      C = C %*% H.inv
      factors =  solve(H.inv) %*% t(factors) %>% t()
    }
    
    if (nchar(var_to_remove_lab) > 0){
      C = add_row_col(mat = C, rows_to_add = as.numeric(var_to_remove), filling_val = 0, total_new_row = nrow(variable_code))
      R = add_row_col(mat = R, rows_to_add = as.numeric(var_to_remove), cols_to_add = as.numeric(var_to_remove),
                      filling_val = 0, total_new_row = nrow(variable_code), total_new_col = nrow(variable_code))
    }
  } else {
    cat('\n    ###### Error in convergence:', conv_mex)
    factors = matrix(NA, nrow = nrow(year_code), ncol = n_factor)
    C = matrix(NA, nrow = nrow(data_work), ncol = n_factor)
    R = matrix(NA, nrow = nrow(data_work), ncol = nrow(data_work))
    A = Q = matrix(NA, nrow = n_factor, ncol = n_factor)
    RSS = TSS = dfm_AICc = dfm_logLik = NA
  }
  if (length(warn_catch) > 0){cat('\n    *** Warnings:', warn_catch)}
  
  # DFM fitted with psychNET package - doesn't return AIC and LogLik
  {
    # # https://rdrr.io/github/guilbran/dynfactoR/f/inst/doc/dynamic-factors.pdf
    # tryDFM = try_DFM(data = data_country, n_factor = n_factor)
    # 
    # if (!is.null(tryDFM$no_error) | !is.null(tryDFM$warn_text)){
    #   
    #   converg = tryDFM$value$out_message[4]
    #   DFM = tryDFM$value$out_res
    #   
    #   year_code = DFM$CALL$pars$other_data %>% select(BEEP, TIME) %>% rename(year = TIME) %>% mutate(BEEP = as.character(BEEP))
    #   factors = DFM$fit$qml
    #   A = DFM$fit$A
    #   C = DFM$fit$C
    #   Q = DFM$fit$Q
    #   R = DFM$fit$R
    #   
    #   if (nchar(var_to_remove_lab) > 0){
    #     C = add_row_col(mat = C, rows_to_add = as.numeric(var_to_remove), filling_val = 0, total_new_row = nrow(variable_code))
    #     R = add_row_col(mat = R, rows_to_add = as.numeric(var_to_remove), cols_to_add = as.numeric(var_to_remove),
    #                     filling_val = 0, total_new_row = nrow(variable_code), total_new_col = nrow(variable_code))
    #   }
    #   
    #   cat('\n ', ifelse(nchar(var_to_remove_lab) > 0, '    ', ''), converg, '-', (country_code %>% filter(ID == country_i))$country) 
    #   
    # } else if(!is.null(tryDFM$error_text)){
    #   cat('\n ####', country_lab, '- ', tryDFM$error_text)
    #   factors = A = C = Q = R = NULL
    # }
  }
  
  # save results
  res_DFM_stats = res_DFM_stats %>% bind_rows(
    data.frame(DFM = 'DFM_univar',
               method = recov_met,
               data = df_type,
               Total_Factors = n_factor,
               country = country_lab,
               AICc = dfm_AICc,
               LogLik = dfm_logLik,
               y_RSS = RSS,
               y_RSS_95 = RSS_95,
               y_RSS_99 = RSS_99,
               y_TSS = TSS,
               y_TSS_95 = TSS_95,
               y_TSS_99 = TSS_99,
               Explain_var = Explain_var,
               Explain_var_99 = Explain_var_99,
               Explain_var_95 = Explain_var_95,
               algo_err_code = ifelse(conv > 0, conv, ''),
               removed_variable = ifelse(length(var_to_remove) > 0, toString(length(var_to_remove)), ''), stringsAsFactors = F)
  )
  
  res_DFM_list[[df_type]][[recov_met]][['DFM_univar']][[paste0(n_factor, '_factors')]][[country_lab]] = list(
    factors = factors %>% `rownames<-`(data_country$TIME),
    A = A,
    C = C %>% `rownames<-`(variable_code$variable),
    Q = Q,
    R = R %>% `rownames<-`(variable_code$variable) %>% `colnames<-`(variable_code$variable),
    dfm_fit = dfm_fit,
    RSS = RSS,
    TSS = TSS,
    dfm_AICc = dfm_AICc,
    dfm_logLik = dfm_logLik
  )
  
  res_DFM_factors = res_DFM_factors %>% bind_rows(
    as.data.frame(factors, stringsAsFactors = F) %>%
      setNames(c(1:n_factor)) %>%
      rownames_to_column(var = 'BEEP') %>%
      left_join(year_code, by = "BEEP") %>%
      select(-BEEP) %>%
      gather('Factor', 'val', -year) %>%
      mutate(Factor = as.numeric(Factor),
             method = recov_met,
             data = df_type,
             country = country_lab,
             Total_Factors = n_factor,
             Var_removed = var_to_remove_lab,
             DFM = 'DFM_univar') %>%
      select(DFM, country, Total_Factors, Factor, everything())
  )
  
  res_DFM_loadings = res_DFM_loadings %>% bind_rows(
    as.data.frame(C) %>%
      setNames(c(1:n_factor)) %>%
      mutate(variable = variable_code$variable) %>%
      gather('Factor', 'loading', -variable) %>%
      mutate(country = country_lab,
             Total_Factors = n_factor,
             method = recov_met,
             data = df_type,
             DFM = 'DFM_univar') %>%
      select(DFM, country, variable, Total_Factors, Factor, everything())
  )
  
  return(list(
    res_DFM_stats = res_DFM_stats,
    res_DFM_list = res_DFM_list,
    res_DFM_factors = res_DFM_factors,
    res_DFM_loadings = res_DFM_loadings
  ))
}

# adjust DFM taking into account all countries together
evaluate_DFM_multivar = function(res_DFM_factors, res_DFM_loadings, res_DFM_stats, res_DFM_list, res_DFM_MAPE,
                                 df_SP_orig, df_type, recov_met, n_factor, VAR_alpha = 0.9, kalm_Q_hat_mode, true_factors = NULL,
                                 res_DFM_list_reload = NULL, multiv_reload_VAR = F, multiv_reload_Kalm = F){
  
  # 1 - Fit a VAR model for x_t factor/states process, stacking all estimated x_t from all countries -> get A* from x_t = A*x_t-1 + N(0, Q*)
  # 2 - Create C* as block diagonal matrix with diag(C*) = [C1, ..., Cn], where Ci is the estimated loading for each country
  # 3 - perform Kalman Filter to filter x_t with A* and C*
  
  # VAR_alpha: penalty for sparseness: 1 is LASSO, 0 is Ridge
  # kalm_Q_hat_mode: select variance matrix for state variable (factors) when filtering with Kalman
  #                  'from VAR': takes residual covariance from VAR model on x_t,  'identity': identity matrix
  # true_factors: if testing Kalman filter's reconstruction ability, provide true factors used to generate y_t in the form of factor_work, i.e.
  #               time on row (e.g. 2010 top, 2017 bottom), pairs country-#factor on column (e.g. "Albania_1", "Albania_2", etc)
  
  # select data
  year_order = sort(unique((res_DFM_factors %>%
                             filter(data == df_type) %>%
                             filter(method == recov_met) %>%
                             filter(DFM == 'DFM_univar') %>%
                             filter(Total_Factors == n_factor))$year))
  variable_order = sort(unique(res_DFM_loadings$variable))
  country_order = sort(unique(res_DFM_factors$country))
  country_fact_order = (expand.grid(country_order, paste0('_', 1:n_factor)) %>%
                          mutate(COL = paste0(Var1, Var2)) %>%
                          select(COL) %>%
                          arrange(COL))$COL
  country_variab_order = (expand.grid(country_order, paste0('_', variable_order)) %>%
                            mutate(COL = paste0(Var1, Var2)) %>%
                            select(COL) %>%
                            arrange(COL))$COL
  
  loading_work = res_DFM_loadings %>%
    filter(data == df_type) %>%
    filter(method == recov_met) %>%
    filter(DFM == 'DFM_univar') %>%
    filter(Total_Factors == n_factor)
  
  factor_work = res_DFM_factors %>%
    filter(data == df_type) %>%
    filter(method == recov_met) %>%
    filter(DFM == 'DFM_univar') %>%
    filter(Total_Factors == n_factor) %>%
    arrange(desc(country)) %>%
    mutate(Factor = paste0(country, '_', Factor)) %>%
    select(Factor, year, val) %>%
    spread(Factor, val) %>%
    arrange(year) %>%
    select(-year) %>%
    # ts(start = c(min(res_DFM_factors$year), 1), frequency = 1)
    select(all_of(country_fact_order)) %>%
    t() %>%
    MARSS::zscore() %>%
    t()
  
  # reload for simulation only
  reload_check = res_DFM_list_reload[[df_type]][[recov_met]][['DFM_multivar']][[paste0(n_factor, '_factors')]][['All']]
  
  # select factors for all countries and estimate a single factor matrix with sparse VAR
  cat('\n** evaluating VAR...')
  if (multiv_reload_VAR & !is.null(reload_check)){
    VAR_fit = reload_check[['VAR_fit']]
    cat(' Reloaded \n')
  } else {
    set.seed(42)
    start_time = Sys.time()
    VAR_fit = fitVAR(factor_work,  # variable (time series) on column
                     p = 1,
                     penalty = 'ENET',
                     nfolds = 5,
                     alpha = VAR_alpha)
    cat(' Done in', round((as.numeric(Sys.time())-as.numeric(start_time)) / 60), 'mins \n')
  }
  
  A_hat = VAR_fit$A[[1]]
  A_sparseness = sum(A_hat!=0) / prod(dim(A_hat))
  Q_hat = VAR_fit$sigma
  VAR_RMSE = sqrt(mean(VAR_fit$residuals ^ 2))
  VAR_CV_RMSE = sqrt(VAR_fit$mse)
  VAR_CV_RMSE_SD = sqrt(VAR_fit$mseSD)
  VAR_MAPE_list = VAR_fit$residuals %>%
    `rownames<-`(year_order) %>%
    as.table() %>%
    as.data.frame(stringsAsFactors = F) %>%
    setNames(c('year', 'variable', 'residual')) %>%
    left_join(VAR_fit$series %>%
                `rownames<-`(year_order) %>%
                as.table() %>%
                as.data.frame(stringsAsFactors = F) %>%
                setNames(c('year', 'variable', 'original')),
              by = c("year", "variable")) %>%
    filter(abs(original) > .Machine$double.eps) %>%
    mutate(MAPE = abs(residual / original)) %>%
    arrange(desc(MAPE))
  VAR_MAPE = mean(VAR_MAPE_list$MAPE)
  VAR_MAPE_max = VAR_MAPE_list$MAPE[1]
  VAR_MAPE_95 = mean((VAR_MAPE_list %>% filter(row_number() > ceiling(0.05 * nrow(VAR_MAPE_list))))$MAPE)
  VAR_MAPE_99 = mean((VAR_MAPE_list %>% filter(row_number() > ceiling(0.01 * nrow(VAR_MAPE_list))))$MAPE)
  res_DFM_MAPE = res_DFM_MAPE %>% bind_rows(
    VAR_MAPE_list %>%
      filter(row_number() < ceiling(0.05 * nrow(VAR_MAPE_list))) %>%
      mutate(method = recov_met,
             data = df_type,
             Total_Factors = n_factor,
             algo = 'VAR') %>%
      select(algo, Total_Factors, everything())
    )

  # create C* as block diagonal matrix
  list_work_C = list_work_R = list()
  lc = 1
  for (con in country_order){
    list_work_C[[lc]] = data.frame(variable = variable_order, stringsAsFactors = F) %>%
      left_join(
        loading_work %>%
          filter(country == con) %>%
          select(variable, Factor, loading) %>%
          spread(Factor, loading),
        by = "variable") %>%
      select(-variable) %>%
      select(all_of(1:n_factor)) %>%
      as.matrix()
    
    list_work_R[[lc]] = res_DFM_list[[df_type]][[recov_met]][['DFM_univar']][[paste0(n_factor, '_factors')]][[con]][['R']]
    lc = lc + 1
  }
  
  C_hat = bdiag(list_work_C) %>% as.matrix()
  R_hat = bdiag(list_work_R) %>% as.matrix()
  
  # perform Kalman Filter to filter x_t with A* and C*
  
  y_df = df_SP_orig %>%   # row: pairs country-variable, col: years
    gather('Variable', 'val', -country, -year) %>%
    spread(year, val) %>%
    arrange(country, Variable) %>%
    # filter(country %in% country_order) %>% # todo: rimuovi, serve solo per il subset di 10 nazioni
    mutate(ref = paste0(country, '_', Variable)) %>%
    select(-country, -Variable) %>%
    select(ref, everything()) 
  
  y_mat = y_df %>%  # standardize
    select(-ref)  %>%
    as.matrix() %>%
    MARSS::zscore()
  
  # check for constant value -> remove yt from Kalman filter
  var_to_remove_const = c()
  if (sum(!is.finite(y_mat)) > 0){
    var_to_remove_const = which(apply(y_mat, 1, function(x) sum(!is.finite(x))) > 0)
  }
  
  # loop Kalman Filter untill convergence
  cat('\n** evaluating Kalman Filter...\n')
  if (multiv_reload_Kalm & !is.null(reload_check)){
    kalm_fit = reload_check[['kalm_fit']]
    var_to_remove_all = C_hat_remov = R_hat_remov = c()
    var_to_remove_lab = ''
    y_mat_work = reload_check[['kalm_y_mat_work']]
    cat(' Reloaded \n')
  } else {
    start_time = Sys.time()
    var_to_remove_lab = ''
    if (length(var_to_remove_const) > 0){
      var_to_remove_all = var_to_remove_const
      cat('\n     *** removed for constant values:\n         ')
      cat(paste0(y_df$ref[var_to_remove_const], sep = ','))
    } else {
      var_to_remove_all = c()
    }
    stop = 0
    while (stop == 0){
      # check variables to remove
      if (length(var_to_remove_all) > 0){
        y_mat_work = y_mat[-var_to_remove_all, ] %>% `rownames<-`(setdiff(country_variab_order, y_df$ref[var_to_remove_all]))
        C_hat_remov = C_hat[-var_to_remove_all, ] %>% `rownames<-`(setdiff(country_variab_order, y_df$ref[var_to_remove_all])) %>% `colnames<-`(country_fact_order)
        R_hat_remov = R_hat[-var_to_remove_all, -var_to_remove_all] %>% `rownames<-`(setdiff(country_variab_order, y_df$ref[var_to_remove_all])) %>% `colnames<-`(setdiff(country_variab_order, y_df$ref[var_to_remove_all]))
        C_hat_work = C_hat_remov
        R_hat_work = R_hat_remov
      } else {
        y_mat_work = y_mat
        C_hat_work = C_hat
        R_hat_work = R_hat
        C_hat_remov = R_hat_remov = NULL
        stop = 1
      }
      # fit Kalman
      if (kalm_Q_hat_mode == 'from VAR'){
        HHt_work = Q_hat
      } else if (kalm_Q_hat_mode == 'identity'){
        HHt_work = diag(1, ncol(factor_work))
      } else {
        HHt_work = diag(1, ncol(factor_work))
        cat('##########################  please enter valid entry for kalm_Q_hat_mode - identity has been selected by default')
      }
      outp = capture.output(
        kalm_fit <- fkf(
          a0 = factor_work[1,],
          P0 = cor(factor_work),
          dt = matrix(0, nrow = ncol(factor_work)),
          ct = matrix(0, nrow = nrow(y_mat_work)),
          Tt = A_hat,
          Zt = C_hat_work,
          HHt = HHt_work,
          GGt = R_hat_work,
          yt = y_mat_work,
          check.input = T
        ),
        type = 'output')
      # check for convergence error in Kalman
      if (sum(kalm_fit$status) > 0){
        var_to_remove_add = unique(rownames(C_hat_work)[kalm_fit$status])
        var_to_remove_all = c(var_to_remove_all, match(var_to_remove_add, y_df$ref))
      } else {
        stop = 1
      }
    } # while
    if(length(var_to_remove_all) > 0){
      var_to_remove_lab = paste0(y_df$ref[var_to_remove_all], collapse = ' | ')
    }
    if (length(var_to_remove_all) != length(var_to_remove_const)){
      cat('\n\n     *** removed for Kalman convergence:')
      cat(paste0('\n         ',y_df$ref[setdiff(var_to_remove_all, var_to_remove_const)]))
    }
    cat('\n  Done in', round((as.numeric(Sys.time())-as.numeric(start_time)) / 60), 'mins \n')
  }
  
  if (sum(kalm_fit$status) > 0){
    cat('###########  Kalman filter did not converge')
    err_code = paste0(country_variab_order[kalm_fit$status], collapse = ' | ')
    kalm_filter_factors = matrix(NA, nrow = nrow(factor_work), ncol = ncol(factor_work))
    kalm_logLik = kalm_AICc = kalm_observation_RMSE = kalm_observation_MAPE = kalm_RSS = kalm_TSS = kalm_Explain_var = 
      kalm_Explain_var_95 = kalm_Explain_var_99 = kalm_factors_test_RMSE = kalm_factors_test_MAPE = kalm_observation_MAPE_max = 
      kalm_observation_MAPE_95 = kalm_observation_MAPE_99 = kalm_RSS_95 = kalm_RSS_99 = kalm_TSS_95 = kalm_TSS_99 = 
      kalm_factors_test_MAPE_95 = kalm_factors_test_MAPE_99 = kalm_factors_test_MAPE_max = NA
  } else {
    err_code = ''
    kalm_filter_factors = t(kalm_fit$att)
    k = 2 * ncol(factor_work) + ncol(factor_work) * nrow(y_mat_work) + nrow(y_mat_work) ^ 2  # number of parameters, i.e. size of T,Z,HH,GG
    kalm_logLik = kalm_fit$logLik  # if NA means error in BLAS and LAPACK - see Details of fkf
    kalm_AICc = (2 * k - 2 * kalm_logLik) + (2 * k^2 + 2 * k) / (ncol(y_mat_work) - k - 1) # corrected AIC
    kalm_filter_observ_error = kalm_fit$vt
    kalm_RSS_val = kalm_filter_observ_error ^ 2
    kalm_TSS_val = (y_mat_work - mean(y_mat_work)) ^ 2
    kalm_RSS = sum(kalm_RSS_val)
    kalm_TSS = sum(kalm_TSS_val)
    ind_95 = kalm_RSS_val <= quantile(kalm_RSS_val, 0.95)
    kalm_RSS_95 = sum(kalm_RSS_val[ind_95])
    kalm_TSS_95 = sum(kalm_TSS_val[ind_95])
    ind_99 = kalm_RSS_val <= quantile(kalm_RSS_val, 0.99)
    kalm_RSS_99 = sum(kalm_RSS_val[ind_99])
    kalm_TSS_99 = sum(kalm_TSS_val[ind_99])
    kalm_Explain_var = 1 - kalm_RSS / kalm_TSS
    kalm_Explain_var_99 = 1 - kalm_RSS_99 / kalm_TSS_99
    kalm_Explain_var_95 = 1 - kalm_RSS_95 / kalm_TSS_95
    kalm_observation_RMSE = sqrt(mean(kalm_filter_observ_error ^ 2))
    kalm_MAPE_list = kalm_filter_observ_error %>%
      `rownames<-`(rownames(y_mat_work)) %>%
      `colnames<-`(colnames(y_mat_work)) %>%
      as.table() %>%
      as.data.frame(stringsAsFactors = F) %>%
      setNames(c('variable', 'year', 'residual')) %>%
      left_join(y_mat_work %>%
                  as.table() %>%
                  as.data.frame(stringsAsFactors = F) %>%
                  setNames(c('variable', 'year', 'original')),
                by = c("variable", "year")) %>%
      filter(abs(original) > .Machine$double.eps) %>%
      mutate(MAPE = abs(residual / original)) %>%
      arrange(desc(MAPE))
    kalm_observation_MAPE = mean(kalm_MAPE_list$MAPE)
    kalm_observation_MAPE_max = kalm_MAPE_list$MAPE[1]
    kalm_observation_MAPE_95 = mean((kalm_MAPE_list %>% filter(row_number() > ceiling(0.05 * nrow(kalm_MAPE_list))))$MAPE)
    kalm_observation_MAPE_99 = mean((kalm_MAPE_list %>% filter(row_number() > ceiling(0.01 * nrow(kalm_MAPE_list))))$MAPE)
    res_DFM_MAPE = res_DFM_MAPE %>% bind_rows(
      kalm_MAPE_list %>%
        filter(row_number() < ceiling(0.05 * nrow(kalm_MAPE_list))) %>%
        mutate(method = recov_met,
               data = df_type,
               Total_Factors = n_factor,
               algo = 'kalm') %>%
        select(algo, Total_Factors, everything())
    )
    if (!is.null(true_factors)){
      kalm_factors_test_RMSE = sqrt(mean((kalm_filter_factors - true_factors) ^ 2))
      kalm_factors_test_MAPE_val = true_factors %>%
        `rownames<-`(1:nrow(true_factors)) %>%
        as.table() %>%
        as.data.frame(stringsAsFactors = F) %>%
        setNames(c('year', 'variable', 'original')) %>%
        left_join(
          kalm_filter_factors %>%
            `rownames<-`(1:nrow(true_factors)) %>%
            `colnames<-`(colnames(true_factors)) %>%
            as.table() %>%
            as.data.frame(stringsAsFactors = F) %>%
            setNames(c('year', 'variable', 'prediction')),
          by = c("year", "variable")
        ) %>%
        mutate(residual = original - prediction) %>%
        filter(abs(original) > .Machine$double.eps) %>%
        mutate(MAPE = abs(residual / original)) %>%
        arrange(desc(MAPE))
      kalm_factors_test_MAPE = mean(kalm_factors_test_MAPE_val$MAPE)
      kalm_factors_test_MAPE_max = kalm_factors_test_MAPE_val$MAPE[1]
      kalm_factors_test_MAPE_95 = mean((kalm_factors_test_MAPE_val %>% filter(row_number() > ceiling(0.05 * nrow(kalm_factors_test_MAPE_val))))$MAPE)
      kalm_factors_test_MAPE_99 = mean((kalm_factors_test_MAPE_val %>% filter(row_number() > ceiling(0.01 * nrow(kalm_factors_test_MAPE_val))))$MAPE)
    } else {
      kalm_factors_test_RMSE = kalm_factors_test_MAPE = kalm_factors_test_MAPE_max = kalm_factors_test_MAPE_95 = kalm_factors_test_MAPE_99 = NA
    }
  }
  
  # save results
  res_DFM_stats = res_DFM_stats %>% bind_rows(
    data.frame(DFM = 'DFM_multivar',
               method = recov_met,
               data = df_type,
               Total_Factors = n_factor,
               country = 'All',
               AICc = kalm_AICc,
               LogLik = kalm_logLik,
               y_RSS = kalm_RSS,
               y_RSS_95 = kalm_RSS_95,
               y_RSS_99 = kalm_RSS_99,
               y_TSS = kalm_TSS,
               y_TSS_95 = kalm_TSS_95,
               y_TSS_99 = kalm_TSS_99,
               Explain_var = kalm_Explain_var,
               Explain_var_99 = kalm_Explain_var_99,
               Explain_var_95 = kalm_Explain_var_95,
               algo_err_code = err_code,
               removed_variable = ifelse(length(var_to_remove_all) > 0, toString(length(var_to_remove_all)), ''),
               kalm_Q_hat_mode = kalm_Q_hat_mode,
               kalm_y_RMSE = kalm_observation_RMSE,
               kalm_y_MAPE_max = kalm_observation_MAPE_max,
               kalm_y_MAPE = kalm_observation_MAPE,
               kalm_y_MAPE_99 = kalm_observation_MAPE_99,
               kalm_y_MAPE_95 = kalm_observation_MAPE_95,
               kalm_x_test_RMSE = kalm_factors_test_RMSE,
               kalm_x_test_MAPE_max = kalm_factors_test_MAPE_max,
               kalm_x_test_MAPE = kalm_factors_test_MAPE,
               kalm_x_test_MAPE_99 = kalm_factors_test_MAPE_99,
               kalm_x_test_MAPE_95 = kalm_factors_test_MAPE_95,              
               VAR_alpha = VAR_alpha,
               VAR_A_sparseness = A_sparseness,
               VAR_RMSE = VAR_RMSE,
               VAR_MAPE_max = VAR_MAPE_max,
               VAR_MAPE = VAR_MAPE,
               VAR_MAPE_99 = VAR_MAPE_99,
               VAR_MAPE_95 = VAR_MAPE_95,
               VAR_CV_RMSE = VAR_CV_RMSE,
               VAR_CV_RMSE_SD = VAR_CV_RMSE_SD,
               stringsAsFactors = F)
  )
  if (is.null(true_factors)){
    res_DFM_stats = res_DFM_stats %>% select(-starts_with('kalm_x_test'))
  } else {
    res_DFM_stats = res_DFM_stats %>% select(AICc, LogLik, starts_with('kalm_x_test'), everything())
  }
  
  res_DFM_list[[df_type]][[recov_met]][['DFM_multivar']][[paste0(n_factor, '_factors')]][['All']] = list(
    factors = kalm_filter_factors %>% `rownames<-`(year_order) %>% `colnames<-`(country_fact_order),
    A_hat = A_hat %>% `rownames<-`(country_fact_order) %>% `colnames<-`(country_fact_order),
    C_hat = C_hat %>% `rownames<-`(country_variab_order) %>% `colnames<-`(country_fact_order),
    C_hat_remov = C_hat_remov,
    Q_hat = Q_hat %>% `rownames<-`(country_fact_order) %>% `colnames<-`(country_fact_order),
    R_hat = R_hat %>% `rownames<-`(country_variab_order) %>% `colnames<-`(country_variab_order),
    R_hat_remov = R_hat_remov,
    VAR_fit = VAR_fit,
    kalm_fit = kalm_fit,
    kalm_y_mat_work = y_mat_work
  )
  
  res_DFM_factors = res_DFM_factors %>% bind_rows(
    kalm_filter_factors %>%
      as.data.frame() %>%
      setNames(country_fact_order) %>%
      mutate(year = year_order) %>%
      gather('Country_Factor', 'val', -year) %>%
      rowwise() %>%
      mutate(country = strsplit(Country_Factor, '_')[[1]][1],
             Factor = strsplit(Country_Factor, '_')[[1]][2]) %>%
      select(-Country_Factor) %>%
      mutate(Factor = as.numeric(Factor),
             method = recov_met,
             data = df_type,
             Total_Factors = n_factor,
             Var_removed = var_to_remove_lab,
             DFM = 'DFM_multivar') %>%
      select(DFM, country, Total_Factors, Factor, everything())
  )
  
  res_DFM_loadings = res_DFM_loadings %>% bind_rows(
    loading_work %>% mutate(DFM = 'DFM_multivar')
  )
  
  return(list(
    res_DFM_stats = res_DFM_stats,
    res_DFM_list = res_DFM_list,
    res_DFM_factors = res_DFM_factors,
    res_DFM_loadings = res_DFM_loadings,
    res_DFM_MAPE = res_DFM_MAPE
  ))
}

# simulates multi-country DFM to choose best VAR_alpha and kalm_Q_hat_mode
simulate_DFM_multivar = function(res_DFM_simulation, reload_list, seed_i, rep, n_var, n_country, n_factor, n_year, A_sp, alpha, kalm_Q,
                                 univ_reload, multiv_reload_VAR, multiv_reload_Kalm){
  
  country_code_sim = data.frame(country = paste('country', formatC(c(1:n_country), width=3, flag='0'), sep = ''),
                                ID = 1:n_country, stringsAsFactors = F)
  variable_code_sim = data.frame(variable = paste('variable', formatC(c(1:n_var), width=2, flag='0'), sep = ''),
                                 CODE = as.character(1:n_var), stringsAsFactors = F)
  country_variable_lab = (expand.grid(country_code_sim$country, paste0('_', variable_code_sim$variable)) %>%
                            mutate(COL = paste0(Var1, Var2)) %>%
                            select(COL) %>%
                            arrange(COL))$COL
  country_fact_lab = (expand.grid(country_code_sim$country, paste0('_', 1:n_factor)) %>%
                        mutate(COL = paste0(Var1, Var2)) %>%
                        select(COL) %>%
                        arrange(COL))$COL

  # reload for simulation only
  if (!is.null(reload_list)){
    res_DFM_list_reload_sim = reload_list$res_DFM_list
  } else {
    res_DFM_list_reload_sim = NULL
  }
  
  # generate matrices
  set.seed(seed_i)
  A_gen = rsparsematrix(n_country * n_factor, n_country * n_factor, density = A_sp) %>% as.matrix()
  
  list_C_gen = list()
  set.seed(seed_i)
  seed_list_C = sample(1:1e4, n_country)
  for (lc in 1:n_country){
    set.seed(seed_list_C[lc])
    list_C_gen[[lc]] = matrix(runif(n_var * n_factor), nrow = n_var, ncol = n_factor)
  }
  C_gen = bdiag(list_C_gen) %>% as.matrix()
  
  set.seed(seed_i)
  R_gen = diag(runif(nrow(C_gen)))
  
  Q_gen = diag(1, nrow(A_gen))
  
  set.seed(seed_i)
  x0 = runif(nrow(A_gen))
  P0 = diag(1, nrow(A_gen))
  
  # generate observations yt and states (factors) xt
  kalman.filter=dse::SS(F = A_gen,
                        Q = Q_gen,
                        H = C_gen,
                        R = R_gen,
                        z0 = x0,
                        P0 = P0
  )
  
  simulate.kalman.filter=simulate(kalman.filter, start = 1, freq = 1, sampleT = n_year)
  yt = simulate.kalman.filter$output  # prime colonne sono variabili per prima nazione, ecc
  xt = simulate.kalman.filter$state
  
  df_SP_orig_sim = yt %>%
    as.data.frame() %>%
    mutate(year = c(1:nrow(yt))) %>%
    gather('country_var', 'val', -year) %>%
    mutate(country_var = rep(country_variable_lab, each = n_year)) %>%
    rowwise() %>%
    mutate(country = strsplit(country_var, '_')[[1]][1],
           variable = strsplit(country_var, '_')[[1]][2]) %>%
    select(-country_var) %>%
    spread(variable, val) %>%
    arrange(country, year)
  
  df_SP_sim = df_SP_orig_sim %>%
    left_join(country_code_sim, by = "country") %>%
    rename(TIME = year) %>%
    select(-country) %>%
    select(ID, everything()) %>%
    setNames(c('ID', 'TIME', 1:nrow(variable_code_sim))) %>%
    arrange(ID, TIME)
  
  true_factors_sim = xt %>%
    # as.data.frame() %>%
    `colnames<-`(country_fact_lab) %>%
    t() %>%
    MARSS::zscore() %>%
    t()
  
  # run algo
  res_DFM_factors_sim = res_DFM_loadings_sim = res_DFM_stats_sim = c()
  res_DFM_list_sim = list()
  for (country_i in 1:n_country){
    
    # Dynamic Factor Model - evaluate for single country (data are standardized)
    DFM_uni_sim = evaluate_DFM_univar(res_DFM_factors = res_DFM_factors_sim,
                                      res_DFM_loadings = res_DFM_loadings_sim,
                                      res_DFM_stats = res_DFM_stats_sim,
                                      res_DFM_list = res_DFM_list_sim,
                                      res_DFM_list_reload = res_DFM_list_reload_sim,
                                      variable_code = variable_code_sim,
                                      country_code = country_code_sim,
                                      df_SP = df_SP_sim,
                                      n_factor, country_i, df_type = 'type0', recov_met = 'met0',
                                      univ_reload = univ_reload,
                                      dfm_eval = 'BFGS',
                                      dfm_max_iter = 2000)
    res_DFM_stats_sim = DFM_uni_sim$res_DFM_stats
    res_DFM_list_sim = DFM_uni_sim$res_DFM_list
    res_DFM_factors_sim = DFM_uni_sim$res_DFM_factors
    res_DFM_loadings_sim = DFM_uni_sim$res_DFM_loadings
  } # country_i
  
  # adjust DFM taking into account all countries together (factors are standardized)
  DFM_multi_sim = evaluate_DFM_multivar(res_DFM_factors = res_DFM_factors_sim,
                                        res_DFM_loadings = res_DFM_loadings_sim,
                                        res_DFM_stats = res_DFM_stats_sim,
                                        res_DFM_list = res_DFM_list_sim,
                                        df_SP_orig = df_SP_orig_sim,
                                        res_DFM_MAPE = c(),
                                        df_type = 'type0', recov_met = 'met0', n_factor = n_factor,
                                        VAR_alpha = alpha,
                                        kalm_Q_hat_mode = kalm_Q,
                                        true_factors = true_factors_sim,
                                        res_DFM_list_reload = res_DFM_list_reload_sim,
                                        multiv_reload_VAR = multiv_reload_VAR,
                                        multiv_reload_Kalm = multiv_reload_Kalm)
  res_DFM_stats_sim = DFM_multi_sim$res_DFM_stats
  
  res_DFM_simulation = res_DFM_simulation %>% bind_rows(
    res_DFM_stats_sim %>%
      filter(DFM == 'DFM_multivar') %>%
      select(-DFM, -method, -data, -Total_Factors, -country) %>%
      mutate(Trial = rep,
             N_variable = n_var,
             N_country = n_country,
             N_factor = n_factor,
             N_year = n_year,
             Sparseness_A = A_sp,
             VAR_alpha = alpha,
             Q_hat_mode = kalm_Q)
  )
  
  return(list(
    res_DFM_simulation = res_DFM_simulation,
    DFM_uni_sim = DFM_uni_sim,
    DFM_multi_sim = DFM_multi_sim
  ))
}

# evaluate index for DFM
evaluate_index_DFM = function(res_index_DFM, res_DFM_best_model, res_DFM_factors, res_DFM_list, res_DFM_loadings, index_1_thresh, index_2_thresh,
                              load_thresh, leading_var, leading_sign, df_type, recov_met, dfm_met, expl_var_to_show=0){
  
  # load best model stats
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
              loading = mean(loading), .groups = 'drop') %>%
    ungroup()
  load_range = c(min(load_ref$load_min), max(load_ref$load_max))
  
  scores_out_raw = scores_out_index = data.frame(COUNTRY = c('INDEX', sort(unique(factors_ref$country))), stringsAsFactors = F)
  load_out = data.frame(VARIABLE = c('VARIABLE', unique(load_ref$variable)), stringsAsFactors = F)
  
  row_list = list()
  i = 1
  for (yr in sort(unique(factors_ref$year))){
    
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
    
    for (fact in 1:min(c(max_factor, 2))){
      
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
      
      row_list[[toString((fact - 1) * min(c(max_factor, 2)) + 1)]][[i]] = ggplotGrob(p)
      row_list[[toString(min(c(max_factor, 2)) * fact)]][[i]] = ggplotGrob(p_load)
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
  png(paste0('./Results/1_Index_plot_thrs', 0, '_', df_type, '_', recov_met, '_', dfm_met, '.png'), width = 22.5 * max_factor, height = 7.5 * uniqueN(factors_ref$year), units = 'in', res=300) # height=60 per serie originale
  grid.draw(g)
  dev.off()
  
  colnames(scores_out_raw) = gsub('\\.1', '', colnames(scores_out_raw))
  write.table(scores_out_raw, paste0('./Results/1_Index_trhs', load_thresh, '_', df_type, '_', recov_met, '_', dfm_met, '_scores_raw.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  colnames(scores_out_index) = gsub('\\.1', '', colnames(scores_out_index))
  write.table(scores_out_index, paste0('./Results/1_Index_trhs', load_thresh, '_', df_type, '_', recov_met, '_', dfm_met, '_scores_index.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  write.table(load_out, paste0('./Results/1_Index_trhs', load_thresh, '_', df_type, '_', recov_met, '_', dfm_met, '_loadings.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  res_index_DFM[[df_type]][[recov_met]][[dfm_met]][['scores_raw']] = scores_out_raw
  res_index_DFM[[df_type]][[recov_met]][[dfm_met]][['scores_index']] = scores_out_index
  res_index_DFM[[df_type]][[recov_met]][[dfm_met]][['loadings']] = load_out
  res_index_DFM[[df_type]][[recov_met]][[dfm_met]][['Expl_Var']] = best_model %>% select(Total_Factors, starts_with('Expl'))
  
  return(res_index_DFM)
}

# space label for z-axis in 3D plot
space_label = function(label_list, offset){
  
  label_list$space = 0
  while (min(label_list$space[-1]) < offset){
    label_list = label_list %>%
      arrange(val) %>%
      mutate(space = c(0, diff(val))) %>%
      mutate(new_val = val)
  for (i in 2:nrow(label_list)){
    if (label_list$space[[i]] <= offset){
      half_space = (offset - label_list$space[[i]]) / 2
      label_list$new_val[[i - 1]] = label_list$val[[i - 1]] - half_space
      label_list$new_val[[i]] = label_list$val[[i]] + half_space
    }
  }
  label_list = label_list %>%
    mutate(val = new_val)
  offset = offset * 0.99
  }
  
  label_list = label_list %>% 
    select(val, label)
  
  return(label_list)
}

# nested facets label
facet_nested <- function(rows = NULL, cols = NULL, scales = "fixed", space = "fixed",
                         shrink = TRUE, labeller = "label_value", as.table = TRUE,
                         switch = NULL, drop = TRUE, margins = FALSE, facets = NULL,
                         nest_line = FALSE, resect = unit(0, "mm"), bleed = FALSE)
{
  if (!is.null(facets)) {
    rows <- facets
  }
  if (is.logical(cols)) {
    margins <- cols
    cols <- NULL
  }
  scales <- match.arg(scales, c("fixed", "free_x", "free_y", "free"))
  free <- list(x = any(scales %in% c("free_x", "free")),
               y = any(scales %in% c("free_y", "free")))
  
  space <- match.arg(space, c("fixed","free_x","free_y","free"))
  space_free <- list(x = any(space %in% c("free_x", "free")),
                     y = any(space %in% c("free_y", "free")))
  
  if (!is.null(switch) && !switch %in% c("both","x","y")) {
    stop("switch must be either 'both', 'x', or 'y'", call. = FALSE)
  }
  
  facets_list <- ggplot2:::grid_as_facets_list(rows, cols)
  n <- length(facets_list)
  if (n > 2L) {
    stop("A grid facet specification can't have more than two dimensions",
         .call = FALSE)
  }
  if (n == 1L) {
    rows <- quos()
    cols <- facets_list[[1]]
  } else {
    rows <- facets_list[[1]]
    cols <- facets_list[[2]]
  }
  labeller <- ggplot2:::check_labeller(labeller)
  ggproto(NULL, FacetNested, shrink = shrink,
          params = list(
            rows = rows,
            cols = cols,
            margins = margins,
            free = free,
            space_free = space_free,
            labeller = labeller,
            as.table = as.table,
            switch = switch,
            drop = drop,
            nest_line = nest_line,
            resect = resect,
            bleed = bleed
          ))
}
FacetNested <- ggplot2::ggproto(
  "FacetNested", ggplot2::
    FacetGrid,
  map_data = function(data, layout, params) {
    # Handle empty data
    if (ggplot2:::empty(data)) {
      return(cbind(data, PANEL = integer(0)))
    }
    # Setup variables
    rows <- params$rows
    cols <- params$cols
    
    vars <- c(names(rows), names(cols))
    margin_vars <- list(intersect(names(rows), names(data)),
                        intersect(names(cols), names(data)))
    
    # Add variables
    data <- reshape2::add_margins(data, margin_vars, params$margins)
    facet_vals <- ggplot2:::eval_facets(c(rows, cols), data, params$plot$env)
    
    # Only set as missing if it has no variable in that direction
    missing_facets <- character(0)
    if (!any(names(rows) %in% names(facet_vals))){
      missing_facets <- c(missing_facets, setdiff(names(rows), names(facet_vals)))
    }
    if (!any(names(cols) %in% names(facet_vals))){
      missing_facets <- c(missing_facets, setdiff(names(cols), names(facet_vals)))
    }
    
    # Fill in missing values
    if (length(missing_facets) > 0) {
      to_add <- unique(layout[missing_facets])
      data_rep <- rep.int(1:nrow(data), nrow(to_add))
      facet_rep <- rep(1:nrow(to_add), each = nrow(data))
      data <- plyr::unrowname(data[data_rep, , drop = FALSE])
      facet_vals <- plyr::unrowname(
        cbind(facet_vals[data_rep, , drop = FALSE],
              to_add[facet_rep, , drop = FALSE])
      )
    }
    
    # Match columns to facets
    if (nrow(facet_vals) == 0) {
      data$PANEL <- NO_PANEL
    } else {
      facet_vals[] <- lapply(facet_vals[], as.factor)
      facet_vals[] <- lapply(facet_vals[], addNA, ifany = TRUE)
      keys <- plyr::join.keys(facet_vals, layout, by = vars[vars %in% names(facet_vals)])
      data$PANEL <- layout$PANEL[match(keys$x, keys$y)]
    }
    data
  },
  compute_layout = function(data, params)
  {
    rows <- params$rows
    cols <- params$cols
    dups <- intersect(names(rows), names(cols))
    if (length(dups) > 0) {
      stop("Facetting variables can only appear in row or cols, not both.\n",
           "Problems: ", paste0(dups, collapse = "'"), call. = FALSE)
    }
    base_rows <- combine_nested_vars(data, params$plot_env, rows, drop = params$drop)
    if (!params$as.table) {
      rev_order <- function(x) factor(x, levels = rev(ggplot2:::ulevels(x)))
    }
    base_cols <- combine_nested_vars(data, params$plot_env, cols, drop = params$drop)
    base <- ggplot2:::df.grid(base_rows, base_cols)
    base <- reshape2::add_margins(base, list(names(rows), names(cols)), params$margins)
    base <- unique(base)
    panel <- plyr::id(base, drop = TRUE)
    panel <- factor(panel, levels = seq_len(attr(panel, "n")))
    rows <- if (!length(names(rows))) {
      1L
    } else {
      plyr::id(base[names(rows)], drop = TRUE)
    }
    cols <- if (!length(names(cols))) {
      1L
    } else {
      plyr::id(base[names(cols)], drop = TRUE)
    }
    panels <- data.frame(PANEL = panel, ROW = rows, COL = cols,
                         base, check.names = FALSE, stringsAsFactors = FALSE)
    panels <- panels[order(panels$PANEL), , drop = FALSE]
    rownames(panels) <- NULL
    panels$SCALE_X <- if (params$free$x) {
      panels$COL
    } else {
      1L
    }
    panels$SCALE_Y <- if (params$free$y) {
      panels$ROW
    } else {
      1L
    }
    panels
  },
  draw_panels = function(panels, layout, x_scales, y_scales, ranges, coord,
                         data, theme, params)
  {
    panel_table <- FacetGrid$draw_panels(panels, layout, x_scales, y_scales,
                                         ranges, coord, data, theme, params)
    
    # Setup strips
    col_vars  <- unique(layout[names(params$cols)])
    row_vars  <- unique(layout[names(params$rows)])
    attr(col_vars, "type")  <- "cols"
    attr(col_vars, "facet") <- "grid"
    attr(row_vars, "type")  <- "rows"
    attr(row_vars, "facet") <- "grid"
    
    # Build strips
    strips <- render_strips(col_vars, row_vars, params$labeller, theme)
    switch_x <- !is.null(params$switch) && params$switch %in% c("both", "x")
    switch_y <- !is.null(params$switch) && params$switch %in% c("both", "y")
    
    # Merging strips
    merge_cols <- apply(col_vars, 2, function(x) any(rle(x)$lengths > 1))
    merge_rows <- apply(row_vars, 2, function(x) any(rle(x)$lengths > 1))
    
    if (any(merge_cols)) {
      if (switch_x) {
        panel_table <- merge_strips(panel_table, strips$x$bottom,
                                    col_vars, switch_x, params, theme, "x")
      } else {
        panel_table <- merge_strips(panel_table, strips$x$top,
                                    col_vars, switch_x, params, theme, "x")
      }
    }
    
    if (any(merge_rows)) {
      if (switch_y) {
        panel_table <- merge_strips(panel_table, strips$y$left,
                                    row_vars, switch_y, params, theme, "y")
      } else {
        panel_table <- merge_strips(panel_table, strips$y$right,
                                    row_vars, switch_y, params, theme, "y")
      }
    }
    panel_table
  }
)

# Helper functions -----------------------------------------------

combine_nested_vars <- function(data, env = emptyenv(), vars = NULL, drop = TRUE) {
  if (length(vars) == 0)
    return(data.frame())
  values <- ggplot2:::compact(plyr::llply(data, ggplot2:::eval_facets, facets = vars,
                                          env = env))
  has_all <- unlist(lapply(values, length)) == length(vars)
  if (!any(has_all)) {
    missing <- lapply(values, function(x) setdiff(names(vars), names(x)))
    missing_txt <- vapply(missing, var_list, character(1))
    name <- c("Plot", paste0("Layer ", seq_len(length(data) - 1)))
    stop("At least one lyaer must contain all faceting variables: ",
         var_list(names(vars)), ".\n", paste0("* ", name, " is missing",
                                              missing_txt, collapse = "\n"),
         call. = FALSE)
  }
  base <- unique(plyr::ldply(values[has_all]))
  if (!drop) {
    base <- ggplot2:::unique_combs(base)
  }
  for (value in values[!has_all]) {
    if (ggplot2:::empty(value))
      next
    old <- base[setdiff(names(base), names(value))]
    new <- unique(value[intersect(names(base), names(value))])
    if (drop) {
      new <- ggplot2:::unique_combs(new)
    }
    old[setdiff(names(base), names(value))] <- rep("", nrow(old))
    base <- rbind(base, ggplot2:::df.grid(old, new))
  }
  if (ggplot2:::empty(base)) {
    stop("Facetting variables must have at least one value",
         call. = FALSE)
  }
  base
}


merge_strips <- function(panel_table, strip, vars, switch, params, theme, orient = c("x","y"))
{
  if (is.null(strip)) {
    return(panel_table)
  }
  n_levels <- nrow(strip[[1]]$layout)
  splitstrip <- lapply(seq_len(n_levels), function(i) {
    switch(orient,
           x = lapply(strip, function(x) x[i, ]),
           y = lapply(strip, function(x) x[, i]))
    
  })
  
  if (params$bleed) {
    merge <- apply(vars, 2, function(x) any(rle(x)$lengths > 1))
  } else {
    merge <- sapply(1:ncol(vars), function(i){
      x <- apply(subset.data.frame(vars, select = seq(i)), 1, paste0, collapse = "")
      return(any(rle(x)$lengths > 1))
    })
  }
  
  if (orient == "y" && !switch) {
    vars <- rev(vars)
    merge <- rev(merge)
  }
  if (orient == "x" && switch) {
    vars <- rev(vars)
    merge <- rev(merge)
    splitstrip <- rev(splitstrip)
  }
  
  sizes <- switch(orient,
                  x = do.call(unit.c, lapply(splitstrip, max_height)),
                  y = do.call(unit.c, lapply(splitstrip, max_width)))
  
  assign("panel_table", panel_table, 1)
  
  grabwhat <- switch(orient,
                     x = grepl("strip-t|strip-b", panel_table$layout$name),
                     y = grepl("strip-r|strip-l", panel_table$layout$name))
  
  pos_y <- unique(panel_table$layout$t[grabwhat])
  pos_x <- unique(panel_table$layout$l[grabwhat])
  panel_pos <- find_panel(panel_table)
  
  if (orient == "x") {
    nudge <- if (pos_y < panel_pos$t) -1 else -1
    panel_table <- panel_table[-pos_y,]
    panel_table <- gtable_add_rows(panel_table, sizes, pos_y + nudge)
    
  } else {
    nudge <- if (pos_x < panel_pos$l) -1 else 0
    panel_table <- panel_table[, -pos_x]
    panel_table <- gtable_add_cols(panel_table, sizes, pos_x + nudge)
  }
  
  for(i in seq_len(n_levels)) {
    if (!merge[i]) {
      panel_table <- gtable_add_grob(
        panel_table, splitstrip[[i]],
        t = pos_y + switch(orient, x = i + nudge, y = 0),
        l = pos_x + switch(orient, x = 0, y = i + nudge),
        z = 2, clip = "on",
        name = paste0("strip-", orient, "-", seq_along(splitstrip[[i]]))
      )
    } else {
      j <- as.numeric(as.factor(vars[,i]))
      ends <- cumsum(rle(j)$lengths)
      starts <- c(1, which(diff(j) != 0) + 1)
      panel_table <- gtable_add_grob(
        panel_table, splitstrip[[i]][starts],
        t = switch(orient, x = pos_y + i + nudge, y = pos_y[starts]),
        b = switch(orient, x = pos_y + i + nudge, y = pos_y[ends]),
        l = switch(orient, x = pos_x[starts], y = pos_x + i + nudge),
        r = switch(orient, x = pos_x[ends],   y = pos_x + i + nudge),
        z = 2, clip = "on",
        name = paste0("strip-", orient, "-", seq_along(splitstrip[[i]][starts]))
      )
      
      if(params$nest_line && any(starts != ends)) {
        insert_here <- which(starts != ends)
        indicator <- linesGrob(
          x = switch(orient,
                     x = unit(c(0, 1), "npc") + c(1, -1) * params$resect,
                     y = if (switch) c(1, 1) else c(0, 0)),
          y = switch(orient,
                     x = if (switch) c(1, 1) else c(0, 0),
                     y = unit(c(0, 1), "npc") + c(1, -1) * params$resect),
          gp = grid::gpar(col = theme$line$colour,
                          lty = theme$line$linetype,
                          lwd = theme$line$size * .pt,
                          lineend = theme$line$lineend))
        panel_table <- gtable_add_grob(
          panel_table, lapply(seq_along(insert_here), function(x) indicator),
          t = switch(orient, x = pos_y + i + nudge,
                     y = pos_y[starts[insert_here]]),
          b = switch(orient, x = pos_y + i + nudge,
                     y = pos_y[ends[insert_here]]),
          l = switch(orient, x = pos_x[starts[insert_here]],
                     y = pos_x + i + nudge),
          r = switch(orient, x = pos_x[ends[insert_here]],
                     y = pos_x + i + nudge),
          z = 3, clip = "on",
          name = "nesting-indicator"
        )
      }
    }
  }
  panel_table
}

# merge rectangle
merge_rectagle = function(d_rect_t){
  
  comb = d_rect_t %>% group_by(algo, factor) %>% summarise(COUNT = n(), .groups = 'drop')
  d_rect = c()
  for (i in 1:nrow(comb)){
    if (comb$COUNT[i] == 1){
      d_rect = d_rect %>% bind_rows(
        d_rect_t %>% filter(algo == comb$algo[i]) %>% filter(factor == comb$factor[i])
      )
    } else {
      tt = d_rect_t %>%
        filter(algo == comb$algo[i]) %>%
        filter(factor == comb$factor[i]) %>%
        mutate(merge = 0)
      cc = 1
      for (j in 2:nrow(tt)){
        if (tt$year[j] - tt$year[j - 1] == 1){
          tt$merge[j - 1] = cc
          tt$merge[j] = cc
        } else {
          cc = cc + 1
        }
      }
      tt_out = tt %>%
        filter(merge == 0) %>%
        bind_rows(tt_out = tt %>%
                    filter(merge != 0) %>%
                    group_by(algo, factor, merge) %>%
                    summarise(year = min(year),
                              xmin = min(xmin),
                              xmax = max(xmax),
                              ymin = min(ymin),
                              ymax = max(ymax)), .groups = 'drop')
      d_rect = d_rect %>% bind_rows(
        tt_out
      )
      
    }
  }
  if ('merge' %in% colnames(d_rect)){d_rect = d_rect %>% select(-merge)}
  
  return(d_rect)
}