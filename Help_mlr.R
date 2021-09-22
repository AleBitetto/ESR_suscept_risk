# models settings, tuned parameters and parameters' space
model_settings = function(task, flag_tuning, save_all, algo_type_work, tuning_strategy, tuning_criterion, regressor_lab, reload_out, inn_cross_val_fold){
  
  lambda_opt = list()
  
  if (flag_tuning == F & length(reload_out) > 0 & !is.null(reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$settings)){
    settings = reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$settings
    cat('\n           *** parameters reloaded')
  } else {
    
    if (is.null(reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$settings$params)){
      cat('\n           *** parameters set to default (all 1) - ready for tuning')
      reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$status = 'default'
      
      ### Random Forest
      {
        if (algo_type_work == "RandomForest"){
          learner_type = "regr.randomForest"
          
          # settings
          params = list(replace = F,
                        importance = T,
                        localImp = F,
                        proximity = T,
                        do.trace = F,
                        keep.forest = T)
          
          # tuned parameters
          params = c(params,
                     ntree = 1,
                     mtry = 1,
                     nodesize = 1)
          
          
          # parameters space
          n_var = ncol(task$env$data) - 1
          n_obs = nrow(task$env$data)
          if (tuning_strategy == 'grid'){
          param_set = makeParamSet(
            makeDiscreteParam("ntree", c(50, 100, 200)),
            makeDiscreteParam("mtry", unique(round(c(0.3, 0.5, 0.7, 1) * n_var))),
            makeDiscreteParam("nodesize", c(1, 5, 10, 20, 40))
          )
          } else if (tuning_strategy == 'bayes'){
            param_set = makeParamSet(
              makeIntegerParam("ntree", 10, 1000),
              makeIntegerParam("mtry", 1, n_var),
              makeIntegerParam("nodesize", 1, n_obs - 1)
            )
          }
          if (n_var == 1){
            param_set$pars$mtry = NULL
            params$mtry = 1
          }
        }
      }
      
      ### Conditional Inference Trees
      {
        if (algo_type_work == "CITree"){
          learner_type = "regr.ctree"
          
          # settings
          n_var = ncol(task$env$data) - 1
          params = list(mincriterion = 0.95,
                        teststat = "quad",
                        stump = F,
                        maxsurrogate = 0, # for missing values
                        nresample = 300,
                        mtry = n_var,
                        maxdepth = 2000)
          
          # tuned parameters
          params = c(params,
                     testtype = "Bonferroni",
                     minsplit = 1,
                     minbucket = 1)
          
          
          # parameters space
          n_obs = nrow(task$env$data)
          if (tuning_strategy == 'grid'){
            param_set = makeParamSet(
              makeDiscreteParam("testtype", c("Bonferroni", "MonteCarlo", "Univariate", "Teststatistic")),
              makeDiscreteParam("minsplit", c(1, 5, 10, 20, 40)),
              makeDiscreteParam("minbucket", c(1, 5, 10, 20, 40))
            )
          } else if (tuning_strategy == 'bayes'){
            param_set = makeParamSet(
              makeDiscreteParam("testtype", c("Bonferroni", "MonteCarlo", "Univariate", "Teststatistic")),
              makeIntegerParam("minsplit", 1, 100),
              makeIntegerParam("minbucket", 1, n_obs - 1)
            )
          }
        }
      }
      
      ### GBM
      {
        if (algo_type_work == "GBM"){
          learner_type = "regr.gbm"
          
          # settings
          params = list(distribution = "gaussian",
                        bag.fraction = 1,
                        keep.data = T)
          
          # tuned parameters
          params = c(params,
                     n.trees = 1,
                     shrinkage = 1,
                     n.minobsinnode = 1,
                     interaction.depth = 1)
          
          # parameters space
          n_var = ncol(task$env$data) - 1
          n_obs = nrow(task$env$data)
          if (tuning_strategy == 'grid'){
          param_set = makeParamSet(
            makeDiscreteParam("n.trees", c(20, 100, 200)),
            makeDiscreteParam("shrinkage", c(0.5, 0.1, 0.01)),
            makeDiscreteParam("n.minobsinnode", c(1, 5)),
            makeDiscreteParam("interaction.depth", sort(unique(c(1, 5, 10, floor(sqrt(n_var))))))
          )
          } else if (tuning_strategy == 'bayes'){
            param_set = makeParamSet(
              makeIntegerParam("n.trees", 10, 1000),
              makeNumericParam("shrinkage", 0.001, 1),
              makeIntegerParam("n.minobsinnode", 1, round(n_obs / 10)),
              makeIntegerParam("interaction.depth", 1, 5)
            )
          }
        }
      }
      
      ### Elastic-Net
      {
        if (algo_type_work == "ElasticNet"){
          learner_type = "regr.cvglmnet"
          
          # settings
          lambda_opt = 0.5
          params = list(alpha = 1,
                        nfolds = max(c(3, inn_cross_val_fold)),
                        standardize = T,   # already standardized in model.feature
                        s = "lambda.1se",  # keep this fixed because the code substitute lambda.1se with chosen lambda.opt
                        type.measure = "mse",
                        type.gaussian = "covariance")
          
          # tuned parameters
          params = c(params)
          
          # parameters space
          n_var = ncol(task$env$data) - 1
          n_obs = nrow(task$env$data)
          if (tuning_strategy == 'grid'){
            param.set = makeParamSet(
              makeDiscreteParam("alpha", values = seq(0, 1, by = 0.1)),
              makeDiscreteParam("type.gaussian", values = c("covariance", "naive"))
            )
          } else if (tuning_strategy == 'bayes'){
            param_set = makeParamSet(
              makeNumericParam("alpha", 0, 1),
              makeDiscreteParam("type.gaussian", values = c("covariance", "naive"))
            )
          }
        }
      }
      
      ### GLM regression
      {
        if (algo_type_work == "GLM"){
          learner_type = "regr.glm"
          
          # settings
          n_var = ncol(task$env$data) - 1
          params = list(family = "gaussian",
                        start = rep(1, n_var + 1))   # initial coefficients
          
          # tuned parameters
          params = c(params,
                     gaussian.link = "identity",
                     epsilon = 1e-8)
          
          # parameters space
          if (tuning_strategy == 'grid'){
            param.set = makeParamSet(
              makeDiscreteParam("epsilon", c(1e-8, 1.1e-8)),
              makeDiscreteParam("gaussian.link", c("identity"))#, "log", "inverse"))
            )
          } else if (tuning_strategy == 'bayes'){
            param_set = makeParamSet(
              makeDiscreteParam("epsilon", c(1e-8, 1.1e-8)),  # useless, used just to allow the search
              makeDiscreteParam("gaussian.link", c("identity"))#, "log", "inverse"))
            )
          }
        }
      }
      
      ### SVM - polynomial
      {
        if (algo_type_work == "SVM-Poly"){
          learner_type = "regr.svm"
          
          # settings
          params = list(type = "eps-regression",
                        kernel = "polynomial",
                        fitted = T,
                        scale = T)
          
          # tuned parameters
          params = c(params,
                     cost = 1,
                     gamma = 1,
                     coef0 = 0,
                     degree = 1)
          
          # parameters space
          n_var = ncol(task$env$data) - 1
          n_obs = nrow(task$env$data)
          data_variance = getTaskData(task) %>% select_if(is.numeric) %>% unlist() %>% sd(na.rm = T)
          if (tuning_strategy == 'grid'){
            param.set = makeParamSet(
              makeDiscreteParam("cost", c(0.0001, 0.001, 0.01, 0.1, 1, 10, 25, 50, 100)),
              makeDiscreteParam("gamma", c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 1 / n_var), 1 / (n_var * data_variance)),
              # makeDiscreteParam("coef0", c(-1, -0.1, -0.01, 0, 0.01, 0.1, 1)),
              makeDiscreteParam("degree", c(2, 3))
            )
          } else if (tuning_strategy == 'bayes'){
            param_set = makeParamSet(
              makeNumericParam("cost", 0.0001, 100),
              makeNumericParam("gamma", 0.0001, 10),
              # makeNumericParam("coef0", -1, 1),
              makeIntegerParam("degree", 2, 3)
            )
          }
        }
      }
      
      ### SVM - RBF
      {
        if (algo_type_work == "SVM-RBF"){
          learner_type = "regr.svm"
          
          # settings
          params = list(type = "eps-regression",
                        kernel = "radial",
                        fitted = T,
                        scale = T)
          
          # tuned parameters
          params = c(params,
                     cost = 1,
                     gamma = 1)
          
          # parameters space
          n_var = ncol(task$env$data) - 1
          n_obs = nrow(task$env$data)
          data_variance = getTaskData(task) %>% select_if(is.numeric) %>% unlist() %>% sd(na.rm = T)
          if (tuning_strategy == 'grid'){
            param.set = makeParamSet(
              makeDiscreteParam("cost", c(0.0001, 0.001, 0.01, 0.1, 1, 10, 25, 50, 100)),
              makeDiscreteParam("gamma", c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 1 / n_var), 1 / (n_var * data_variance))
            )
          } else if (tuning_strategy == 'bayes'){
            param_set = makeParamSet(
              makeNumericParam("cost", 0.0001, 100),
              makeNumericParam("gamma", 0.0001, 10)
            )
          }
        }
      }
      
      ### Single layer neural network
      {
        if (algo_type_work == "SingleNN"){
          learner_type = "regr.nnet"
          
          # settings
          params = list(maxit = 200,
                        skip = T,
                        MaxNWts = 5000,
                        trace = F)
          
          # tuned parameters
          params = c(params,
                     size = 2,
                     # skip = T,
                     decay = 0)
          
          # parameters space
          n_var = ncol(task$env$data) - 1
          if (tuning_strategy == 'grid'){
            param.set = makeParamSet(
              makeDiscreteParam("size", values = seq(1, round(n_var * 0.7), by = round(n_var / 20))),
              # makeDiscreteParam("skip", c(TRUE, FALSE)),
              makeDiscreteParam("decay", values = seq(0, 1, by = 0.01))
            )
          } else if (tuning_strategy == 'bayes'){
            param_set = makeParamSet(
              makeIntegerParam("size", 1, ifelse(n_var > 1, round(n_var *0.8), 1)),
              # makeDiscreteParam("skip", c(TRUE, FALSE)),
              makeNumericParam("decay", 0, 10)
            )
          }
        }
      }
      
      ### Fully Connected Neural Network
      {
        if (algo_type_work == "FCNN"){
          learner_type = "regr.h2o.deeplearning"
          
          # settings
          params = list(seed = 66,
                        epochs = 30,
                        adaptive_rate = T,
                        stopping_rounds = 5,
                        stopping_metric = "AUTO",
                        fast_mode = T)   #quiet_mode
          
          # tuned parameters
          params = c(params,
                     list(loss = "Quadratic",
                          hidden = c(100, 100)))
          
          # parameters space
          n_var = ncol(task$env$data) - 1
          n_obs = nrow(task$env$data)
          if (tuning_strategy == 'grid'){
            param.set = makeParamSet(
            )
          } else if (tuning_strategy == 'bayes'){
            param_set = makeParamSet(
              makeDiscreteParam("loss", c("Quadratic", "Automatic")),
              makeIntegerVectorParam("hidden", len = 2, lower = 1, upper = n_var)#, trafo = function(x) x[x != 0])
            )
          }
        }
      }
      
      ### Spline regression
      {
        if (algo_type_work == "Spline"){
          learner_type = "regr.crs"
          
          # settings
          n_var = ncol(task$env$data) - 1
          params = list(random.seed = 66,
                        segments = rep(1, n_var),
                        complexity = "degree-knots")
          
          # tuned parameters
          params = c(params,
                     list(degree = rep(2, n_var)))
          
          # parameters space
          
          if (tuning_strategy == 'grid'){
            param.set = makeParamSet(
            )
          } else if (tuning_strategy == 'bayes'){
            param_set = makeParamSet(
              makeIntegerVectorParam("degree", len = n_var, lower = 1, upper = 2)
            )
          }
        }
      }
      
      ### MARS Multivariate Adaptive Regression Splines
      {
        if (algo_type_work == "MARS"){
          learner_type = "regr.earth"
          
          # settings
          n_var = ncol(task$env$data) - 1
          params = list(trace = 0,
                        nk = min(200, max(20, 2 * n_var)) + 1,
                        nfold = 0)
          
          # tuned parameters
          params = c(params,
                     list(degree = 2))
          
          # parameters space
          
          if (tuning_strategy == 'grid'){
            param.set = makeParamSet(makeIntegerParam("degree", 1, 3)
            )
          } else if (tuning_strategy == 'bayes'){
            param_set = makeParamSet(
              makeIntegerParam("degree", 1, 3)
            )
          }
        }
      }
      
    } else {
      cat('\n           *** parameters reloaded:', reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$status)
      reload_set = reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$settings
      learner_type = reload_set$learner_type
      params = reload_set$params
      param_set = reload_set$param_set
    } # default parameters
    
    settings = list(learner_type = learner_type,
               params = params,
               param_set = param_set,
               lambda_opt = lambda_opt)
    
    if (save_all){
      cat('\n           *** exporting parameters:', reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$status)
      reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$settings = settings
    }
  } # flag_tuning
  
  return(list(settings = settings,
              reload_out = reload_out))
}

# tuning function
model_tuning = function(algo_type_work, tuning_strategy, learner, task, param_set, out_cross_val_fold, inn_cross_val_fold, tuning_criterion){
  
  start_time = Sys.time()
  par_names = names(param_set$pars)
  par_type = unlist(lapply(learner$par.vals, typeof))[par_names]
  
  # define optimization strategy
  if (tuning_strategy == 'grid'){
    control = makeTuneControlGrid()
  } else if (tuning_strategy == 'bayes'){
    control = makeTuneControlMBO()
    # control$mbo.control$iters = 2
    # control = makeTuneControlIrace(budget = 200)
  }
  
  # define INNER and OUTER resampling strategy
  inner = makeResampleDesc(method = "CV",
                           iters = inn_cross_val_fold,
                           predict = "both"
                           # stratify.cols = c('country')
  )
  
  outer = makeResampleDesc(method = "CV",
                           iters = out_cross_val_fold,
                           predict = "both"
                           # stratify.cols = c('country')
  )
  
  # define performance measures
  eval(parse(text=paste0(c("meas = list(", tuning_criterion, ", setAggregation(", tuning_criterion, ", test.sd), setAggregation(", tuning_criterion,
                           ", train.", tuning_criterion, "), setAggregation(", tuning_criterion, ", train.sd))"), collapse = "")))
  
  # make learner wrapper
  learner_t = makeTuneWrapper(learner = learner,
                              resampling = inner,
                              par.set = param_set,
                              control = control,
                              measures = meas,
                              show.info = F)
  
  # tuning
  set.seed(1, "L'Ecuyer-CMRG")
  parallelStart(mode = "socket", cpus = detectCores(), show.info = F)
  parallel::clusterSetRNGStream(iseed = 1)
  res = resample(learner = learner_t,
                 task = task,
                 resampling = outer,
                 extract = getTuneResult,
                 show.info = F,
                 models = T)
  parallelStop()
  
  # get tuning results
  res_summ = getNestedTuneResultsOptPathDf(res)
  if (algo_type_work == "FCNN"){    # todo: per i parametri vettoriali bisognerebbe fare un controllo a parte, anche per la parte sotto dove applichi il tipo
    res_summ = res_summ %>%         #       e anche nell'applicazione dei parametri nell'addestrare i modelli
      unite("hidden", starts_with("hidden"), remove = T, sep = ",")
  }
  colnames(res_summ) = gsub(paste0(tuning_criterion, "."), "", colnames(res_summ))
  res_summ_avg = res_summ %>%
    rename_at(vars(contains(tuning_criterion)), funs(sub(tuning_criterion, 'mean', .))) %>%  # added for rmse
    group_by_(.dots = par_names) %>%
    summarise_at(.vars = c("test.mean", "train.mean", "test.sd", "train.sd"),
                 .funs = funs(mean, min, max)) %>%
    ungroup() %>%
    mutate(DELTA_TRAIN_TEST = abs(train.mean_mean - test.mean_mean),
           FULL_SET_PERF = 0) %>%
    select_(.dots = c(par_names, "FULL_SET_PERF", "DELTA_TRAIN_TEST", "test.mean_mean",
                      "train.mean_mean", "test.sd_mean", "test.mean_min", "test.mean_max",
                      "train.sd_mean", "train.mean_min", "train.mean_max")) %>%
    mutate_at(.vars = par_names, .funs = as.character)
  if (algo_type_work != "FCNN"){
    for (p in names(par_type)){
      res_summ_avg = res_summ_avg %>% mutate_at(.vars = p, .funs = funs(!!paste0("as.", par_type[p])))
    }
  }
  
  # add performance on full set
  if (algo_type_work == "ElasticNet"){
    res_summ_avg_t = c()
  }
  if (algo_type_work == "FCNN"){
    params_t = learner$par.vals
    
    # hidden_model_list = res_summ %>%    # needed to reload trained model ?
    #   group_by(iter) %>%
    #   arrange(dob) %>%
    #   filter(row_number() == n()) %>%
    #   group_by(hidden) %>%
    #   filter(row_number() == 1)
    # res_summ_avg = res_summ_avg %>%
    #   left_join(hidden_model_list %>% select(hidden, iter), by = "hidden") %>%
    #   filter(!is.na(iter))
    for (i in 1:nrow(res_summ_avg)){
      hidden_vector = strsplit(res_summ_avg[i, 'hidden'] %>% unlist(), ',')[[1]] %>% as.numeric() %>% as.vector()
      hidden_vector = hidden_vector[hidden_vector != 0]
      params_t[['hidden']] = hidden_vector
      params_t['quiet_mode'] = T
      if (length(params_t[['hidden']]) > 0 ){
        learner_t = makeLearner(cl = learner$id,
                                par.vals = params_t,
                                predict.type = 'response',
                                fix.factors.prediction = T)
        oo = capture.output(model_t <- train(learner_t, task), silent = T)
        oo = capture.output(prediction_t <- predict(model_t, task = task), silent = T)
        res_summ_avg$FULL_SET_PERF[i] = performance(prediction_t, eval(parse(text=tuning_criterion)))
      } else {
        res_summ_avg$FULL_SET_PERF[i] = NA
      }
    }
    res_summ_avg = res_summ_avg %>%
      filter(!is.na(FULL_SET_PERF))
  } else {
    for (i in c(1:nrow(res_summ_avg))){
      params_t = learner$par.vals
      
      params_t[names(param_set$pars)] = as.list(res_summ_avg[i, par_names])
      
      learner_t = makeLearner(cl = learner$id,
                              par.vals = params_t,
                              predict.type = 'response',
                              fix.factors.prediction = T)
      
      model_t = train(learner_t, task)
      if (algo_type_work != "ElasticNet"){
        prediction_t = predict(model_t, task = task)
        res_summ_avg$FULL_SET_PERF[i] = performance(prediction_t, eval(parse(text=tuning_criterion)))
      } else {
        cv_res = model_t$learner.model$cvm   # get cv-performance metric for each lamba tested
        for (ll in 1:length(cv_res)){
          lam = model_t$learner.model$lambda[ll]   # get corresponding lambdas
          model_t$learner.model$lambda.1se = lam
          prediction_t = predict(model_t, task = task)
          res_summ_avg_t = res_summ_avg_t %>%
            bind_rows(
              res_summ_avg[i, ] %>%
                mutate(lambda = lam,
                       DELTA_TRAIN_TEST = -1,
                       FULL_SET_PERF = as.numeric(performance(prediction_t, eval(parse(text=tuning_criterion)))),
                       "test.mean_mean" = cv_res[i],
                       "train.mean_mean" = cv_res[i],
                       "test.sd_mean" = -1,
                       "test.mean_min" = -1,
                       "test.mean_max" = -1,
                       "train.sd_mean" = -1,
                       "train.mean_min" = -1,
                       "train.mean_max" = -1)
            )
        }
      }
    }
  }
  if (algo_type_work == "ElasticNet"){
    res_summ_avg = res_summ_avg_t %>%
      select(alpha, type.gaussian, lambda, everything())
  }
  
  # order according to tuning_criterion optimal value
  eval(parse(text = paste0("optim_crit = ", tuning_criterion, "$minimize")))
  if (optim_crit){
    res_summ_avg = res_summ_avg %>% arrange(FULL_SET_PERF)
  } else {
    res_summ_avg = res_summ_avg %>% arrange(desc(FULL_SET_PERF))
  }
  
  if (sum(is.na(res_summ$error.message)) < nrow(res_summ)){
    cat("\n\n           ########## Error on", algo_type_work, "tuning, check tune_res$res_summ\n\n")
  }
  
  # collapse parameters name and values in single string
  res_summ_log = res_summ_avg
  if (algo_type_work == "ElasticNet"){par_names = c(par_names, 'lambda')}
  for (p in par_names){
    res_summ_log = res_summ_log %>% mutate_at(.vars = p, .funs = funs(paste0(p, "=", .)))
  }
  out_cv_test_size = max(unlist(lapply(getResamplingIndices(res)$test.inds, length)))
  inn_cv_test_size = suppressWarnings(max(unlist(lapply(getResamplingIndices(res, inner = T)[[1]]$test.inds, length))))
  res_summ_log = cbind(data.frame(ALGO = algo_type_work, PERF = tuning_criterion, N_OBS = task$task.desc$size,
                                  OUTER_CV = out_cross_val_fold, OUTER_TEST_SIZE = out_cv_test_size,
                                  INNER_CV = inn_cross_val_fold, INNER_TEST_SIZE = inn_cv_test_size, stringsAsFactors = F),
                       res_summ_log %>%
                         unite_("PARAMETERS", par_names, sep = ",") %>%
                         select(PARAMETERS),
                       res_summ_avg %>% select_(.dots = paste0("-", par_names)))
  
  cat(' Done in', round((as.numeric(Sys.time())-as.numeric(start_time)) / 60), 'mins ')
  
  return(list(res = res,
              res_summ = res_summ,
              res_summ_avg = res_summ_avg,
              res_summ_log = res_summ_log))
  
}

# performance of final model with resampling
model_perf = function(algo_type_work, model, learner, task, prediction, out_cross_val_fold, inn_cross_val_fold, tuning_criterion, param_set){
  
  # define resampling strategy
  rdesc = makeResampleDesc(method = "CV",
                           iters = out_cross_val_fold,
                           predict = "both"
  )
  
  # define performance measures
  eval(parse(text=paste0(c("meas = list(", tuning_criterion, ", setAggregation(", tuning_criterion, ", test.sd), setAggregation(", tuning_criterion,
                           ", train.mean), setAggregation(", tuning_criterion, ", train.sd))"), collapse = "")))
  
  # resample
  if (algo_type_work != "ElasticNet"){
    set.seed(1, "L'Ecuyer-CMRG")
    if (algo_type_work == "FCNN"){
      oo = capture.output(perf <- resample(learner = learner,
                                           task = task,
                                           resampling = rdesc,
                                           measures = meas,
                                           show.info = F,
                                           models = T), silent = T)
    } else {
      perf = resample(learner = learner,
                      task = task,
                      resampling = rdesc,
                      measures = meas,
                      show.info = F,
                      models = T)
    }
    
    train_val = perf$measures.train[[tuning_criterion]]
    test_val = perf$measures.test[[tuning_criterion]]
    stat = as.data.frame(t(perf$aggr))
    
    if (algo_type_work == "FCNN"){
      best_par = paste0('loss=', learner$par.vals['loss'], ',hidden=', paste0(learner$par.vals[['hidden']], collapse = ","))
    } else {
      best_par = unlist(learner$par.vals[names(param_set$pars)])
      best_par = paste0(paste(names(best_par), best_par, sep = "="), collapse = ",")
    }
  } else {
    best_lambda = model$learner.model$lambda.1se
    best_lambda_ind = which.min(abs(model$learner.model$lambda - best_lambda))
    train_val = test_val = model$learner.model$cvm[best_lambda_ind]
    best_par = unlist(learner$par.vals[names(param_set$pars)])
    best_par = paste0(paste0(names(best_par), "=", best_par, collapse = ","), ',lambda=', best_lambda)
    stat = data.frame(test.sd = -1, train.sd = -1)
  }
  
  
  res_perf = data.frame(ALGO = algo_type_work, PERF = tuning_criterion, 
                        N_OBS = task$task.desc$size, OUTER_CV = out_cross_val_fold, INNER_CV = inn_cross_val_fold,
                        PARAMETERS = best_par,
                        FULL_SET_PERF = performance(prediction, eval(parse(text=tuning_criterion))),
                        DELTA_TRAIN_TEST = abs(mean(test_val) - mean(train_val)),
                        TEST_mean = mean(test_val),
                        TRAIN_mean = mean(train_val),
                        TEST_sd = stat %>% select(ends_with("test.sd")) %>% as.numeric,
                        TEST_min = min(test_val),
                        TEST_max = max(test_val),
                        TRAIN_sd = stat %>% select(ends_with("train.sd")) %>% as.numeric,
                        TRAIN_min = min(train_val),
                        TRAIN_max = max(train_val),
                        stringsAsFactors = F)
  
  return(res_perf)
}

# wrapper for model estimation
threshold_sensitivity_fit = function(flag_tuning, force_tuning, reload_final_perf, save_all, inn_cross_val_fold, out_cross_val_fold, tuning_criterion, tuning_strategy,
                                     df_type, recov_met, fit_met, algo_type, algo_type_work, var_target, reload_out, index_1_set, index_2_set,
                                     res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
                                     regr_to_test, index_1_thresh = NULL, index_2_thresh = NULL,
                                     df_work_orig = NULL, df_work_index = NULL, df_work_rank = NULL, df_work_raw = NULL){
  
  # select dataset
  if (regr_to_test == 'original'){
    df_work = df_work_orig
    regressor_lab = 'x_original'
    index_1_thresh_out = index_2_thresh_out = 'None'
  } else if (regr_to_test == 'index'){
    df_work = df_work_index
    regressor_lab = ifelse(is.null(index_2_thresh), paste0('x_Ind1_', index_1_thresh), paste0('x_Ind1_', index_1_thresh, '-Ind2_', index_2_thresh))
    index_1_thresh_out = as.character(index_1_thresh)
    index_2_thresh_out = as.character(index_2_thresh)
    if (length(index_2_thresh_out) == 0){index_2_thresh_out = 'None'}
  } else if (regr_to_test == 'rank_index'){
    df_work = df_work_rank
    regressor_lab = 'x_rank_index'
    index_1_thresh_out =  paste0(levels(cut(1, breaks = c(-Inf,index_1_set,Inf))), collapse = '')
    index_2_thresh_out =  paste0(levels(cut(1, breaks = c(-Inf,index_2_set,Inf))), collapse = '')
  } else if (regr_to_test == 'raw_index'){
    df_work = df_work_raw
    regressor_lab = 'x_raw_index'
    index_1_thresh_out = index_2_thresh_out = 'None'
  }
  suppressWarnings(parallelStop())

  # todo: mettere country e year? randomforest non supporta più di 53 classi
  # df_work = df_work %>%
  #   mutate(country = as.factor(country),
  #          year = as.integer(year))
  obs_lab = df_work %>% select(country, year)
  df_work = df_work %>%
    select(-country, -year) %>%
    mutate_all(function(x) { attributes(x) <- NULL; x })
  if (algo_type_work %in% c("SingleNN", "FCNN", "Spline")){
    df_work = df_work %>%
      clever_scale(c('TARGET')) %>%
      data.frame()
  }
  target_var_stats = paste0(df_work$TARGET %>% mean(na.rm = T) %>% round(4), '±', df_work$TARGET %>% sd(na.rm = T) %>% round(4))

  # define task for outer resampling
  task = makeRegrTask(id = "Regression",
                      data = df_work,
                      target = "TARGET",
                      fixup.data = "warn", 
                      check.data = T
  )
  
  # define model settings
  set_out = model_settings(task, flag_tuning, save_all, algo_type_work, tuning_strategy, tuning_criterion, regressor_lab, reload_out, inn_cross_val_fold)
  settings = set_out$settings
  reload_out = set_out$reload_out
  learner_type = settings$learner_type
  params = settings$params
  param_set = settings$param_set
  
  # set learner
  learner = makeLearner(cl = learner_type,
                        par.vals = params,
                        predict.type = 'response',
                        fix.factors.prediction = T)

  # tune parameters and set best parameters
  err_lab = ""
  suppressWarnings(rm(tune_res))
  if ((flag_tuning == T & reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$status != 'tuned') | force_tuning){
    cat('\n           --- Tuning...')
    if (algo_type_work == "MARS"){
      err = try(tune_res <- model_tuning(algo_type_work, tuning_strategy, learner, task, param_set, out_cross_val_fold, inn_cross_val_fold, tuning_criterion), silent = T)
      if (class(err) == "try-error" & regr_to_test != 'original'){    # take best tuning params (degree) from original variable set
        tune_res = reload_out[[algo_type_work]][['x_original']][[tuning_criterion]]$tuning
        cat('\n         ############## error in fitting. Used tuning from original variables ##############')
        err_lab = "Error in fitting. Used tuning from original variables"
      }
    } else {
      tune_res = model_tuning(algo_type_work, tuning_strategy, learner, task, param_set, out_cross_val_fold, inn_cross_val_fold, tuning_criterion)
    }
    reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$status = 'tuned'
    } else if (reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$status == 'tuned'){
    # reload tuned list
    tune_res = reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$tuning
    if (is.null(tune_res)){
      cat('\n           ############ tune_res not found')
    } else {
      cat('\n           --- Reloaded tuned list')
    }
  }
  
  # set best parameters
  tune_res_summ = tune_res$res_summ_avg
  if (algo_type_work == "FCNN"){
      hidden_vector = strsplit(tune_res_summ[1, 'hidden'] %>% unlist(), ',')[[1]] %>% as.numeric() %>% as.vector()
      hidden_vector = hidden_vector[hidden_vector != 0]
      params[['hidden']] = hidden_vector
      params['quiet_mode'] = T
      params['loss'] = tune_res_summ[1, 'loss'] %>% unlist()
  } else {
    params[names(param_set$pars)] = as.list(tune_res_summ[1, names(param_set$pars)])
  }

  learner = makeLearner(cl = learner_type,
                        par.vals = params,
                        predict.type = 'response',
                        fix.factors.prediction = T)
  
  # train model
  cat('\n           --- Training final model')
  reload_model = reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$model
  if (!is.null(reload_model) & reload_final_perf == T){
    model = reload_model$model
    cat('    - reloaded')
  } else {
    if (algo_type_work == "FCNN"){
      oo = capture.output(model <- train(learner, task), silent = T)
    } else {
      model = train(learner, task)
    }
  }
  if (algo_type_work != "ElasticNet"){
    if (algo_type_work == "FCNN"){
      oo = capture.output(prediction <- predict(model, task = task), silent = T)
    } else {
      prediction = predict(model, task = task)
    }
  } else {
    best_lambda = tune_res_summ$lambda[1]
    model$learner.model$lambda.1se = best_lambda
    prediction = predict(model, task = task)
  }
  
  # assess final performance - average only on outer cross-validation
  cat('\n           --- Assessing final performances')
  if (!is.null(reload_model) & reload_final_perf == T){
    perf_res = reload_model$performance
    cat('    - reloaded')
  } else {
    perf_res = model_perf(algo_type_work, model, learner, task, prediction, out_cross_val_fold, inn_cross_val_fold, tuning_criterion, param_set)
  }
  
  # store results
  res_thresh_sensitivity_list = res_thresh_sensitivity_list %>% bind_rows(
    cbind(data.frame(data = df_type,
          method = recov_met,
          fit_method = fit_met,
          target_var = var_target,
          target_var_avg = target_var_stats,
          regressor_type = regr_to_test,
          index_1_thresh = index_1_thresh_out,
          index_2_thresh = index_2_thresh_out,
          algo_main = algo_type,
          error = err_lab, stringsAsFactors = F),
          tune_res$res_summ_log %>%
            mutate(Best_params = ifelse(row_number()==1, 'YES', '')) %>%
            select(Best_params, everything()))
  )
  
  res_thresh_sensitivity_best = res_thresh_sensitivity_best %>% bind_rows(
    cbind(data.frame(data = df_type,
          method = recov_met,
          fit_method = fit_met,
          target_var = var_target,
          target_var_avg = target_var_stats,
          regressor_type = regr_to_test,
          index_1_thresh = index_1_thresh_out,
          index_2_thresh = index_2_thresh_out,
          algo_main = algo_type,
          error = err_lab, stringsAsFactors = F),
          perf_res)
  )
  
  res_thresh_sensitivity_residual = res_thresh_sensitivity_residual %>% bind_rows(
    cbind(data.frame(data = df_type,
          method = recov_met,
          fit_method = fit_met,
          target_var = var_target,
          regressor_type = regr_to_test,
          index_1_thresh = index_1_thresh_out,
          index_2_thresh = index_2_thresh_out,
          algo_main = algo_type,
          algo = algo_type_work,
          perform = tuning_criterion, stringsAsFactors = F),
          obs_lab %>% bind_cols(
            prediction$data %>%
              select(-id) %>%
              mutate(residual = truth - response)
          ))
  )
  
  # save results
  if (save_all){
    reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$tuning = tune_res
    settings$params = params
    reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$settings = settings
    reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$model = list(task = task,
                                                                                   learner = learner,
                                                                                   model = model,
                                                                                   prediction = prediction,
                                                                                   performance = perf_res)
    cat('\n           *** exporting results:', reload_out[[algo_type_work]][[regressor_lab]][[tuning_criterion]]$status)
  }
  
  return(list(res_thresh_sensitivity_list = res_thresh_sensitivity_list,
              res_thresh_sensitivity_best = res_thresh_sensitivity_best,
              res_thresh_sensitivity_residual = res_thresh_sensitivity_residual,
              reload_out = reload_out))
}
