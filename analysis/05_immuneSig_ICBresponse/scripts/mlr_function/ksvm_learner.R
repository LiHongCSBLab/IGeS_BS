
## c-SVM hypermeter tuning

csvm_tuning = function(traindata){
  
  set.seed(123)
  mod_ksvm = makeLearner("classif.ksvm", 
                          predict.type = "prob", 
                          fix.factors.prediction = TRUE)

  # mod_ksvm = makeUndersampleWrapper(mod_ksvm,  usw.rate = 1/(N_num/P_num))

  # set model with parameters that need to be tuned
  
  # Define SVM learner parameter set ----
  svm_ps <- ParamHelpers::makeParamSet(
    # ParamHelpers::makeDiscreteParam("type",
    #                                 values = c("C-svc", "C-bsvc")),
    ParamHelpers::makeNumericParam("C",
                                   lower = -5,
                                   upper = 5,
                                   trafo = function(x) 10^{x} # ,
                                   #requires = quote(type %in% c("C-svc", "C-bsvc"))
                                   ),
    
    ParamHelpers::makeDiscreteParam("kernel", values = c("rbfdot",
                                                         "vanilladot",
                                                         # "polydot",
                                                         "laplacedot",
                                                         "besseldot")),
    ParamHelpers::makeNumericParam("sigma",
                                   lower = -5,
                                   upper = 5,
                                   trafo = function(x) 10^{x},
                                   requires = quote(kernel %in% c("rbfdot",
                                                                  "laplacedot",
                                                                  "besseldot"))),
    ParamHelpers::makeIntegerParam("degree",
                                   lower = 1L,
                                   upper = 11L,
                                   requires = quote(kernel %in% c("polydot",
                                                                  "besseldot"))),
    ParamHelpers::makeNumericParam("scale",
                                   lower = -5,
                                   upper = 5,
                                   trafo = function(x) 10^{x},
                                   requires = quote(kernel %in% c("polydot"))),
    ParamHelpers::makeNumericParam("offset",
                                   lower = -3,
                                   upper = 3,
                                   trafo = function(x) 2^{x},
                                   requires = quote(kernel %in% c("polydot"))),
    ParamHelpers::makeIntegerParam("order",
                                   lower = 0L,
                                   upper = 6L,
                                   requires = quote(kernel == "besseldot"))
  )
  
  

  mbo_ctrl <- mlrMBO::makeMBOControl(impute.y.fun = function(x, y, opt.path, ...) -0.5) # This is the worst AUC
  mbo_ctrl <- mlrMBO::setMBOControlTermination(mbo_ctrl, 
                                               iters = 10 )
  surrogate_lrn <- mlr::makeImputeWrapper(mlr::makeLearner("regr.ranger", 
                                          predict.type = "se", 
                                          replace = FALSE),
                                          classes = list(numeric = mlr::imputeMax(2),
                                                       factor = mlr::imputeConstant("__miss__")))
  ctrl <- mlr:::makeTuneControlMBO(learner = surrogate_lrn,
                                 mbo.control = mbo_ctrl)

  # ctrl = makeTuneControlGrid()

  # rdesc = makeResampleDesc("CV", iters = 10L)
  rdesc = makeResampleDesc("CV", stratify = T, iters = 10L)

  # parameters tuning
   svm_tuning = tuneParams(mod_ksvm, 
                    task = traindata, 
                    resampling = rdesc,
                    par.set = svm_ps, 
                    control = ctrl,
                    measures = list(auc, setAggregation(auc, test.sd)),
                    show.info = TRUE)
  
  
  return(svm_tuning)
}


## nu-SVM hypermeter tuning

nusvm_tuning = function(traindata){
  
  set.seed(123)
  mod_ksvm = makeLearner("classif.ksvm", predict.type = "prob", fix.factors.prediction = TRUE)
  mod_ksvm = makeUndersampleWrapper(mod_ksvm,  usw.rate = 1/(N_num/P_num))

  # set model with parameters that need to be tuned
  
  # Define SVM learner parameter set ----
  svm_ps <- ParamHelpers::makeParamSet(

    ParamHelpers::makeDiscreteParam("type",
                                    values = c("nu-svc")),

    ParamHelpers::makeNumericParam("nu",
                                 lower = 0,
                                 upper = 1,
                                 requires = quote(type == "nu-svc")),,
    
    ParamHelpers::makeDiscreteParam("kernel", values = c("rbfdot",
                                                         "vanilladot",
                                                         "polydot",
                                                         "laplacedot",
                                                         "besseldot")),
    ParamHelpers::makeNumericParam("sigma",
                                   lower = -5,
                                   upper = 2,
                                   trafo = function(x) 10^{x},
                                   requires = quote(kernel %in% c("rbfdot",
                                                                  "laplacedot",
                                                                  "besseldot"))),
    ParamHelpers::makeIntegerParam("degree",
                                   lower = 1L,
                                   upper = 5L,
                                   requires = quote(kernel %in% c("polydot",
                                                                  "besseldot"))),
    ParamHelpers::makeNumericParam("scale",
                                   lower = -5,
                                   upper = 5,
                                   trafo = function(x) 10^{x},
                                   requires = quote(kernel %in% c("polydot"))),
    ParamHelpers::makeNumericParam("offset",
                                   lower = -3,
                                   upper = 3,
                                   trafo = function(x) 2^{x},
                                   requires = quote(kernel %in% c("polydot"))),
    ParamHelpers::makeIntegerParam("order",
                                   lower = 0L,
                                   upper = 6L,
                                   requires = quote(kernel == "besseldot"))
  )
  
  

  mbo_ctrl <- mlrMBO::makeMBOControl(impute.y.fun = function(x, y, opt.path, ...) -0.5) # This is the worst AUC
  mbo_ctrl <- mlrMBO::setMBOControlTermination(mbo_ctrl, 
                                               iters = 10 )
  surrogate_lrn <- mlr::makeImputeWrapper(mlr::makeLearner("regr.ranger", 
                                          predict.type = "se", 
                                          replace = FALSE),
                                          classes = list(numeric = mlr::imputeMax(2),
                                                       factor = mlr::imputeConstant("__miss__")))
  ctrl <- mlr:::makeTuneControlMBO(learner = surrogate_lrn,
                                 mbo.control = mbo_ctrl)

  # ctrl = makeTuneControlGrid()

  rdesc = makeResampleDesc("CV", iters = 10L)

  # parameters tuning
  svm_tuning = tuneParams(mod_ksvm, 
                    task = traindata, 
                    resampling = rdesc,
                    par.set = svm_ps, 
                    control = ctrl,
                    measures = list(auc, setAggregation(auc, test.sd)),
                    show.info = TRUE)
  
  return(svm_tuning)
}
