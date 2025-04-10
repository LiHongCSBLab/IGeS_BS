## elastic net hypermeter tuning

xgboost_tuning <- function(traindata){
  
  set.seed(123)
  
  mod_xgboost = makeLearner("classif.xgboost", 
                            predict.type = "prob", 
                            fix.factors.prediction = TRUE,
                            print_every_n = 1000L,
                            nthread = 1L,
                            eval_metric='logloss',                          
                            #save_name = "/dev/null",
                            config = list(on.learner.error = "warn"))
  
  # mod_xgboost = makeUndersampleWrapper(mod_xgboost,  usw.rate = 1/(N_num/P_num))

  xgboost_ps <- ParamHelpers::makeParamSet(
    ParamHelpers::makeDiscreteParam("booster",
                                    values = c("gbtree")),
    ParamHelpers::makeNumericParam("eta",
                                   lower = 0,
                                   upper = 1,
                                   trafo = function(x) x/1e4,
                                   requires = quote(booster %in% c("gbtree"))),
    ParamHelpers::makeNumericParam("gamma",
                                   lower = 0,
                                   upper = 10,
                                   requires = quote(booster %in% c("gbtree"))),
    ParamHelpers::makeIntegerParam("max_depth",
                                   lower = 1L,
                                   upper = 14L,
                                   requires = quote(booster %in% c("gbtree"))),
    ParamHelpers::makeNumericParam("min_child_weight",
                                   lower = 0L,
                                   upper = 7L,
                                   trafo = function(x) 2^{x},
                                   requires = quote(booster %in% c("gbtree"))),
    ParamHelpers::makeNumericParam("subsample",
                                   lower = 0,
                                   upper = 1,
                                   requires = quote(booster %in% c("gbtree"))),
    ParamHelpers::makeNumericParam("colsample_bytree",
                                   lower = 0,
                                   upper = 1,
                                   requires = quote(booster %in% c("gbtree"))),
    ParamHelpers::makeNumericParam("colsample_bylevel",
                                   lower = 0,
                                   upper = 1,
                                   requires = quote(booster %in% c("gbtree"))),
    ParamHelpers::makeNumericParam("lambda",
                                   lower = -10,
                                   upper = 10,
                                   trafo = function(x) 2^{x},
                                   requires = quote(booster %in% c("gblinear"))),
    ParamHelpers::makeNumericParam("lambda_bias",
                                   lower = 0,
                                   upper = 1e5,
                                   trafo = function(x) x/1e4,
                                   requires = quote(booster %in% c("gblinear"))),
    ParamHelpers::makeNumericParam("alpha",
                                   lower = 0,
                                   upper = 1e5,
                                   trafo = function(x) x/1e4,
                                   requires = quote(booster %in% c("gblinear"))),
    ParamHelpers::makeNumericParam("base_score",
                                   lower = 0,
                                   upper = 1),
    ParamHelpers::makeIntegerParam("nrounds",
                                   lower = 1L,
                                   upper = 5000L)
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

  rdesc = makeResampleDesc("CV", stratify = T, iters = 10L)


  
  xgboost_tuning = tuneParams(mod_xgboost, 
                              task = traindata, 
                              resampling = rdesc,
                              par.set = xgboost_ps, 
                              control = ctrl,
                              measures = list(auc, setAggregation(auc, test.sd)),
                              show.info = TRUE)
  
  return(xgboost_tuning)
}