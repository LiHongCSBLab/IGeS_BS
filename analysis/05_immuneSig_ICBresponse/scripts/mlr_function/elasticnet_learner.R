## elastic net hypermeter tuning

elasticnet_tuning <- function(traindata){
  
  set.seed(1234)
  
  mod_elasticnet = makeLearner("classif.cvglmnet", 
                               predict.type = "prob", 
                               fix.factors.prediction = TRUE,
                               nlambda = 1000L,
                               lambda.min.ratio = 1e-5,
                               nfolds = 5,
                               config = list(on.learner.error = "warn"))
  
  # mod_elasticnet = makeUndersampleWrapper(mod_elasticnet,  
  #                              usw.rate = 1/(N_num/P_num))

  elasticnet_ps <- ParamHelpers::makeParamSet(
    ParamHelpers::makeLogicalParam("standardize"),
    ParamHelpers::makeDiscreteParam("s",
                                    values = c("lambda.1se", "lambda.min")),
    ParamHelpers::makeNumericParam("alpha",
                                   lower = 0,
                                   upper = 1)
    )

  
  # mbo_ctrl <- mlrMBO::makeMBOControl(impute.y.fun = function(x, y, opt.path, ...) -0.5) # This is the worst AUC
  # mbo_ctrl <- mlrMBO::setMBOControlTermination(mbo_ctrl, 
  #                                              iters = 10 )
  # surrogate_lrn <- mlr::makeImputeWrapper(mlr::makeLearner("regr.ranger", 
  #                                         predict.type = "se", 
  #                                         replace = FALSE),
  #                                         classes = list(numeric = mlr::imputeMax(2),
  #                                                      factor = mlr::imputeConstant("__miss__")))
  # ctrl <- mlr:::makeTuneControlMBO(learner = surrogate_lrn,
  #                                mbo.control = mbo_ctrl)

  ctrl = makeTuneControlGrid()
  rdesc = makeResampleDesc("CV", stratify = T, iters = 10L)


  elasticnet_tuning = tuneParams(mod_elasticnet, 
                   task = traindata, 
                   resampling = rdesc,
                   par.set = elasticnet_ps, 
                   control = ctrl,
                   measures = list(auc, setAggregation(auc, test.sd)),
                   show.info = TRUE)
 

  
  return(elasticnet_tuning)
}