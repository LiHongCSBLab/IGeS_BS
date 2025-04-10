## elastic net hypermeter tuning

naiveBayes_tuning <- function(traindata){
  
  set.seed(123)
  
  mod_naiveBayes = makeLearner("classif.naiveBayes", 
                               predict.type = "prob", 
                               fix.factors.prediction = TRUE)
  
  naiveBayes_ps <- ParamHelpers::makeParamSet(
  ParamHelpers::makeNumericParam("laplace",
                                 lower = 0,
                                 upper = 10)
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

  # rdesc = makeResampleDesc("CV", iters = 10L)
  rdesc = makeResampleDesc("CV", stratify = T, iters = 10L)


  naiveBayes_tuning = tuneParams(mod_naiveBayes, 
                   task = traindata, 
                   resampling = rdesc,
                   par.set = naiveBayes_ps, 
                   control = ctrl,
                   measures = list(auc, setAggregation(auc, test.sd)),
                   show.info = TRUE)
 
  return(naiveBayes_tuning)
}