## logistic regression hypermeter tuning

glm_tuning <- function(traindata){
  
  
  set.seed(123)
  
  mod_glm = makeLearner("classif.binomial", 
                        predict.type = "prob", 
                        fix.factors.prediction = TRUE)
  
  mod_glm = makeUndersampleWrapper(mod_glm,  usw.rate = 1/(N_num/P_num))

  glm_ps <- ParamHelpers::makeParamSet(
             ParamHelpers::makeDiscreteParam("link",
                                  values = c("logit",
                                             "probit",
                                             "cloglog"))
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

  glm_tuning = tuneParams(mod_glm, 
                   task = task, 
                   resampling = rdesc,
                   par.set = glm_ps, 
                   control = ctrl,
                   measures = list(auc, setAggregation(auc, test.sd)),
                   show.info = TRUE)

  return(glm_tuning$x)
}