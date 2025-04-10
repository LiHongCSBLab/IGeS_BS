## logistic regression hypermeter tuning

nnet_tuning <- function(traindata, len){
  
  
  set.seed(123)
  
  mod_nnet = makeLearner("classif.nnet", 
                        predict.type = "prob", 
                        fix.factors.prediction = TRUE)
  
  mod_nnet = makeUndersampleWrapper(mod_nnet,  usw.rate = 1/(N_num/P_num))

  nnet_ps <- ParamHelpers::makeParamSet(
              ParamHelpers::makeIntegerParam('size',
                                              lower = 1L,
                                              upper = 5L),
             ParamHelpers::makeIntegerParam("rang", 
                                              lower = 1L,
                                              upper = 5L,
                                              trafo = function(x) {x}*0.1),
             ParamHelpers::makeNumericParam("decay", 
                                              lower = 0,
                                              upper = 1),
             ParamHelpers::makeIntegerParam("maxit", 
                                              lower = 1L,
                                              upper = 5L,
                                              trafo = function(x) {x}*100)
             #ParamHelpers::makeDiscreteParam("softmax",
             #                     values = c( 'linout', 'entropy', 'softmax', 'censored'))
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

  # ctrl = generateGridDesign(nnet_ps)

  nnet_tuning = tuneParams(mod_nnet, 
                   task = task, 
                   resampling = rdesc,
                   par.set = nnet_ps, 
                   control =  ctrl, 
                   measures = list(auc, setAggregation(auc, test.sd)),
                   show.info = TRUE)

  return(nnet_tuning$x)
}