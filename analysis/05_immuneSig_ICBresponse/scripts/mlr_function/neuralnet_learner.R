## logistic regression hypermeter tuning

neuralnet_tuning <- function(traindata, len){
  
  
  set.seed(123)
  
  mod_neuralnet = makeLearner("classif.neuralnet", 
                        predict.type = "prob", 
                        fix.factors.prediction = TRUE)
  
  mod_neuralnet = makeUndersampleWrapper(mod_neuralnet,  usw.rate = 1/(N_num/P_num))

  neuralnet_ps <- ParamHelpers::makeParamSet(
              ParamHelpers::makeIntegerVectorParam('hidden',
                                                  len = len,
                                                  lower = 1L,
                                                  upper = 10L),
             ParamHelpers::makeIntegerParam("stepmax", 
                                                  lower = 1L,
                                                  upper = 5L,
                                                  trafo = function(x) {x}*10^5),
             ParamHelpers::makeDiscreteParam("algorithm",
                                  values = c( 'backprop', 'rprop+', 'rprop-', 'sag', 'slr')),
             ParamHelpers::makeDiscreteParam("act.fct",
                                  values = c( 'logistic', 'tanh'))

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

  # ctrl = generateGridDesign(neuralnet_ps)

  neuralnet_tuning = tuneParams(mod_neuralnet, 
                   task = task, 
                   resampling = rdesc,
                   par.set = neuralnet_ps, 
                   control =  ctrl, 
                   measures = list(auc, setAggregation(auc, test.sd)),
                   show.info = TRUE)

  return(neuralnet_tuning)
}
