## logistic regression hypermeter tuning

nnTrain_tuning <- function(traindata, len){
  
  
  set.seed(123)
  
  mod_nnTrain = makeLearner("classif.nnTrain", 
                        predict.type = "prob", 
                        fix.factors.prediction = TRUE)
  
  mod_nnTrain = makeUndersampleWrapper(mod_nnTrain,  usw.rate = 1/(N_num/P_num))

  nnTrain_ps <- ParamHelpers::makeParamSet(

             ParamHelpers::makeIntegerParam('hidden',
                                             lower = 1L,
                                             upper = 20L,
                                             trafo = function(x)c({x})),
             ParamHelpers::makeNumericParam("learningrate", 
                                              lower = 0,
                                              upper = 1,
                                              trafo = function(x) {x}*0.1),
             ParamHelpers::makeNumericParam("momentum", 
                                              lower = 0,
                                              upper = 1,
                                              trafo = function(x) {x}*0.1),
             ParamHelpers::makeIntegerParam("numepochs", 
                                              lower = 1L,
                                              upper = 10L),
             ParamHelpers::makeIntegerParam("batchsize", 
                                              lower = 5L,
                                              upper = 20L,
                                              trafo = function(x) {x}*10),

             ParamHelpers::makeDiscreteParam("activationfun",
                                  values = c( "sigm", "linear", "tanh"))
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

  # ctrl = generateGridDesign(nnTrain_ps)

  nnTrain_tuning = tuneParams(mod_nnTrain, 
                   task = task, 
                   resampling = rdesc,
                   par.set = nnTrain_ps, 
                   control =  ctrl, 
                   measures = list(auc, setAggregation(auc, test.sd)),
                   show.info = TRUE)

  return(nnTrain_tuning$x)
}