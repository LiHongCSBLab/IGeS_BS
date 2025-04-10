nohup Rscript --vanilla metaPlusFilter/ITSssgsea_CombatRef_17datasets_finalmodel_metaPlusFilter.R > finalmodel_metaPlusFilter.log 2>&1&
nohup Rscript --vanilla metaPlusFilter/ITSssgsea_CombatRef_17datasets_finalmodel_ridge_metaPlusFilter.R > finalridge_metaPlusFilter.log 2>&1&
nohup Rscript --vanilla ITSssgsea_CombatRef_17datasets_finalmodel.R > finalmodel.log 2>&1&
nohup Rscript --vanilla ITSssgsea_CombatRef_17datasets_ridge_glmnet.R > finalridge.log 2>&1&

nohup Rscript --vanilla metaPlusFilter/ITSssgsea_CombatRef_17datasets_CV_metaPlusFilter.R > CV_metaPlusFilter.log 2>&1&
nohup Rscript --vanilla ITSssgsea_CombatRef_17datasets_CV.R > CV.log 2>&1&
