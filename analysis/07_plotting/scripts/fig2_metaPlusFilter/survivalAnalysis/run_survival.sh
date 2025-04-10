#!/bin/bash

source /home/fengfangyoumin/miniconda3/bin/activate r_env

Rscript --vanilla Braun_dataset.R > braun.log 2>&1&
Rscript --vanilla Gide_dataset.R > gide.log 2>&1&
Rscript --vanilla GSE126044.R > GSE126044.log 2>&1&
Rscript --vanilla GSE135222.R > GSE135222.log 2>&1&
Rscript --vanilla Hugo_dataset.R > hugo.log 2>&1&
Rscript --vanilla IMvigor210_dataset.R > IMvigor210.log 2>&1&
Rscript --vanilla Nathanson_dataset.R > nathanson.log 2>&1&
Rscript --vanilla PMID_29301960.R > PMID_29301960.log 2>&1&
Rscript --vanilla Riaz_dataset.R > riaz.log 2>&1&
Rscript --vanilla Zhao_dataset.R > zhao.log 2>&1&

Rscript --vanilla GSE176307.R > GSE176307.log 2>&1&
Rscript --vanilla Lauss.R > lauss.log 2>&1&
Rscript --vanilla Liu.R > liu.log 2>&1&
