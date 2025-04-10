Run deconv.py to deconvolute gene expression for Glioblastoma (GBM).

In linux/mac deconv.py can be run from terminal as:
python deconv.py

Deconvolution requires three input files: RNA-seq data, tumor purity and copy number data.

The desired output is in output folder.

In order to run deconvolution with bootstrap issue the following command:
python deconv_bootstrap.py

The file deconv_bootstrap.csv.bz2 is the output from all bootstrap analysis.

Dependencies: 
The script depends on pandas and scipy libraries which can be installed in linux/mac by running:
sudo pip install pandas scipy

In ubuntu/debian, pip can be installed by running the following command in the terminal:
sudo apt install python-pip
In case of Red Hat/CentOS, the command is:
sudo yum install python-pip


