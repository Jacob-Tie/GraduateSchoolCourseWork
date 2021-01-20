# File Descriptions
## PlinkToZarr.sh
This file turns a set of .bed .bim and .fam files into a zarr data file which python can read much more quickly.
## MagmaCommands.sh
This file outputs the genes each SNP in the data is directly associated with (is within its direct transcription region). This by no means captures all of the SNPs, and only covers around half. Better techniques for gene mapping may be available.
## TabNetTest.py/sh
This is where I showcase a highly custom neural network model which, while beneficial for some things like vastly reducing the time to run of these neural networks, was unable to get results better than a baseline solution. The basic structure is as follows:

We start with an input of SNP values. Each SNP value is assumed to be exactly within a single gene, so a mapping can be made from SNPs to genes (multiple SNPs can point to the same gene). This allows us to conceptualize a sparse matrix of weights for a neural net, where each input is thought of as a SNP and each output is thought of as a gene such that the only possible connections from SNP to gene are for SNPs that can be exactly mapped to that gene. This is a form of biologically minded dimensionality reduction, and depending on the complexity of the gene different "interpretive" neural networks can be applied (for instance some genes may have several hundred SNPs connected to them, so we can create a larger neural network in between the input and this gene if it is helpful, see small, medium and large blocks (all classes) in the code), but we found that direct connections vs more complex networks had no difference in predictive power. Any neural network can be placed after this, and we attempted to use Google's TabNet ([1]) (using the implementation found here: https://github.com/dreamquark-ai/tabnet). We believe that hyperparameter optimization might be crucial for this method, and given time and computational restrictions we were unable to finish this project, but the idea seems promising enough to work with future projects. NOTE: the .sh file is merely for job submission.
## Sources
[1] arXiv:1908.07442 [cs.LG]
