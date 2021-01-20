# Quality Control on UKBIOBANK Data
## File Descriptions
### Create_Test_Train_Split.py
This file creates two output files, one with the ids of the individuals in the training set and one with the individuals in the test set (note that the test set will be split into test and validation later). These files are of the format the plink expects (so the ids are repeated on each line). Note: there is no job submission file because I ran this in an interactive session.
### Replicate_QC.sh
Replicates basic quality control steps detailed by [1]. I account for covariates in the scripts which run the machine learning algorithms, so this is just the initial selection of individuals, and making sure all of the SNPs have high call rates. The most important output from this file is that of the GWAS output. This will give a .assoc.linear file which will give the single SNP p-values for a given trait (in our case height).
### ChooseTopNSNPs.py
In this file I output a .txt file with the top 50,000 SNPs associated with height (using the p-value to determine which SNPs are most related, this is indicated by [1]). Note: there is no job submission file because I ran this in an interactive session.
### Output_Final_Bed_Train.sh
This will output the final .bed files using only the SNPs found in ChooseTopNSNPs.py, and only the individuals in the train set who passed the quality control steps detailed before.
### Output_Final_Bed_Test.sh and Remove_Files.sh
This will output the final .bed files using only the SNPs found in ChooseTopNSNPs.py, and only the individuals in the test set. Note: in most of the other files of this type (which use the base UKBIOBANK data) I remove the created .bed and .bim files that we gunzip at the start. This runs into a problem for this script since it runs so quick that it will sometimes remove files before some of the jobs start, so I have a separate file Remove_Files.sh to remove those gunzipped files after use.
### plinkToZarrQC.py and .sh
This will create a zarr array of the data we just created in the Output_Final_Bed files. This will use 8 bit integers, and will greatly increase the speed that python reads in the data. The .sh file is for job submission. Note: as written this will only create zarr files for the train set, but a simple change of file paths will give the test set as well (which I would recommend to keep file clutter low).
##### Sources
[1] Lello, L., Avery, S. G., Tellier, L., Vazquez, A. I., Campos, G. D., & Hsu, S. D. (2017). Accurate Genomic Prediction Of Human Height. doi:10.1101/190124
