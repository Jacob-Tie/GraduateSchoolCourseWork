# Getting Started in Computational Genetics in Python
## Introduction:
I am in the process of wrapping up coursework for a degree in applied mathematics which involved a cumulative project in which I explored genomic data. I believe that I have learned enough of the tricks of working with genomic data to provide a post about ‘what I wish I knew’. In this spirit I will aim this article at people who have studied computer science from a generalist perspective but are now trying to dive into genomic data. To start, I will give a brief introduction to working with cluster-based computing since most genomic data is too large in scale to work with on a PC. I will specifically focus slurm-based computer systems (https://en.wikipedia.org/wiki/Slurm_Workload_Manager) since it is used in 60% of the top 500 supercomputers, so there is a good chance you’ll encounter it.
## Overview of this article
We start with some workflow ideas, such as working on a remote computer (using jupyterhub, and scheduling jobs via Slurm). Then we discuss details specific to genomic data, and in particular the care we must take since the datasets are so large.  
  
TL;DR: use slurm, conda (for package management), plink, zarr (not xarray), and pytorch’s iterable datasets data too large for memory.
## Job Submission and Scheduling with Slurm
As a ‘step zero’, get comfortable with the documentation on your computing cluster! These clusters are usually run by professionals and have good documentation. For example, here at the University of Colorado Boulder, I became very friendly with the website https://curc.readthedocs.io/en/latest/. What I say below may not apply to your computing resource, so verify with the documentation.  
  
Another important part of reading this documentation would be to learn what non-command line tools are available to you. Much of the work I had done in data science before graduate school was done in Jupyter notebooks which made prototyping much easier than using Python scripts. Some institutions have portals which allow you to access Jupyter notebooks, MS Visual Studio, or other interactive coding environments (though these might have limited computational resources), which may help with small scale tests or just getting a feel for the data you are working with.  
  
Aside from this, the first thing that I wish I knew at the start of this degree would have been how to use slurm, a job scheduling program used in situations where multiple users may be vying for computational resources. At its core, slurm is a job manager which takes many requests all with different computational needs and priority levels within the system and outputs the most efficient queue to execute those computations. There are many tricks and rules of best practice for how many resources to ask for (as in processors vs nodes vs RAM per node), but when starting this project, the most helpful thing was a brief explanation of how job submission worked, and an example script for how to submit a job (this documentation helped me a lot https://curc.readthedocs.io/en/latest/running-jobs/batch-jobs.html).  
  
Job submission itself is easy once you understand that you will not be submitting a python, R, Matlab, etc. script directly to the job scheduler. Instead, you will submit a bash file (https://en.wikipedia.org/wiki/Bash_(Unix_shell)) which will carry a set of console commands which will then execute a script that you want to execute (ie you could call python3 python_script.py within the bash script to have it execute a python script). Built into this bash file are also flags which will tell the computer how many computational resources your job needs. The most efficient way to understand this is to just look at a full example of a typical job script:  

    #!/bin/bash
    #SBATCH --qos=preemptable
    #SBATCH --time=04:00:00
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --job-name=plinkToZarrQC
    #SBATCH --output=plinkToZarrQC.%j.out



    source PATH_TO_ANACONDA
    conda PATH_TO_YOUR_ENVIRONMENT

    python3 plinkToZarrQC.py

  
This is an example of a job script which will open an anaconda environment then run a script called plinkToZarrQC.py. The header is the most important part from slurm’s perspective, so I will quickly explain each flag:  
  
#!/bin/bash – This line is required as the heading of every job submission file. This tells linux what shell to use (for most shells you will use \bin\bash) (see https://en.wikipedia.org/wiki/Shebang_(Unix) for details).  
  
#SBATCH --qos=preemptable – This sets the “Quality of Service” for your job and varies by institution. Think of a QoS as a queue or a different line at the grocery store: each cluster usually has several different queues. I have set this to preemptable because at my institution there is a low priority queue that can be used for less important jobs, but you should check what qos to apply at your institution. Usually, your account will only have access to a few queues; for example, I have access to the generic “preemptable” queue and my department’s queue, but I will not have access to, say, a high-priority queue that a research group may have purchased.  
  
#SBATCH --time=04$:$00$:$00 – This sets the maximum time that the job can run. If the job takes longer than this time (4 hours in this case) then it will be automatically canceled. So, this should be an overestimate! Keep in mind, though, that if you request a very long walltime, your job may have a lower priority and take longer in the queue.  
  
#SBATCH --nodes=1 – This tells Slurm how many computational nodes we are requesting. This can be thought of as the number of individual computers that you would like to work on your task. Were you to request more nodes then you would have access to the computational resources of more than one computer. For instance, if all the nodes on your computer cluster only had access to 1 GB of RAM, but you needed 4GB of RAM you would need to ask to use 4 nodes. Keep in mind that this would also require specific programming to actually make use for 4 nodes.  
  
#SBATCH --ntasks=1 – This tells Slurm how many processors we are requesting. This should be 1 if your program does not take advantage of multithreading. Be aware that even if your code is not multithreaded, you may call libraries (numpy, pytorch, etc.) that are multithreaded.  
  
#SBATCH --job-name=plinkToZarrQC – This tells Slurm what to call your job while it’s running and while it is in the queue.
#SBATCH --output=plinkToZarrQC.%j.out – This tells Slurm where to write the console output. For instance, if we run a python script in this file that prints “hello” to the console this is where that output would be written and saved to. Note: this file is created as soon as your job is submitted but is sometimes not populated until after the job has finished running. Also note that the %j refers to the job number.  
  
As a brief aside, the clusters I used in graduate school operated under a ‘modules’ system (https://modules.readthedocs.io/en/latest/module.html), which allows for modular configurations of what is and isn’t loaded into the environment. One byproduct of this is that using anaconda environments (https://docs.anaconda.com/anaconda/user-guide/getting-started/) is very easy. If you have the ability to set up one of these environments on your system, I would definitely consider it since this step makes setting up the libraries that I will mention later much easier than having to rely on your system’s administrators to set these up for you. If anaconda isn’t an option for you here is an article with some alternatives: https://medium.com/@krishnaregmi/pipenv-vs-virtualenv-vs-conda-environment-3dde3f6869ed.  
  
The other major trick that you can do with slurm is to submit so called ‘batch jobs’. These are multiple jobs all being submitted at the same time (but independently from each other), so if your problem can be parallelized effectively you can drastically cut down on execution time. In genetics and genomics applications, chromosomes are an example of a natural way to parallelize your scripts, because each of the 23 chromosomes (22 autosomes and one sex chromosome) can be run independently. Here is an example of one such job:  

    #!/bin/bash
    #SBATCH --qos=preemptable
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=8
    #SBATCH --time=9:00:00
    #SBATCH --job-name=ReplicateQC
    #SBATCH --output=./slurm_out/ReplicateQC.%A-%a.out
    #SBATCH --array=1-22

    gunzip -c PATH_TO_ZIPPED_UKBIOBANK_DATA_SPLIT_BY_CHROM_${SLURM_ARRAY_TASK_ID}.bed.gz > OUTPUT_BED.bed
    gunzip -c PATH_TO_ZIPPED_UKBIOBANK_DATA_SPLIT_BY_CHROM_${SLURM_ARRAY_TASK_ID}_v2.bim.gz > OUTPUT_BIM.bim
  
The only new information that is here is the flag #SBATCH --array=1-22, one for each of the 22 autosomes of the human genome. This will create 22 jobs labeled 1 through 22 within the system. Which job you are in specifically can be accessed through the {SLURM_ARRAY_TASK_ID} variable, so what this is doing to taking data from the UKBioBank that is split into 22 files (by chromosome), and it is unzipping them all simultaneously. Note that this is not the best use of this procedure (this script does not take long to run as a single job with a for loop), but it does illustrate the power of batch job submission.  
  
In general, this should be enough to get started with slurm job submission (to submit the job script bash file you will need to run a command that loads slurm like module load slurm, and then **sbatch <job_file.sh>** will submit the job to slurm from the console), so now you should be ready to experiment with jobs on whatever system you need to.  
  
Another thing to consider is that these datasets can be hard to store on hard drives due to their size. Many shared computing structures will give you a personal partition, but even moderately sized genomic datasets may take all or most of this space. If this is a problem than one potential solution is scratch drives which are just shared hard drives that system admins will occasionally wipe, but these typically have much more space than personal drives and often have faster i/o properties and are not backed up; any kind of temporary data should be written to the scratch drive.
## Into the Specifics: Preprocessing Data with Plink
Computational geneticists have developed many tools that may save you substantial time writing your own scripts if you know where to look. Possibly the most widely used of these tools is plink (https://www.cog-genomics.org/plink/2.0/). Plink is a means of working a particular set of file types which you may have already encountered when looking for genomic datasets. Most often, these are a trio of files in the form of .bed, .bim, and .fam files, though additional file formats (https://www.cog-genomics.org/plink/2.0/input) can be utilized by plink2.  
  
The .bed file contains information in binary (2bit binary to be precise) which represents what specific nucleotide each SNP that your dataset keeps track of is. For instance, if your dataset has SNP rs1815739 then the information of whether this is an A, T, C, or G for all of the individuals in the dataset will be found in this file. This means that the .bed file is typically the largest of the three.  
  
The .bim file contains metadata about the SNPs including chromosome codes, variant identifier, position, base-pair coordinate, and which allele it is on.  
  
The .fam file contains metadata about the individuals in the dataset including family ID, within-family ID, within-family ID of the father/mother, sex code, and phenotype values.  
  
Plink (found here: https://www.cog-genomics.org/plink/) is not only a standard way of storing genomic data/metadata, but it also provides a means of doing basic manipulations and extensive analyses of the data as well as genome wide association studies (GWAS) which find which SNPs in your data are important to the phenotype of interest. For instance, in many papers you may see preprocessing steps like “remove SNPs with minor allele frequencies <.1%” or “Remove individuals with missing call rates >10%”. Both of these steps can be done easily in plink. To see how let us look at an example script that does some preprocessing steps that are common in plink:  

    plink1.9 --bed OUTPUT_BED.bed\
                              --bim OUTPUT_BIM.bim\
                              --fam PATH_TO_FAM_FILE.fam\
                              --pheno PATH_TO_ANY_OTHER_TARGET_VARS_TO_INCLUDE.txt\
                              --geno .03\
                              --mind .1\
                              --maf 0.001\
                              --keep PATH_TO_INDIV_IN_TRAIN_SET.txt\
                              --linear\
                              --write-snplist\
                              --out DESIRED_PATH_TO_OUTPUT_${SLURM_ARRAY_TASK_ID}

  

Plink1.9 is a stand in for wherever you store plink on your computer (this should be a file path). The bed flag simply tells plink where the unzipped .bed file is. The bim flag is like the bed flag. The fam flag tells plink where the .fam file is. Finally, the pheno flag tells plink to look at an additional file of phenotypes. For my setup there was a separate file with the height phenotype in it. This file is a txt file that with at least three columns, which should be structured as follows:  
  
iid fid new phenotype (Ex: 111111 111111 1.6)  
Note: The iid (individual ID) and the fid (family ID) are often treated as the same, when no family structure is present or known within your dataset.  
  
This flag also tells plink which phenotype you would like it to do any GWAS studies on. Now that plink knows what data we are working with we can get to the manipulations:  
  
--geno: this flag tells plink to only consider SNPs with missing call rates less than .03 (3 percent).  
--mind: this flag tells plink to only keep individuals with less than 10 percent missing information in their SNP values.  
--maf: this flag tells plink to remove SNPs with minor allele frequencies at .1 percent   
--keep: this flag makes sure that, for all these calculations, that plink is only using the training set. This takes as its argument a text file of the form: iid fid.  
--linear: tells plink to perform a single marker GWAS on the data that is left. This is important since many papers will prompt you to “only select the top N SNPs from the dataset by p-value”.  
--write-snplist: tells plink to write a file that contains a list of all the SNPs that made it through these steps.  
--out: tells plink where to write its output (this will be multiple files, the most important of which will be one ending in .assoc.linear, this is the results of the GWAS which will contain p-values for whether or not a particular SNP had an effect on the phenotype you are interested in studying).  
  
There is also a flag (not in this script) for only including specific SNPs in the data which can be used to filter out only the top N SNPs by p-value. In general, you will first run a command like this to get p-values for each SNP, then you will create a txt file with the SNP names of the top N SNPs by p-value, then run another plink call that will throw out everything except these SNPs. This encapsulates most of the basics of plink and brings us to the end of what can be done at the command line. The rest of this article will be focused on working in python with genomic data. Since I believe python is a relatively common language my code examples will become sparser in the subsequent section, but I will provide a link to a github with a full project worked out at the end.
## Working with Genomic Data in Python
Tldr: Never make your own datatype or data reader functions! Use libraries like pandas-plink to read the data in, convert to zarr arrays to make reading data faster, and use PyTorch’s iterable datasets to make piecemeal reading data into RAM faster for datasets that cannot fit into memory.  
  
The first thing that I wanted to mention was how to properly separate test and train sets for your machine learning models. If you have not worked with a lot of large datasets before than you might make many copies of the data, or even just use seeded random numbers to make the test train split, but this is not ideal since random seeds can change across libraries and making copies of the data can now be hard. For this reason, it is best practice to keep master files with just the iids of the training and set and just the iids of the testing. This ensures that you never get strange results or accidently analyze different datasets for different libraries.  
  
Since we are working with plink data types, we will need a special way of reading in plink data into python. Fortunately, there is such a library for translating plink data into something usable by python: pandas-plink (https://pypi.org/project/pandas-plink/). There is, however, a problem with using this. Pandas plink uses xarrays to store the panda’s data (assuming it is too large to load into memory), which skirts the problem of large data, but presents a new difficulty: loading data into memory. The problem is that since xarrays are basically just python datatypes that point to an area in the hard drive without loading into RAM they will have to sift through memory every time that you want to load any portion of the data from the hard drive into RAM. This is especially slow since the way xarrays are structured do not lend themselves to almost any procedure that might speed this process up. Therefore, I advocate for the use of Zarr arrays (https://zarr.readthedocs.io/en/stable/tutorial.html). (which was suggested to me by Richard Border https://scholar.google.com/citations?user=fQhvPM8AAAAJ&hl=en) These are arrays that use and HDF5 (https://en.wikipedia.org/wiki/Hierarchical_Data_Format) like hierarchy for storing their data, and practically speaking are much faster than dask or xarray. Here is a useful script for changing plink data into zarr arrays (which is a modified version of code that Richard Border showed me when teaching me how to use zarr arrays):  
  
    import pandas_plink as pp
    from dask.diagnostics import ProgressBar
    for i in range(1, 23):
            G = pp.read_plink1_bin("PATH_TO_BED_TRAIN" + str(i) + '.bed')
            G = G.astype('int8')
            G = G.to_dataset()
            with ProgressBar():
                G.compute()
            with ProgressBar():
                G.to_zarr('PATH_TO_TRAINING_DATA'+str(i))
  
Note: we store the integers as int8 since this is the lowest precision that python can use; for large datasets this is a crucial step. This datatype can be used in a stepwise fashion so you only have to load in a chunk of the data at a time, and they can even be easily applied to certain python iterables like PyTorch’s iterable dataset class (here is a discussion on the PyTorch forums with an example: https://discuss.pytorch.org/t/example-for-torch-utils-data-iterabledataset/101175/3), but I will discuss options for large data later. This is about all there is to loading data into python, now I would like to go over some specific preprocessing techniques that you are liable to see.
## Preprocessing with Python
The last part of preprocessing that might be important to your applications is the processing of the target variable. Since many phenotypes have significant covariates which can confuse algorithms which only look at genomic data it may be important to abstract away many of these covariates from your target phenotype. There are many ways in which to do this, but they all involve some sort of standardization based around the covariate. Take height for example. This has two very significant covariates in the form of age and sex, each of which should be accounted for. [1] is recent paper in the business of predicting height using genomic data which offers a scheme for accounting for these two covariates.  
  
[1] accounts for sex by normalizing men and women in the data separately. This works by splitting the data into self-identified men and women and then normalizing each of these splits as if it were its own dataset (by taking each entry and subtracting the mean and dividing by the standard deviation). This may look like a strange step since we are taking two different sections of the same distribution and assuming that they are the same after the normalization, but since our goal is to isolate the effects of the genome rather than these covariates this change should be somewhat intuitive.  
  
Age is accounted for in [1] by subtracting a linear trend trained on the phenotype with age as the only predictor (note that this is done after normalizing based off of sex). This means that they first train a linear model to predict height based on age, then feed each individual’s age into the model and subtract this model’s output from the phenotype.  
  
There is no standard way to account for covariates, but this is a crucial step for many models so it is something you should consider in your projects. These two steps are often useful for many phenotypes but depending on the nature of your problem you may have to apply different methods to account for many covariates, but this is, again, a good starting point.
This brings me to the final topic of this discussion which is just a list of some ways to deal with large data in python, and some problems you may encounter.
## Big Data in Python, Approaches:
Depending on the nature of your computational resources and dataset size this may or may not be a problem that occupies much of your time. If your whole dataset fits easily into memory then this is a non-issue; however if not then there are a few ways you can deal with this problem. You may also find that your dataset fits into RAM and works fine when you run on a CPU, but when you move to a GPU, the GPU’s lower memory suddenly requires these big-memory considerations.  
  
The first scenario that I will discuss is if your data fits into memory in a low precision form but cannot fit if stored in higher precision. I encountered this with the UKBioBank data since the computers that I was working with could easily fit ~80GB of data into RAM (which was about the size of the section of the UKBioBank data I was working with) but could not fit it when I wanted to normalize it. This is because normalization in python (for most implementations including scikit learn’s implementation) will automatically upscale the precision.  
  
To mitigate this problem, we can always work with batched methods where you keep a ‘master’ dataset at low precision but convert it to high precision in batches when feeding it into the model. This works especially well with deep learning since any method that uses neural networks can be batched. There are even ways of doing this with linear regression, though you would have to implement them by hand for things like lasso regression (see proximal gradient descent for lasso regression if you are interested https://en.wikipedia.org/wiki/Proximal_gradient_methods_for_learning). This, however, is a band aid to a much more likely problem: what do I do if I cannot fit this dataset into memory?  
  
This problem is much more difficult, but I will advocate for using PyTorch’s iterable datasets, even if your problem does not require deep learning or anything that PyTorch is typically used for. This is because PyTorch has already implemented a multiprocess loading scheme for this type of process. To illustrate why this is important let me describe the way that chunking works in Zarr for this type of problem.  
  
If a dataset is too large to be stored into memory, then Zarr will ‘chunk’ it, i.e., break it up into manageable sections that Zarr can easily find when it needs to. This is a good first step (now we only need to load in a single chunk at a time to do most processing jobs), but if you are running Zarr arrays with only a single process this may run into a problem. Consider trying to take the mean of your large data matrix. This will require you to visit each data point in turn, so Zarr will load chunk after chunk and slowly compute the mean. If a single process is doing this then it might take a very long time since a single process will load the data in, delete it from RAM when it is no longer needed, then find and load in the next set of data. If instead we were running multiple processes then we would have a process that would sit in memory right where it needed to be, then, it can load its data in parallel to another process. This is a parallelized implementation of data streaming, so you should get all the speed benefits from parallelized processes simply by using this class.
One thing to note is that if you do use this approach you will be working with PyTorch tensors, but since these can be easily converted to numpy it should work for general data loading problems. These are the main ways that you can deal with large data in Python, but you could also make your own parallelized implementation.  
  
I left this as a small section since I had almost no need for this in my project, but it is possible to do parallelization in Python, but it is not easy. Python uses a global interpreter lock which prevents many tasks from being performed simultaneously, though there are libraries that do support some tasks in multiprocessing. In general, I would try to avoid coding your own parallelized methods in Python if possible since it is more trouble than in other languages, but the functionality is still there if you need it.
## A Fully Worked Example
Here is a link to a github with a couple of examples that I worked on throughout this project: https://github.com/Jacob-Tie/GraduateSchoolCourseWork/tree/master/Capstone_Project. 
## Sources Cited
[1] Lello, L., Avery, S. G., Tellier, L., Vazquez, A. I., Campos, G. D., & Hsu, S. D. (2017). Accurate Genomic Prediction Of Human Height. doi:10.1101/190124