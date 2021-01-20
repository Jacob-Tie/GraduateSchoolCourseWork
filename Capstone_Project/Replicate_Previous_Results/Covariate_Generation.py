#Note: Not all of these import statements are required, this is just my default list of imports
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import time
from sklearn.model_selection import train_test_split
import zarr
import csv
import progressbar

#Load the Zarr files created in QC
def loadDataDict():
    output = {}
    for i in range (1, 23):
        output[str(i)] = zarr.open('PATH_TO_TRAINING_DATA'+str(i), mode='r')
    return output
data = loadDataDict()

#read in an aligned list of covariates, gives dictionaries with keys as the iids, and values as the gender or age values
with open("COVARIATES_FILE.txt") as fp:
    reader=csv.reader(fp,delimiter='\t')
    age = {}
    sex = {}
    for row in reader:
        age[row[1]] = row[3]
        sex[row[1]] = row[15]
fp.close()


#Need to reorder these so they match the order of our Zarr array:
ageOrdered = {}
sexOrdered = {}

indvWithoutData = []

iids = np.array(data['2'].iid)

#Some individuals will not have data (which would throw an exception in np.float), we will need to keep track of these individuals
for i in iids:
    try:
        ageOrdered[i] = np.float(age[i])
        sexOrdered[i] = np.float(sex[i])
    except:
        indvWithoutData.append(i)
#Find the index numbers for the individuals without covariate data
#Use binary search to speed up this process:
#----------------
#from https://www.geeksforgeeks.org/binary-search-bisect-in-python/
#----------------
from bisect import bisect_left 
"""
This function preforms binary search
@param a: this array we are preforming binary search on
@param x: The value we are trying to find

Returns: -1 if item not found and the index of the item if it is found
"""
def BinarySearch(a, x): 
    i = bisect_left(a, x) 
    if i != len(a) and a[i] == x: 
        return i 
    else: 
        return -1
indvWithoutData = sorted(indvWithoutData)

relIndiv = np.ones(iids.shape[0], dtype = bool)
for i in range(iids.shape[0]):
    if BinarySearch(indvWithoutData, iids[i]) != -1:
        relIndiv[i] = False

#We will also format the target:
#Load in the target
target = {}
with open("PATH_TO_TARGET.fam") as fp:
    reader=csv.reader(fp,delimiter=' ')
    for row in reader:
        target[row[0]] = np.float(row[5])
fp.close()

#-9 corresponds to a missing height value, remove these from the relevant rows
x = np.array(list(target.values()))
for i in range(x.shape[0]):
    if i == -9:
        relIndiv[i] = False

#Save all of the covariate's and relevant individual's .npy files for later use:
np.save('PATH_TO_AGE_COVARIATE_TRAIN.npy', ageOrdered)
np.save('PATH_TO_GENDER_COVARIATE_TRAIN.npy', sexOrdered)
np.save('PATH_TO_RELEVANT_ROWS_TRAIN.npy', relIndiv)

#-----------------------------------------------
#Now do the test set covs:
#-----------------------------------------------
def loadDataDict():
    output = {}
    for i in range (1, 23):
        output[str(i)] = zarr.open('PATH_TO_TEST_ZARR'+str(i), mode='r')
    return output
data = loadDataDict()

#read in an aligned list of covariates
with open("PATH_TO_COVARIATES.txt") as fp:
    reader=csv.reader(fp,delimiter='\t')
    #keeps track of progress bar progress
    i = 0
    age = {}
    sex = {}
    for row in reader:
        age[row[1]] = row[3]
        sex[row[1]] = row[15]
        i += 1
fp.close()

#Need to reorder these so they match the order of our Zarr array:
ageOrdered = {}
sexOrdered = {}

indvWithoutData = []

iids = np.array(data['2'].iid)

for i in iids:
    try:
        ageOrdered[i] = np.float(age[i])
        sexOrdered[i] = np.float(sex[i])
    except:
        indvWithoutData.append(i)

#from https://www.geeksforgeeks.org/binary-search-bisect-in-python/
#----------------
from bisect import bisect_left 
"""
This function preforms binary search
@param a: this array we are preforming binary search on
@param x: The value we are trying to find

Returns: -1 if item not found and the index of the item if it is found
"""
def BinarySearch(a, x): 
    i = bisect_left(a, x) 
    if i != len(a) and a[i] == x: 
        return i 
    else: 
        return -1
indvWithoutData = sorted(indvWithoutData)

relIndiv = np.ones(iids.shape[0], dtype = bool)
for i in range(iids.shape[0]):
    if BinarySearch(indvWithoutData, iids[i]) != -1:
        relIndiv[i] = False

#Load in the target
target = {}
with open("PATH_TO_TEST_TARGET.fam") as fp:
    reader=csv.reader(fp,delimiter=' ')
    for row in reader:
        target[row[0]] = np.float(row[5])
fp.close()

#-9 corresponds to a missing height value, remove these from the relevant rows
x = np.array(list(target.values()))
for i in range(x.shape[0]):
    if i == -9:
        relIndiv[i] = False

np.save('PATH_TO_AGE_COVARIATE_TEST.npy', ageOrdered)
np.save('PATH_TO_GENDER_COVARIATE_TEST.npy', sexOrdered)
np.save('PATH_TO_RELEVANT_ROWS_TEST.npy', relIndiv)