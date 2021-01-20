import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import time
from sklearn.model_selection import train_test_split
import zarr
import csv
from pytorch_tabnet.tab_network import TabNetNoEmbeddings
from torch.nn.utils import clip_grad_norm_

"""
loadDataDict -- Loads each chromosome's data into a dictionary. The keys are strings of the chromosome's number (starting at '1'), and the values are zarr arrays
"""
def loadDataDict(name):
    output = {}
    for i in range (1, 23):
        output[str(i)] = zarr.open(name+str(i), mode='r')
    return output
data = loadDataDict('PATH_TO_TRAINING_DATA') #Contains training data
test_data = loadDataDict("PATH_TO_TESTING_DATA") #Contains data on the hold out sets, will be split into test and validation sets later in this script

#output.txt keeps track of progress as the program is running
f = open("./output.txt", "w")
f.write("1\n")
f.close()


#Load train Covariates, see the python script on covariate generation
age = np.load('PATH_TO_AGE_COVARIATE_TRAIN.npy',allow_pickle='TRUE').item() 
age = np.array(list(age.values()))
sex = np.load('PATH_TO_GENDER_COVARIATE_TRAIN.npy',allow_pickle='TRUE').item()
sex = np.array(list(sex.values()))
relRows = np.load('PATH_TO_RELEVANT_ROWS_TRAIN.npy',allow_pickle='TRUE')

#Load test Covariates
age_test = np.load('PATH_TO_AGE_COVARIATE_TEST.npy',allow_pickle='TRUE').item()
age_test = np.array(list(age_test.values()))
sex_test = np.load('PATH_TO_GENDER_COVARIATE_TEST.npy',allow_pickle='TRUE').item()
sex_test = np.array(list(sex_test.values()))
relRows_test = np.load('PATH_TO_RELEVANT_ROWS_TEST.npy',allow_pickle='TRUE')


#Read the Target Variable:

#Load in the target as a dictionary, perhaps it would be better to read it in directly as a numpy array, but this helps with debugging
target = {}
with open("PATH_TO_FAM_FILE_CREATED_IN_QC") as fp:
    reader=csv.reader(fp,delimiter=' ')
    for row in reader:
        target[row[0]] = np.float(row[5])
fp.close()

target = np.array(list(target.values()))
#subset so that we only have values for which there are valid target values, and complete cases for the covariates
target = target[relRows]
print("The mean of the train target is: " + str(np.mean(target)))
print("The standard deviation of the train target is: " + str(np.std(target)))

#Load in the test target, with a similar procedure that was used for the train set
test_target = {}
with open("PATH_TO_FAM_FILE_CREATED_IN_QC_TEST_SET") as fp:
    reader=csv.reader(fp,delimiter=' ')
    for row in reader:
        test_target[row[0]] = np.float(row[5])
fp.close()

test_target = np.array(list(test_target.values()))

test_target = test_target[relRows_test]
print("The mean of the test target is: " + str(np.mean(test_target)))
print("The standard deviation of the test target is: " + str(np.std(test_target)))


#load the Zarr into a set of numpy arrays (meta data is not retained in these arrays). These arrays are in RAM so be sure that enough mem is allocated:
relData = {}
for chrom in data.keys():
    relData[chrom] = np.array(data[chrom].genotype[:,:])

relData_test = {}
for chrom in test_data.keys():
    relData_test[chrom] = np.array(test_data[chrom].genotype[:,:])
"""
This function will take a dictionary whose keys point to data and it will return a single dataset that represents the concatination of all the data

@param relData: the dictionary containing the data

Returns: a numpy array that has all the data joined column wise starting with the first key in the dictionary

Pseudo-Code:
First we find the dimensions of the final output array
for each chromosome make sure that there are valid data points (invalid data points are at -10)
if valid data is contained in chromosome then we know that the number of rows in the output concatinated array must be the rows of this chromosome's genotype array (this should be constant across all chrom)
if valid data is contained in chromosome then we must add the number of columns (number of SNPs) in this chromosome to what is the total number of columns in the output

create an output matrix of the correct size
Read in all of the data from the valid chromosomes to the output
return output
"""
def concatData(relData):
    rows = 0
    cols = 0
    for chrom in relData:
        if np.all(relData[chrom] != -10):
            rows = relData[chrom].shape[0]
            cols += relData[chrom].shape[1]
    output = np.zeros((rows,cols), dtype = np.int8)
    colStart = 0
    for chrom in relData:
        if np.all(relData[chrom] != -10):
            output[:, int(colStart):colStart+int(relData[chrom].shape[1])] = relData[chrom]
            colStart += relData[chrom].shape[1]
    return output

relData = concatData(relData)
relData_test = concatData(relData_test)

relData = relData[relRows, :]
relData_test = relData_test[relRows_test, :]

#Keep track of progress
f = open("./output.txt", "a")
f.write("Shape of the test data and test target: " + str(relData_test.shape) + " " + str(test_target.shape))
f.write("2\n")
f.close()

#eliminate rows with ages that have very few individuals represented. Note that these age cut offs were found via EDA earlier in this project, see covariate generation:
relRows = np.ones(relData.shape[0], dtype = bool)
for i in range(age.shape[0]):
    if age[i] <= 40 or age[i] >= 70:
        relRows[i] = False

#Subset the target, data, and covariate objects with the rows that we know are relevant
target = target[relRows]
relData = relData[relRows,:]
sex = sex[relRows]
age = age[relRows]

#Confirm that all of the data is of the correct shape
print("The shape of the relevant data/target/covariate arrays is:")
print(relData.shape)
print(target.shape)
print(sex.shape)
print(age.shape)

#Used to test a later method, left in because one could easily change the target and since it doesn't add much to the computional time it was left in
saved_target = np.copy(target)

#Now perform adjustments for covariates suggested by the paper:
#create Z-scores based on sex, there is a male Z-score and a female Z-score as recommended by the paper we are replicating
#Keep track of the male and female mean/standard deviations since we will need them to normalize the test set as well
mean = []
std = []
k = 0
for i in np.unique(sex):
    temp = np.zeros(age.shape[0], dtype = bool)
    for j in range(sex.shape[0]):
        if sex[j] == i:
            temp[j] = True
    mean.append(np.mean(target[temp]))
    std.append(np.std(target[temp]))
    target[temp] = (target[temp]-mean[k])/std[k]
    k+=1

#Do the same for test:
k = 0
for i in np.unique(sex_test):
    temp = np.zeros(age_test.shape[0], dtype = bool)
    for j in range(sex_test.shape[0]):
        if sex_test[j] == i:
            temp[j] = True
    test_target[temp] = (test_target[temp]-mean[k])/std[k]
    k+=1


#fit linear regression based on age and standardize (subtract this trend from the targets)
from sklearn.linear_model import LinearRegression
#NOTE: NO NEED TO SEED THIS SINCE STOCASTIC GRAD DESC NOT USED
reg = LinearRegression().fit(age.reshape(-1,1), target)

for i in range(target.shape[0]):
    temp = reg.predict(age[i].reshape(-1,1))
    target[i] = target[i] - temp

#Do the same for the test set:
for i in range(test_target.shape[0]):
    temp = reg.predict(age_test[i].reshape(-1,1))
    test_target[i] = test_target[i] - temp



#we will need to normalize reldata for each batch
from sklearn.preprocessing import Normalizer
transformer = Normalizer().fit(relData)

#Now we can train the model:
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

from sklearn import linear_model
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
lambdas = [.00000000085]
f = open("./output.txt", "a")
f.write("3\n")
f.close()

#-------------
#A useful transformation: back to centimeters
#-------------
def backToCm(y, mean, std, reg, age, sex):
    out = np.copy(y)
    for i in range(y.shape[0]):
        temp = reg.predict(age[i].reshape(-1,1))
        out[i] = out[i] + temp
    k = 0
    for i in np.unique(sex):
        temp = np.zeros(age.shape[0], dtype = bool)
        for j in range(sex.shape[0]):
            if sex[j] == i:
                temp[j] = True
        out[temp] = out[temp]*std[k]+mean[k]
        k+=1
    return out

#test this to make sure it is working as intended
print("Test of the backToCm method: " + str(np.sum(np.isclose(backToCm(target, mean, std, reg, age, sex), saved_target))))

#--------------
#try a baseline solution:
#--------------
baseline_predictor = np.zeros(target.shape[0])
print("The MSE of the baseline prediction (all zeros) is: " + str(mean_squared_error(target, baseline_predictor)))
f = open("./output.txt", "a")
f.write("Baseline prediction for train data in standardized units: "+str(mean_squared_error(target, baseline_predictor))+"\n")
f.write("Baseline prediction for train data in cm: " + str(mean_squared_error(backToCm(target, mean, std, reg, age, sex), backToCm(baseline_predictor, mean, std, reg, age, sex)))+"\n")
f.close()
print("Baseline prediction for train data in standardized units: "+str(mean_squared_error(target, baseline_predictor)))
print("Baseline prediction for train data in cm: " + str(mean_squared_error(backToCm(target, mean, std, reg, age, sex), backToCm(baseline_predictor, mean, std, reg, age, sex))))

baseline_predictor = np.zeros(test_target.shape[0])
print("Baseline prediction for test data in standardized units: "+str(mean_squared_error(test_target, baseline_predictor)))
print("Baseline prediction for train data in cm: " + str(mean_squared_error(backToCm(test_target, mean, std, reg, age_test, sex_test), backToCm(baseline_predictor, mean, std, reg, age_test, sex_test))))



#--------------
#train the real model
#--------------
relData_test = transformer.transform(relData_test)

#Now we can train the model:
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

from sklearn.linear_model import SGDRegressor
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
batch_size = 64
num_steps = int(np.ceil(relData.shape[0]/batch_size))
epochs = 50
f = open("./output.txt", "a")
f.write("3\n")
f.close()

#--------------
#try a baseline solution:
#--------------
baseline_predictor = np.zeros(target.shape[0])
print("The MSE of the baseline prediction (all zeros) is: " + str(mean_squared_error(target, baseline_predictor)))
f = open("./output.txt", "a")
f.write("Baseline prediction: "+str(mean_squared_error(target, baseline_predictor))+"\n")
f.close()


#--------------
#train the real model
#--------------
for i in lambdas:
    clf = SGDRegressor(penalty = 'l1', alpha=i, eta0 = .00000001, l1_ratio=1, random_state = 0)
    for j in range(epochs):
        for k in range(num_steps):
            if k != num_steps - 1:
                batch_x = transformer.transform(relData[k*batch_size:(k+1)*batch_size,:])
                batch_y = target[k*batch_size:(k+1)*batch_size]
            else:
                batch_x = transformer.transform(relData[k*batch_size:relData.shape[0],:])
                batch_y = target[k*batch_size:target.shape[0]]
            clf.partial_fit(batch_x, batch_y)
        f = open("./output.txt", "a")
        f.write(str(j)+"\n")
        f.close()
        pred = clf.predict(batch_x)
        mse = mean_squared_error(batch_y, pred)
        print("The MSE of the last batch was: " + str(mse) + " for epoch " + str(j+1) + " was " + str(mse))

#--------------
#Print model evaluation statistics
#--------------
#see how many of the coefficients of the model were zero, and how many were close to zero:
params = clf.coef_
numZero = 0
total = 0
for l in params:
    total += 1
    if l == 0:
        numZero += 1
f = open("./output.txt", "a")
f.write("Prop of zero parameters: " + str(numZero/total)+"\n")
f.close()

numZero = 0
total = 0
for l in params:
    total += 1
    if np.abs(l) < .000001:
        numZero += 1
f = open("./output.txt", "a")
f.write("Prop of close to zero parameters: " + str(numZero/total)+"\n")
f.close()

#Get training statistics, mse in standardized units, mse in cm, and r-squared values
pred = clf.predict(relData)
r = r2_score(target, pred)
mse = mean_squared_error(target, pred)
mse_cm = mean_squared_error(backToCm(target, mean, std, reg, age, sex), backToCm(pred, mean, std, reg, age, sex))
print("The correlation for lambda " + str(lambdas[0]) + " was " + str(r) + " on the train set")
print("The MSE for lambda " + str(lambdas[0]) + " was " + str(mse) + " on the train set")
print("The MSE in cm for this lambda was: " + str(mse_cm))

#------------------
#Evaluate on the test set
#------------------

#First split into test and validation
from sklearn.model_selection import train_test_split
indicies = list(range(relData_test.shape[0]))
X_test, X_valid, y_test, y_valid = train_test_split(indicies, indicies, test_size=0.5, random_state=1)

#This somewhat convoluted solution allows us to convert the testing and validation back to cm
test_bool_vec = np.zeros(relData_test.shape[0], dtype = bool)
valid_bool_vec = np.zeros(relData_test.shape[0], dtype = bool)

for i in X_test:
    test_bool_vec[i] = True
for i in X_valid:
    valid_bool_vec[i] = True
X_test = relData_test[test_bool_vec,:]
X_valid = relData_test[valid_bool_vec,:]
y_test = test_target[test_bool_vec]
y_valid = test_target[valid_bool_vec]
age_test = age_test[test_bool_vec]
sex_test = sex_test[test_bool_vec]

#Get test statistics, mse in standardized units, mse in cm, and r-squared values
pred = clf.predict(X_test)
r = r2_score(y_test, pred)
mse = mean_squared_error(y_test, pred)
mse_cm = mean_squared_error(backToCm(y_test, mean, std, reg, age_test, sex_test), backToCm(pred, mean, std, reg, age_test, sex_test))
print("The correlation for lambda " + str(lambdas[0]) + " was " + str(r) + " on the test set")
print("The MSE for lambda " + str(lambdas[0]) + " was " + str(mse) + " on the test set")
print("The MSE in cm for this lambda was: " + str(mse_cm))