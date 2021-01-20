import csv

#Get all of the iids in the data
iids = []
with open("PATH_TO_UKBIOBANK", 'r') as f:
    reader = csv.reader(f, delimiter=' ')
    j = 0
    for row in reader:
        iids.append(row[0])
f.close()

from sklearn.model_selection import train_test_split
import numpy as np
#split the iids into test and train sets (Note that the test set will eventually be split further into test and validation)
X_train, X_test = train_test_split(np.array(iids), test_size=0.1, random_state=42)

#Create a file with all the the individuals in the train set (repeated twice so plink can read it)
file1 = open("PATH_TO_INDIV_IN_TRAIN_SET.txt","w") 
for i in range(X_train.shape[0]):
    string = X_train[i] + " " + X_train[i] + "\n"
    file1.write(string)
file1.close()

#Create a file with all the the individuals in the test set (repeated twice so plink can read it)
file1 = open("PATH_TO_INDIV_IN_TEST_SET.txt","w") 
for i in range(X_test.shape[0]):
    string = X_test[i] + " " + X_test[i] + "\n"
    file1.write(string)
file1.close()