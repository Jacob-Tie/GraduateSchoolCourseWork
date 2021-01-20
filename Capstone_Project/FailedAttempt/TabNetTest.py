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

print("Is there a GPU?")
print(torch.cuda.is_available())

#Load the data as a dictionary
def loadDataDict():
    output = {}
    for i in range (1, 23):
        output[str(i)] = zarr.open('PATH_TO_ZARR_'+str(i), mode='r')
    return output
data = loadDataDict()

"""
This function will output a dictionary containing what SNPs are contained in what Genes further organized by chromosome

Returns: a dictionary that leads to other dictionaries. The first key will be the chromosome, then this leads to a dictionary of genes -> SNPs
"""
def GetGeneConnections():
    result = {}
    #iterate through chromosomes
    for i in range(1,23):
        print("Working on loading connections from chromosome " + str(i))
        connections = {}
        #open the log file we created using magma
        with open('MAGMA_OUTPUT_PATH_' + str(i) + '.genes.annot', 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            #j will allow us to skip the first two rows of the file (which we won't need)
            j = 0
            #read file line by line
            for row in reader:
                #Flag is for debugging purposes
                flag = 0
                if j > 1:
                    temp = []
                    #iterate through the row
                    for k in range(len(row)):
                        #first two entries are meta data (not required)
                        if k > 1:
                            temp.append(row[k])
                    #open the file with gene -> SNP relations
                    with open("PATH_TO_MAGMA_GENE_LOC", 'r') as fp:
                        reader2 = csv.reader(fp, delimiter='\t')
                        #read this line by line
                        for rows in reader2:
                            #row[0] (from our log file) corresponds to a gene number, and rows[0] (from the second file) should also be gene number
                            if rows[0] == row[0]:
                                #rows[5] contains the gene name (to be used by gene ontology) 
                                connections[rows[5]] = temp
                                flag = 1
                                break
                    fp.close()
                    #debugging
                    if flag == 0:
                        print("Something went wrong: gene name not found")
                j+=1
        f.close()
        #add the connections in the chromosome to our result dictionary
        result[str(i)] = connections
    return result
connections = GetGeneConnections()
print("Verifying each chromosome was found in our connections dictionary:")
print(connections.keys())
print("Verifying a chromosome contains genes:")
print("Number of Genes found in Chromosome 1: " + str(len(connections['1'].keys())))

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
"""
This function will sort all of the list of SNPs in our dictionary of connections

@param connections: a dictionary created by GetGeneConnections()
Returns: a dictionary in the same format as GetGeneConnections() outputs, but with all of the SNP lists sorted
"""
def sortAllSNPs(connections):
    result = {}
    for i in connections.keys():
        print("Sorting chromosome " + i + "s SNPs")
        temp = connections[i]
        for j in temp.keys():
            temp[j] = sorted(temp[j])
        result[i] = temp
    return result
connections = sortAllSNPs(connections)

"""
This function will find the indicies of relavent SNPs in a dataset:
@param relSNPs: a sorted list of strings corresponding to the SNPs that we want to find
@param data: a dictionary containing the data organized by chromosome (the key "1" should be the data from the first chromosome, etc)

Returns: a dictionary with the same keys as data that correspond to the indicies in that part of the data with the relevant SNPs
"""
def findIndOfRelSNPs(relSNPs, data):
    #Initialize dictionary for result
    result = {}
    #for each key in data
    for chrom in data:
        #tempData will store a list of SNP names in [chrom] key of the data
        tempData = data[chrom].snp
        #list of indicies of SNPs found in 'chrom'
        foundInd = []
        #for all the SNPs in this slice of the data
        for i in range(0, len(tempData)):
            #preform binary search on relSNPs
            boolFound = BinarySearch(relSNPs, tempData[i])
            #if ith snp in tempData is found in relSNPs append i
            if boolFound != -1:
                foundInd.append(int(i))
        #add what we found to the result with the same key as data
        result[chrom] = foundInd
    return result
"""
This function will subset a zarr array into a numpy array.
@param relInd: a dictionary containing the relevent indicies from each chromosome (the indicies we are intersted in)
@param data: a dictionary of zarr arrays containing the complete dataset
@param n: the first row of the zarr array to return
@param m: the number of rows to return
Returns: a dictionary of numpy arrays that are subsetted as described above
"""
# V2
def fastSubsetData(relInd, data, n, m):
    result = {}
    for chrom in data:
        #temp = np.zeros(data[chrom].genotype.shape[1], dtype = bool)
        #temp[relInd[chrom]] = True
        if len(relInd[chrom]) != 0:
            result[chrom] = data[chrom].genotype.get_orthogonal_selection((list(range(n,m+n)), relInd[chrom]))
        else:
            result[chrom] = np.array([-10])
    return result 

"""
This function will take a dictionary whose keys point to data and it will return a single dataset that represents the concatination of all the data

@param relData: the dictionary containing the data

Returns: a numpy array that has all the data joined column wise starting with the first key in the dictionary
"""
def concatData(relData):
    rows = 0
    cols = 0
    for chrom in relData:
        if np.all(relData[chrom] != -10):
            rows = relData[chrom].shape[0]
            cols += relData[chrom].shape[1]
    output = np.zeros((rows,cols))
    colStart = 0
    for chrom in relData:
        if np.all(relData[chrom] != -10):
            output[:, int(colStart):colStart+int(relData[chrom].shape[1])] = relData[chrom]
            colStart += relData[chrom].shape[1]
    return output

"""
Gets target values from a file that we have previously created

@param n: the number of values to get
@param m: starting value
"""
#Note that this should be replaced by a better function to read target files
def fastGetTarget(n,m):
    #read in the relevant SNP data
    intNaVals = []
    target = []
    with open("PATH_TO_TARGET") as f:
        #create an iterable object, NOTE: the file is spaced by tabs not commas
        reader=csv.reader(f)
        i = 0
        flag = 0
        #iterate through the file
        for row in reader:
            temp = str(row)
            temp = temp[2:]
            temp = temp[:-2]
            if i >= m:
                if flag == 1 and temp != "==========":
                    intNaVals.append(int(temp))
                if flag == 2 and temp != "==========":
                    target.append(float(temp))
                if temp == "==========":
                    flag += 1
                if flag == 3:
                    break
                    
            i+=1
    f.close()
    target = np.array(target)
    return target, intNaVals

def eliminateNaRows(relData, intNaVals):
    for chrom in relData:
        temp = np.ones(relData[chrom].shape[0], dtype = bool)
        temp[intNaVals] = False
        relData[chrom] = relData[chrom][temp,:]
    return relData

"""
This function will load the whole Dataset into RAM from a given subset of columns.

@param data: the Zarr file (or more accuratly a dictionary of Zarr files) that we want to load into memory
@param relInd: a dictionary containing the columns of each Zarr data that we want to load
@param oneHot: an instance of the myOHE class which will specify how to preprocess the data

Returns: the whole dataset loaded as 8 bit integers
"""
def loadWholeDataset(data, relInd):
    chunk = 10000
    n = data['1'].genotype.shape[0]
    steps = int(np.ceil(n/chunk))  
    for i in range(0,steps):
        if i == 0:
            relData = fastSubsetData(relInd, data, chunk*i, chunk)
            relData = concatData(relData)
            relData = relData.astype("int8")
            outputData = relData
            
        elif chunk*(i+1) < n:
            relData = fastSubsetData(relInd, data, chunk*i, chunk)
            relData = concatData(relData)
            relData = relData.astype("int8")
            outputData = np.vstack((outputData, relData))
            
        else:
            relData = fastSubsetData(relInd, data, chunk*i, n-chunk*i)
            relData = concatData(relData)
            relData = relData.astype("int8")
            outputData = np.vstack((outputData, relData))
    return outputData

# For reproducability
torch.manual_seed(1)
"""
Define a series of blocks of varying size (for experimentation) to see if bigger networks sitting between the SNP and gene layers are important
"""
class smallBlock(nn.Module):
    def __init__(self, nInput):
        super(smallBlock, self).__init__()
        self.nInput = nInput
        self.fc1 = nn.Linear(self.nInput, 1)
        self.leReLU = nn.LeakyReLU(.1)
    def forward(self, x):
        x = self.leReLU(self.fc1(x))
        return x
    
class mediumBlock(nn.Module):
    def __init__(self, nInput):
        super(mediumBlock, self).__init__()
        self.nInput = nInput
        self.fc1 = nn.Linear(self.nInput, 10)
        self.fc2 = nn.Linear(10, 1)
        self.leReLU = nn.LeakyReLU(.1)
    def forward(self, x):
        x = self.leReLU(self.fc2(self.leReLU(self.fc1(x))))
        return x
    
class largeBlock(nn.Module):
    def __init__(self, nInput):
        super(largeBlock, self).__init__()
        self.nInput = nInput
        self.fc1 = nn.Linear(self.nInput, 100)
        self.fc2 = nn.Linear(100, 50)
        self.fc3 = nn.Linear(50, 25)
        self.fc4 = nn.Linear(25, 10)
        self.fc5 = nn.Linear(10, 1)
        self.leReLU = nn.LeakyReLU(.1)
    def forward(self, x):
        x = self.leReLU(self.fc1(x))
        x = self.leReLU(self.fc2(x))
        x = self.leReLU(self.fc3(x))
        x = self.leReLU(self.fc4(x))
        x = self.leReLU(self.fc5(x))
        return x
#See the readme for a description of what is happening here
class Net(nn.Module):
    def __init__(self, connections, batch_size):
        super(Net, self).__init__()
        self.connections = connections
        self.SNPtoGene = nn.ModuleList()
        self.numSmall = 0
        self.numMed = 0
        self.numLarge = 0
        for i in self.connections.keys():
            for j in self.connections[i].keys():
                n = len(self.connections[i][j])
                if n <= 10:
                    self.SNPtoGene.append(smallBlock(n))
                    self.numSmall += 1
                elif n <= 100:
                    self.SNPtoGene.append(smallBlock(n))
                    self.numMed += 1
                else:
                    self.SNPtoGene.append(smallBlock(n))
                    self.numLarge += 1
        self.numGenes = 0
        for i in self.connections.keys():
            self.numGenes += len(self.connections[i].keys())
        self.tab = TabNetNoEmbeddings(self.numGenes, 1, n_steps = 4)
        self.batch_size = batch_size
        self.leReLU = nn.LeakyReLU(.2)
    def forward(self, x):
        start_gene = 0
        start_SNP = 0
        currentModel = 0
        temp = torch.zeros(self.batch_size, self.numGenes)
        temp = temp.cuda()
        for i in self.connections.keys():
            for j in self.connections[i].keys():
                n = len(self.connections[i][j]) 
                temp[: , start_gene:start_gene+1] = self.leReLU(self.SNPtoGene[currentModel](x[:, start_SNP:start_SNP + n]))
                start_gene += 1
                start_SNP += n
                currentModel += 1
        x = temp
        x = self.tab(x)
        return(x)
    def printStats(self):
        print("The number of small models is " + str(self.numSmall))
        print("The number of medium models is " + str(self.numMed))
        print("The number of large models is " + str(self.numLarge))
        print("The number of genes is " + str(self.numGenes)) 
    def getMasks(self, x):
        return self.tab.forward_masks(x)
    def getNumGenes(self):
        return self.numGenes

def findIndOfRelSNPDict(relSNPs, data):
    #Initialize dictionary for result
    result = {}
    SNPID = {}
    #for each key in data
    for chrom in data:
        print("Finding relevent indicies for Chromosome " + chrom)
        #tempData will store a list of SNP names in [chrom] key of the data
        tempData = data[chrom].snp
        tempDataSort = sorted(tempData)
        mapping = {}
        for i in range(0, len(tempData)):
            #found is the index in the sorted array that the ith element of tempData was found in
            found = BinarySearch(tempDataSort, tempData[i])
            if found == -1:
                print("something went wrong with binary search")
            else:
                mapping[found] = i
        #list of indicies of SNPs found in 'chrom'
        foundInd = []
        SNPs = []
        #for all the SNPs in this slice of the data
        for i in relSNPs[chrom].keys():
            for j in relSNPs[chrom][i]:
                found = BinarySearch(tempDataSort, j)
                if found != -1:
                    foundInd.append(found)
                    SNPs.append(j)
        for i in range(0, len(foundInd)):
            foundInd[i] = mapping[foundInd[i]]
        #add what we found to the result with the same key as data
        result[chrom] = foundInd
        SNPID[chrom] = SNPs
    return result, SNPID

relInd, SNPs = findIndOfRelSNPDict(connections, data)

#relData = fastSubsetData(relInd, data, 0, 10000)
#target, intNaVals = fastGetTarget(10000, 0)
#relData = eliminateNaRows(relData, intNaVals)

#relData = concatData(relData)
#relData = relData.astype('int8')
#relData, target = loadWholeDataset(data, relInd)
#Load data from dictionaries:
def readCovariates():
    return np.load('PATH_TO_AGE_COV',allow_pickle='TRUE').item(), np.load('PATH_TO_GEND_COV',allow_pickle='TRUE').item()
age, sex = readCovariates()

#Convert the strings of age and sex to numbers:
def convertDictValsToNum(x):
    for i in x.keys():
        x[i] = int(x[i])
    return x
age = convertDictValsToNum(age)
sex = convertDictValsToNum(sex)

target = np.load('PATH_TO_TARGET',allow_pickle='TRUE').item() #Has 'NA' values

naIndices = []
numValid = 0
for i in range(data['1'].genotype.shape[0]):
    if target[i] == 'NA':
        naIndices.append(i)
    else:
        numValid +=1
        target[i] = float(target[i])

nonNaTarget = torch.zeros(numValid)
j = 0
for i in range(data['1'].genotype.shape[0]):
    if target[i] != 'NA':
        nonNaTarget[j] = target[i]
        j+=1

mask = torch.ones(data['1'].genotype.shape[0], dtype = bool)
mask[naIndices] = False
def convertDictToNp(x):
    out = np.zeros(len(x.keys()))
    for i in x.keys():
        out[i] = x[i]
    return out
age = convertDictToNp(age)
sex = convertDictToNp(sex)
age = age[mask]
sex = sex[mask]

#exclude rows with low numbers of people
relRows = torch.ones(age.shape[0], dtype = bool)
for i in torch.unique(torch.from_numpy(age)):
    temp = torch.zeros(age.shape[0], dtype = bool)
    if i < 40 or i > 70:
        for j in range(age.shape[0]):
            if age[j] == i:
                relRows[j] = False
age = age[relRows]
sex = sex[relRows]
nonNaTarget = nonNaTarget[relRows]

#Preform the adjustments
for i in torch.unique(torch.from_numpy(sex)):
    temp = torch.zeros(age.shape[0], dtype = bool)
    for j in range(sex.shape[0]):
        if sex[j] == i:
            temp[j] = True
    mean = torch.mean(nonNaTarget[temp]).item()
    std = torch.std(nonNaTarget[temp])
    nonNaTarget[temp] = (nonNaTarget[temp]-mean)/std

from sklearn.linear_model import LinearRegression
reg = LinearRegression().fit(age.reshape(-1,1), nonNaTarget)

for i in range(nonNaTarget.shape[0]):
    temp = reg.predict(age[i].reshape(-1,1))
    nonNaTarget[i] = nonNaTarget[i] - temp

relData = loadWholeDataset(data, relInd)
relData = relData[mask,:]
relData = relData[relRows, :]
target = nonNaTarget
del(mask)
del(relRows)
print("We have successfully loaded the data")


from sklearn.model_selection import train_test_split
relData, X_test, nonNaTarget, y_test = train_test_split(relData, nonNaTarget, test_size=0.05, random_state=0)
X_test, X_valid, y_test, y_valid = train_test_split(X_test, y_test, test_size = .5, random_state = 0)

from sklearn.preprocessing import Normalizer
transformer = Normalizer().fit(relData)


batch_size = 2**10
model = Net(connections, batch_size)
if torch.cuda.is_available():
    model.cuda()
MSEloss = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=.02)
maxEpochs = 5

epoch = 0
batchesPerEpoch = int(np.ceil(relData.shape[0]/batch_size))
while epoch < maxEpochs:
    start=time.time()
    running_loss = 0
    temp = torch.randperm(relData.shape[0])
    for k in range(batchesPerEpoch-1):
        optimizer.zero_grad()
        if (k+1)*batch_size < relData.shape[0]:
            batch = torch.from_numpy(transformer.transform(relData[k*batch_size:(k+1)*batch_size, :]))
            batch_target = target[k*batch_size:(k+1)*batch_size]
        else:
            batch = torch.from_numpy(transformer.transform(relData[k*batch_size:relData.shape[0], :]))
            batch_target = target[k*batch_size:relData.shape[0]]
        if torch.cuda.is_available():
            batch = batch.cuda()
            batch_target = batch_target.cuda()
        output = model(batch.float())
        loss = MSEloss(output[0].float(), batch_target.float().view(output[0].shape))+.001*output[1]
        clip_grad_norm_(model.parameters(), 1)
        loss.backward()
        optimizer.step()
        running_loss += loss.item() - .001*output[1]
        if k == batchesPerEpoch - 2:
            print("The loss for epoch " + str(epoch) + " was " + str(running_loss/k))
    epoch += 1
    if epoch == maxEpochs:
        if k == batchesPerEpoch - 1:
            print("The final mean squared error was " + str(loss.item()))
print("The run time of this test was: " + str(time.time()-start))

del(relData)
del(target)
torch.save(model.state_dict(), "./TabNetV2")

model.eval()
numTestBatches = int(np.ceil(X_test.shape[0]/batch_size))
pred = torch.zeros((numTestBatches-1)*batch_size)
running_loss = 0
feature_importances_ = torch.zeros(model.getNumGenes())
#print("Number of Genes:")
#print(feature_importances_.shape)
for k in range(numTestBatches-1):
    batch = torch.from_numpy(transformer.transform(X_test[k*batch_size:(k+1)*batch_size,:]))
    batch_target = y_test[k*batch_size:(k+1)*batch_size]
    batch = batch.cuda()
    batch_target = batch_target.cuda()
    
    output = model(batch.float())
    loss = MSEloss(output[0].float(), batch_target.float().view(output[0].shape))
    running_loss += loss.item()
    pred[k*batch_size:(k+1)*batch_size]=output[0].view(batch_size).cpu()
    #x,y = model.getMasks(batch)
    #x = x.cpu()
    #feature_importances_ += x.sum(dim=0)
    if k == batchesPerEpoch - 2:
            print("The loss for the test set was " + str(running_loss/k))

pred = pred.detach().cpu().numpy()
print("The correlation of our predictions with the ground truth is: " + str(np.corrcoef(pred, test_y[0:(numTestBatches-1)*batch_size])))

#Feature Attribution:
#feature_importances_ = feature_importances_ / torch.sum(feature_importances_)
print("features:")
print(feature_importances_.shape)