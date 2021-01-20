import csv
import numpy as np
chroms = {}
for i in range(1,23):
    pVals = {}
    #Choose the file to read the p-vales from
    file_name = "DESIRED_PATH_TO_OUTPUT_"+str(i)+".assoc.linear"
    with open(file_name, 'r') as f:
        #Tab delimiter doesn't actually separate this properly, will do manually
        reader = csv.reader(f, delimiter='\t')
        #Don't read the header
        j = 0
        
        #read each row
        for row in reader:
            #Need to manually parse the row
            parsed_row = []
            #keep track of the previous character
            prev_char = ' '
            #the string we are reading
            element_string = ''
            #tell us if we are at the end of the row, need to this to make sure it reads the last column
            k = 0
            for char in row[0]:
                #If we have relevant information
                if char != ' ':
                    element_string += char
                #If we have read to the end of the relevant information
                if prev_char != ' ' and char == ' ':
                    parsed_row.append(element_string)
                    element_string = ''
                #If we are at the end of the line and we have unappended data
                if k == len(row[0])-1 and prev_char != '':
                    parsed_row.append(element_string)
                k+=1
                prev_char = char
            #if not in the header
            if j != 0:
                pVals[parsed_row[1]] = np.float64(parsed_row[8])
            j = 1
    f.close()
    chroms[i] = pVals

#create a sorted list of p-values to quickly find correct SNPs
sortedList = []
for i in chroms.keys():
    sortedList.append(sorted(chroms[i].values()))
import progressbar
#get the lowest N p-value SNP names, separated by chromosome:
N = 50000
#list of the top N SNPs
topNSNPs = {}
for i in chroms.keys():
    topNSNPs[i] = []

#Find N SNPs that are the top of each list
with progressbar.ProgressBar(max_value=N) as bar:
    for i in range(N):
        #keep track of the index of the minimum value (this just finds what chromosome has the next min SNP value)
        minInd = 0
        #For each chromosome
        for j in range(1,len(sortedList)):
            #compare sortedList of chromosome J at the 0th index (this will be the min in the list)
            if sortedList[j][0] < sortedList[minInd][0]:
                minInd = j
        #Find the SNP name
        #Note: the chromosome number in the normal naming scheme is one greater than the index we found the min p-val at
        for j in chroms[minInd+1].keys():
            if chroms[minInd+1][j] == sortedList[minInd][0]:
                topNSNPs[minInd+1].append(j)
                del chroms[minInd+1][j]
                break
        del sortedList[minInd][0]
        #The chromosome of the minInd must now be iterated since we have now gotten that SNP
        listOfInd[minInd] += 1
        bar.update(i)

"""
Verify two facts about the output: N SNPs were found and no repeated SNPs
"""
verifyN = 0
for i in topNSNPs.keys():
    verifyN += len(topNSNPs[i])
print("The number of SNPs found is: " + str(verifyN))

isDistinct = 0
for i in topNSNPs.keys():
    if len(topNSNPs[i]) == len(list(np.unique(topNSNPs[i]))):
        isDistinct += 1
if isDistinct == len(topNSNPs.keys()):
    print("All SNPs found are distinct")

#Write everything to files:
for i in topNSNPs.keys():
    file1 = open("PATH_TO_RELEVANT_SNPS_"+str(i)+".txt","w")
    for j in range(len(topNSNPs[i])):
        string = topNSNPs[i][j] + "\n"
        file1.write(string)
    file1.close()