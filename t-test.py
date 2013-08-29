"""
This script will read in some gene expression data, and perform a t-test to find differences between groups
"""

# Lee Lancashire, August 2013
 
import numpy as np
from scipy.stats import ttest_ind
import csv as csv

csv_file_object = csv.reader(open('test.csv', 'rb')) #Load in the training csv file
# contains 244 rows (244 samples) X 729 cols (728 variables + 1 output) 
header = csv_file_object.next() #Skip the fist line as it is a header
train_data=[] #Creat a variable called 'train_data'
for row in csv_file_object: #Skip through each row in the csv file
    train_data.append(row[0:]) #adding each row to the data variable

train_data = np.array(train_data) #Then convert from a list to an array

# transpose
#train_data = train_data.transpose()
# get dimensions
train_data.shape

# find data belonging to each class
group0 = train_data[:, 728] == '0'
group0 = train_data[group0, :728].astype(np.float) # convert from string (as read in from file), to float
group1 = train_data[:, 728] == '1'
group1 = train_data[group1, :728].astype(np.float)

# for each variable- performa a two-sample t-test
# null hypothesis: the two groups have the same mean
# this test assumes the two groups have the same variance...
# (can be checked with tests for equal variance)
# independent groups: e.g., how boys and girls fare at an exam
# dependent groups: e.g., how the same class fare at 2 different exams
# t_statistic, p_value = ttest_ind(group1, group2)

for x in xrange(0,728):
	t_gene1, p_gene1 = ttest_ind(group0[:, x], group1[:, x])
	print t_gene1
	print p_gene1
