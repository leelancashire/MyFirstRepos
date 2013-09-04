# -*- coding: utf-8 -*-
"""
This script will read in some gene expression data, then some mutation data, and 
search for genes whose mutation status may cause drastic over or under expression of a list of pathways

"""

# Lee Lancashire, August 2013
 
import numpy as np
#from scipy.stats import ttest_ind
from scipy.stats import fisher_exact
import csv as csv
import heatmap # heatmap.py file
## R stuff
import rpy2.robjects as robjects
r = robjects.r
# allow conversion of numpy objects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from rpy2.robjects.packages import importr
utils = importr("utils")
rstats = importr("stats")
limma = importr("limma")
surv = importr("survival")
#gf = importr("genefilter")

# R commands can be used with a preceding print, eg length, etc..
# print(fl)
# print(r.length(fl))

""" Read Data """

# Read input data
csv_file_object = csv.reader(open('/Users/llancashire/Dropbox/ThomsonReuters/Current Projects/omics integration/source_dat/gene_dat.csv', 'rb')) #Load in the training csv file
# contains 244 rows (244 samples) X 12042 cols  
train_header = csv_file_object.next() #Skip the fist line as it is a header
geneNames = robjects.StrVector(train_header)
train_data=[] #Creat a variable called 'train_data'
for row in csv_file_object: #Skip through each row in the csv file
    train_data.append(row[0:]) #adding each row to the data variable- start at 1 if you want to avoid rownames

train_data = np.array(train_data) #Then convert from a list to an array
train_data = train_data.astype(np.float) # convert from string to float

# Read endpoint data
csv_file_object = csv.reader(open('/Users/llancashire/Dropbox/ThomsonReuters/Current Projects/omics integration/source_dat/mut_dat.csv', 'rb')) #Load in the endpoint csv file
# contains 244 rows (244 samples)
endpoint_header = csv_file_object.next() #Skip the fist line as it is a header
endpoint_data=[] #Creat a variable called 'train_data'
for row in csv_file_object: #Skip through each row in the csv file
    endpoint_data.append(row[0:]) #adding each row to the data variable

endpoint_data = np.array(endpoint_data) #Then convert from a list to an array
endpoint_data = endpoint_data.astype(np.float) # convert from string to float

# read clinical data
csv_file_object = csv.reader(open('/Users/llancashire/Dropbox/ThomsonReuters/Current Projects/TCGA/Processed Data/survival_dat.csv', 'rb'))
surv_header = csv_file_object.next() #Skip the fist line as it is a header
surv_data=[] #Creat a variable called 'train_data'
for row in csv_file_object: #Skip through each row in the csv file
    surv_data.append(row[0:]) #adding each row to the data variable

surv_data = np.array(surv_data) #Then convert from a list to an array
surv_data = surv_data.astype(np.int) # convert from string to float

# transpose
# train_data = train_data.transpose()
# get dimensions
train_data.shape
nrow = train_data.shape[0]
ncol = train_data.shape[1]


# standardize training data; mean = 0, std = 1
for x in xrange(0, ncol):
	mean = np.mean(train_data[:, x])
	std = np.std(train_data[:, x])
	train_data[:, x] = (train_data[:, x] - mean) / std

# set up data
m = robjects.r['matrix'](train_data.transpose(), nrow = ncol) # data must be genes in rows, samples in cols, for limma

""" Loop through each mutation, finding the DEGs associated with mutation using limma"""

for i in xrange(0, len(endpoint_header)):
	print "Variable: " + str(i) + " of " + str(ncol) 
	# endpoint
	endpoint = endpoint_data[:, i]
	group0 = endpoint == 0
	# ensure we have at least 2 samples in each class
	num_grp0 = sum(group0 == True)
	num_grp1 = sum(group0 == False)
	if num_grp0 > 1 and num_grp1 > 1:
	    fl = robjects.FactorVector(endpoint) # factor for R limma
	    robjects.globalenv["description"] = fl
	    fmla = robjects.Formula('~ description + 0')
            design = rstats.model_matrix(fmla)  
            design.colnames = robjects.StrVector(['Norm','Mut']) 
            # print(design)
            # robjects.globalenv["design"] = design
            fit = limma.lmFit(m, design)
            contMat = robjects.IntVector([-1, 1])   
            fit2 = limma.contrasts_fit(fit, contMat)
            fit2 = limma.eBayes(fit2)   
            corrGenes = limma.decideTests(fit2, adjust_method='fdr', p_value=0.01)
            tT = limma.topTable(fit2, adjust='fdr', sort_by="B", number=ncol, genelist=geneNames)
            # print(r.head(tT))
            # loop through corrGenes, and find the DEGs
            DEGs = []
            for x in xrange(0, ncol):
	      	if corrGenes[x] != 0.0:
		    DEGs.append(x)
            # create expression matrix just containing DEGs
	    DEGs_data = train_data[:, DEGs]
	    numDEGs = DEGs_data.shape[1]
	    # create 2 binary expression matrices; 1. Up-regulated genes; 2. Down-regulated genes
  	    # For 1; if z-score > 2, value = 1, 0 otherwise
  	    # For 2; if z-score < -2, value = 1, 0 otherwise
  	    up = np.array(DEGs_data, copy=True)
	    up[up <2] = 0.0
	    up[up >2] = 1.0
	    down = np.array(DEGs_data, copy=True)
	    down[down >-2] = 0.0
	    down[down <-2] = 1.0
	    # for each gene's up and down table, calculate fisher exact
	    #Given a 2x2 table:
 	    #+---+---+
            #| a | b |
            #+---+---+
            #| c | d |
            #+---+---+
 	    # represented by a list of lists::
 	    # [[a,b],[c,d]]
 	    #
    	    #this calculates the sum of the probability of this table and all others
    	    #more extreme under the null hypothesis that there is no association between
    	    #the categories represented by the vertical and horizontal axes.
    	    #a = IS differentially expressed and IS mutated
    	    #b = NOT differentially expressed and IS mutated
    	    #c = IS differentially expressed and NOT mutated
    	    #d = NOT differentially expressed and NOT mutated
    	    #
    	    # storage
	    fisher_res_up = []
	    fisher_res_down = []
	    sig_genes_up = []
	    sig_genes_down = []
	    for x in xrange(0, numDEGs):
		  a_up = 0
	          b_up = 0
		  c_up = 0
		  d_up = 0
		  a_down = 0
		  b_down = 0
		  c_down = 0
		  d_down = 0
		  for y in xrange(0,len(endpoint)):
			 if up[y, x] == 1 and endpoint[y] == 1.0:
				    a_up += 1
			 elif up[y, x] == 0 and endpoint[y] == 1.0:
				    b_up += 1
			 elif up[y, x] == 1 and endpoint[y] == 0.0:
				    c_up += 1
			 else:
				    d_up +=1
		  for y in xrange(0,len(endpoint)):
			 if down[y, x] == 1 and endpoint[y] == 1.0:
				    a_down += 1
			 elif down[y, x] == 0 and endpoint[y] == 1.0:
				    b_down += 1
			 elif down[y, x] == 1 and endpoint[y] == 0.0:
				    c_down += 1
			 else:
				    d_down +=1
		  # calc fisher p-val
		  fisher_up = fisher_exact([[a_up, b_up],[c_up, d_up]])[1]
		  fisher_down = fisher_exact([[a_down, b_down],[c_down, d_down]])[1]
		  if fisher_up < 0.05/numDEGs:
			 sig_genes_up.append(x)
		  if fisher_down < 0.05/numDEGs:
			 sig_genes_down.append(x)
		  fisher_res_up.append(fisher_up)
		  fisher_res_down.append(fisher_down)
	   # Heatmap time
	   # generate a heatmap for the significant genes
	   # get list of sig genes
	    s = [DEGs[z] for z in sig_genes_up]
	    s.extend(DEGs[z] for z in sig_genes_down)
	   #s = sig_genes_up
	   #s.extend(sig_genes_down)
	    if len(s) > 1:
	         # create matrix
		 dat = train_data[:, s]
		# create vector to label as mutated or non-mutated
		 mut_status = []
		 for x in xrange(0,len(endpoint)):
		    if endpoint[x] == 0:
		        mut_status.append("-")
		    else:
			mut_status.append("MUTATION")
		# create vector of biomarker names
		 col_names = [train_header[z] for z in s]
		 """heatmap(x, row_header, column_header, row_method,
	        column_method, row_metric, column_metric,
	        color_gradient, filename)
		"""
		 fname = 'gene-heatmaps/' + endpoint_header[i] + '.png'
		 heatmap.heatmap(x = dat.transpose(), row_header = col_names, column_header = mut_status, \
	    	 row_method = 'ward', column_method = 'ward', row_metric = 'euclidean', column_metric = 'euclidean', \
			 color_gradient = 'red_white_blue', filename = fname)
	       # Survival analysis time
	       # Calculate survival scores the significant genes
	         time = robjects.IntVector(surv_data[:, 0])
	         event = robjects.IntVector(surv_data[:, 1])
	         for x in xrange(0, len(s)):
	             y = [getGrps(dat[l, x], endpoint[l]) for l in range(len(endpoint))]
	         y = robjects.StrVector(y)
	         fit1 = surv.Surv(time, event)
	         robjects.globalenv["fit1"] = fit1
	         robjects.globalenv["y"] = y
                 fmla = robjects.Formula("fit1 ~ y")
	         fit2 = surv.survfit(fmla)
	             
	         r.plot(fit2, lty=1:4, col=1:4, conf_int = FALSE, xlab="Survival in Days", ylab = "Proportion Survived")
                r.legend(1500, .9, names(fit$strata), lty = 1:4, col=1:4 )    
	             
	             
	             fmla = robjects.Formula("time, event ~ y")
	             surv_data[:, 0], surv.data[:, 1]) ~ y
	       
	       

fit <- survfit(Surv(time, event) ~ y)


# function to split vector into 4 groups, and return the group
def getGrps(gene, mut):
    if gene > 0 and mut == 1:
        grp = "Mutated with high expression"
    elif gene > 0 and mut == 0:
        grp = "Not Mutated with high expression"
    elif gene < 0 and mut == 1:
        grp = "Mutated with low expression"
    else:
        grp = "Not Mutated with low expression"
    return grp


"""
t-test code
# find data belonging to each class
	group0 = endpoint == 0
	# ensure we have at least 2 samples in each class
	num_grp0 = sum(group0 == True)
	num_grp1 = sum(group0 == False)
	if num_grp0 > 1 and num_grp1 > 1:
		group0 = train_data[group0, :]
		group1 = endpoint == 1
		group1 = train_data[group1, :]
		# for each variable- performa a two-sample t-test
		# tstatistic, p_value = ttest_ind(group1, group2)
		# set up storage for pvalues
		tscore = []
		pval = []
		DEGs = []
		for x in xrange(0, ncol):
			t_gene, p_gene = ttest_ind(group0[:, x], group1[:, x])
			tscore.append(t_gene)
			pval.append(p_gene)
			if p_gene < (0.05/ncol): # find indeces of genes with p-val < 0.05/n (boneferroni correction)
				DEGs.append(x)
"""
				
