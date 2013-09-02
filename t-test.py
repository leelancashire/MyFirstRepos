    """
    /* A numerically stable one-pass algorithm is used, see
       http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm 
       here: s ~ mean, ss ~ M2.
       Work through the large matrix x in the order in which it is in memory (column-wise) -
       in the hope that this may speed up getting it through the CPU. */
    switch(which) {
	case 0:  /* by row */
	"""

i=0
nrgrp = 2
nr = nrow
nc = ncol
nt = nc # num tests
fac = endpoint_data[:, i]

# storage for intermediate quantities
s = np.empty((2, nt))
ss = np.empty((2, nt))


x = train_data

for grp in xrange(0, nrgrp):
	for i in xrange(0, nt):
		s[grp][i] = ss[grp][i] = 0 


n = np.empty((2))
for grp in xrange(0, nrgrp):
	n[grp] = 0
for i in xrange(0, nr):
	grp = fac[i]
	n[grp]+= 1
	for j in xrange(0, nc):
		z = x[i, j]
		delta = z - s[grp][j]
                newmean = s[grp][j] + delta/n[grp]
                s[grp][j]  = newmean
                ss[grp][j] += delta*(z-newmean)
                        
# calc statistic
df = n[0]+n[1]-2;
dm = []
statistic = []
p = []
factor = np.sqrt((df) * n[0] * n[1] / (n[0]+n[1]))
for i in xrange(0, nt):
    z = ss[0][i] + ss[1][i]
    dm.append(s[0][i] - s[1][i])
    statistic.append(factor * dm[i] / np.sqrt(z))
    p.append(2*(1-stats.t.sf(statistic[i], df)))