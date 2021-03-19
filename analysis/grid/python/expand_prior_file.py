import numpy as np
import sys

prior = np.loadtxt(sys.argv[1],unpack=True)
try:
    x,y,m,r,e1,e2,f,n,zb = prior
except:
    x,y,m,r,e1,e2,f = prior

new_prior = np.empty((prior.shape[0]+2,prior.shape[1]))
#new_prior = prior.copy()
new_prior[:-2,:] = prior

for i in range(len(prior)):
	if r[i]==0 and f[i]==0 and e1[i]==0 and e2[i]==0:
		if first: 
			first = False
			first_non_dupl = i-1000
		new_prior[2:-2,i] = prior[2:,first_non_dupl+m[i]]
	else:
		first = True

new_prior[-2,:] = new_prior[-3,:]
new_prior[-1,:] = np.linspace(0,len(x)-1,len(x)).transpose()

if prior.shape[0]==7:
    np.savetxt(sys.argv[2],zip(*new_prior),['%8.2f','%8.2f','%5.2f','%5.3f','%6.3f','%6.3f','%4.2f', '%4.2f', '%4.2f'],delimiter=' ')
else:
    np.savetxt(sys.argv[2],zip(*new_prior),['%8.2f','%8.2f','%5.2f','%5.3f','%6.3f','%6.3f','%4.2f','%4.2f','%4.2f','%4.2f','%d'],delimiter=' ')
