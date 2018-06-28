import csv
from math import sqrt
from sys import exit
import pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import matplotlib.mlab as mlab
from mpl_toolkits.basemap import Basemap
import os
from scipy.stats import norm
import matplotlib as mpl
from matplotlib.patches import Polygon
import matplotlib.cm as cm
import random

###############################################
# Functions
###############################################
def stdDev(x):
    '''function to compute standard deviation'''
    xAvg=np.mean(x)
    xOut=0.
    for k in range(len(x)):
        xOut=xOut+(x[k]-xAvg)**2
    xOut=xOut/(k+1)
    xOut=sqrt(xOut)
    return xOut

def Variance(x):
    '''function to compute the variance (std dev squared)'''
    xAvg=np.mean(x)
    xOut=0.
    for k in range(len(x)):
        xOut=xOut+(x[k]-xAvg)**2
    xOut=xOut/(k+1)
    return xOut

def SumOfSquares(x):
    '''function to compute the sum of squares'''
    xOut=0.
    for k in range(len(x)):
        xOut=xOut+x[k]**2
    return xOut

def corr(x,y):
    ''' function to find the correlation of two arrays'''
    xAvg=np.mean(x)
    Avgy=np.mean(y)
    rxy=0.
    n=min(len(x),len(y))
    for k in range(n):
        rxy=rxy+(x[k]-xAvg)*(y[k]-Avgy)
    rxy=rxy/(k+1)
    stdDevx=stdDev(x)
    stdDevy=stdDev(y)
    rxy=rxy/(stdDevx*stdDevy)
    return rxy

###############################################

wddata='/Users/lilllianpetersen/iiasa/data/'
wdfigs='/Users/lilllianpetersen/iiasa/figs/'

fdescription=open(wddata+'ETH_2004_average_household_expendetures/asceII_data/ETH_2004_E_P_data_description.csv','r')
fdata=open(wddata+'ETH_2004_average_household_expendetures/asceII_data/ASCII/ETH_2004_E_P.dat','r')

indexids=[]
indexDescriptions=[]
charecterBegin=np.zeros(shape=(33),dtype=int)
charecterEnd=np.zeros(shape=(33),dtype=int)
i=-1
for line in fdescription:
	i+=1
	line=line.replace("'",'')
	line=line.replace('\n','')
	tmp=line.split(',')
	indexids.append(tmp[0])
	charecterBegin[i]=int(tmp[1])-1
	charecterEnd[i]=int(tmp[2])
	indexDescriptions.append(tmp[3])

for i in range(len(indexids)):
	if indexids[i]=='region' or indexids[i]=='hid_p' or indexids[i]=='rururb':
		vars()[indexids[i]]=np.zeros(shape=(21595),dtype=int)
	else:
		vars()[indexids[i]]=np.zeros(shape=(21595))

totalExp=np.zeros(shape=(10,2,4672))
totalExpMask=np.ones(shape=(10,2,4672))
totalExpAvgDev=np.zeros(shape=(10,2,2))
k=-1
n=-1
for line in fdata:
	n+=1
	k+=1
	for i in range(len(indexids)):
		if i!=3 and i!=4:
			vars()[indexids[i]][k]=float(line[charecterBegin[i]:charecterEnd[i]])
	exit()

	rururb[k]=rururb[k]-1 # 0=rural, 1=urban
	region[k]=region[k]-1 #0='Tigray', 1='Afar', 2='Amhara', 3='Oromiya', 4='Somale', 5='Benishangul-Gumuz', 6='SNNPRG', 7='Harari', 8='Addis Ababa', 9='Dire Dawa'
	if region[k]!=region[k-1]:
		n=0
	totalExp[region[k],rururb[k],n]=hhexp_n[k]
	if totalExp[region[k],rururb[k],n]!=0:
		totalExpMask[region[k],rururb[k],n]=0 # 0=good

totalExp=np.ma.masked_array(totalExp,totalExpMask)

regionNames=['Tigray', 'Afar', 'Amhara', 'Oromiya', 'Somale', 'BG', 'SNNPRG', 'Harari', 'Addis Ababa', 'Dire Dawa']

for r in range(10):
	for ru in range(2):
		totalExpAvgDev[r,ru,0]=np.ma.mean(totalExp[r,ru])
		totalExpAvgDev[r,ru,1]=np.ma.std(totalExp[r,ru])

fig, ax = plt.subplots()
ax1 = ax.twinx()
ax.grid(True)

x=np.arange(10)-.18
ydata=np.ma.compressed(totalExpAvgDev[:,1,0])
ax.bar(x,ydata,width=.35,color='g',label='urban')
ax.legend(loc='upper left')

x=np.arange(10)+.18
ydata=np.ma.compressed(totalExpAvgDev[:,0,0])
ax1.bar(x,ydata,.35,color='b',label='rural')
ax1.legend(loc='upper right')

plt.xticks(np.arange(10),regionNames)
for tick in ax.get_xticklabels():
	    tick.set_rotation(30)
ax1.set_ylim([0,13000])
ax.set_ylim([0,13000])
ax1.set_yticks([])
ax.set_ylabel('Average Household Expenditures, Local Currency')

plt.title('Ethiopian Average Household Expenditures by Region')
plt.savefig(wdfigs+'eth_household_expenditures.pdf')

	


