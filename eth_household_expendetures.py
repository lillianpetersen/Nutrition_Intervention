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
		print indexids[i]
	else:
		vars()[indexids[i]]=np.zeros(shape=(21595))
		print indexids[i]

totalExp=np.zeros(shape=(11,2,4672))
totalExpMask=np.ones(shape=(11,2,4672))
totalExpDev=np.zeros(shape=(11,2))
avgExpRU=np.zeros(shape=(11,2))
avgExp=np.zeros(shape=(11))
regions=np.load(wdvars+'regions.npy')
regionNametoIndex={}
regionNames=['Tigray', 'Afar', 'Amhara', 'Oromia', 'Somali', 'Benishangul-Gumuz', 'SNNP', 'Harari', 'Addis Ababa', 'Dire Dawa']

i=-1
for reg in regions:
	i+=1
	regionNametoIndex[reg]=i

k=-1
n=-1
for line in fdata:
	n+=1
	k+=1
	for i in range(len(indexids)):
		if i!=3 and i!=4:
			vars()[indexids[i]][k]=float(line[charecterBegin[i]:charecterEnd[i]])

	rururb[k]=rururb[k]-1 # 0=rural, 1=urban
	region[k]=region[k]-1 # Local numbering system
	regionName=regionNames[region[k]]
	## switch to my global numbering system ##
	if region[k]!=region[k-1]:
		n=0
	ireg=regionNametoIndex[regionName]

	totalExp[ireg,rururb[k],n]=hhexp_n[k]
	if totalExp[ireg,rururb[k],n]!=0:
		totalExpMask[ireg,rururb[k],n]=0 # 0=good

totalExp=np.ma.masked_array(totalExp,totalExpMask)

perUrban=np.load(wdvars+'percentUrban.npy')
perUrban=perUrban/100.
perRural=1-perUrban
for r in range(11):
	for ru in range(2):
		avgExpRU[r,ru]=np.ma.mean(totalExp[r,ru])
		totalExpDev[r,ru]=np.ma.std(totalExp[r,ru])

avgExp[:]=(avgExpRU[:,0]*perRural)+(avgExpRU[:,1]*perUrban)

np.save(wdvars+'avgExpRURuralUrban',avgExpRU)


Mask=np.zeros(shape=(roadsByArea.shape),dtype=bool)
for i in range(len(roadsByArea)):
	if np.isnan(roadsByArea[i])==True or np.isnan(avgExpRU[i,0])==True:
		Mask[i]=True 

plt.clf()
fig, ax = plt.subplots(figsize=(12, 6))
ax.grid(True)

x=np.arange(10)
ydata1=np.ma.compressed(np.ma.masked_array(avgExpRU[:,1],Mask))
ydata2=np.ma.compressed(np.ma.masked_array(avgExpRU[:,0],Mask))
ydata3=np.ma.compressed(np.ma.masked_array(avgExp[:],Mask))
ax.bar(x-.26,ydata1,width=.25,color='g',label='urban')
ax.bar(x,ydata2,width=.25,color='b',label='rural')
ax.bar(x+.26,ydata3,width=.25,color='k',label='total')

ax.legend(loc='upper left')
#
#x=np.arange(11)+.18
#ydata=np.ma.compressed(avgExpRU[:,0])
#ax1.bar(x,ydata,.35,color='b',label='rural')
#ax1.legend(loc='upper right')

plt.xticks(np.arange(11),regions[np.invert(Mask)])
for tick in ax.get_xticklabels():
	    tick.set_rotation(20)
#ax1.set_ylim([0,13000])
ax.set_ylim([4000,13000])
#ax1.set_yticks([])
ax.set_ylabel('Average Household Expenditures, Local Currency')

plt.title('Ethiopian Average Household Expenditures by Region')
plt.savefig(wdfigs+'eth_household_expenditures_new.pdf')

###########################################################
# Correlations
###########################################################

###### Roads ######
roadsByArea=np.load(wdvars+'roadsByArea.npy')

Mask=np.zeros(shape=(roadsByArea.shape),dtype=bool)
for i in range(len(roadsByArea)):
	if np.isnan(roadsByArea[i])==True or np.isnan(avgExpRU[i,0])==True:
		Mask[i]=True 
	
Mask[0]=True
roadsByAreaM=np.ma.compressed(np.ma.masked_array(roadsByArea,Mask))
avgExpM=np.ma.compressed(np.ma.masked_array(avgExp[:],Mask))

corrRoads=corr(roadsByAreaM,avgExpM)
slope,bInt=np.polyfit(roadsByAreaM,avgExpM,1)
yfit=slope*roadsByAreaM+bInt

plt.clf()
plt.plot(roadsByAreaM,avgExpM,'*b',roadsByAreaM,yfit,'-g')
plt.xlabel('km of roads/sq. km')
plt.ylabel('Average Household Expenditures')
plt.title('Roads and Household Expenditures, Corr = '+str(round(corrRoads,2)))
plt.savefig(wdfigs+'roads_expenditures',dpi=700)

###### Flights ######
flightsByRegion=np.load(wdvars+'flightsByRegion.npy')
flightsByArea=np.load(wdvars+'flightsByArea.npy')
airportsByRegion=np.load(wdvars+'airportsByRegion.npy')
airportsByArea=np.load(wdvars+'airportsByArea.npy')
Mask[0]=True
flightsByRegionM=np.ma.compressed(np.ma.masked_array(flightsByRegion,Mask))
flightsByAreaM=np.ma.compressed(np.ma.masked_array(flightsByArea,Mask))
airportsByRegionM=np.ma.compressed(np.ma.masked_array(airportsByRegion,Mask))
airportsByAreaM=np.ma.compressed(np.ma.masked_array(airportsByArea,Mask))
avgExpM=np.ma.compressed(np.ma.masked_array(avgExp[:],Mask))

x=flightsByRegionM
corrFlights=corr(x,avgExpM)
slope,bInt=np.polyfit(x,avgExpM,1)
yfit=slope*x+bInt
plt.clf()
plt.plot(x,avgExpM,'*b',x,yfit,'-g')
plt.xlabel('Number of Originating Flights')
plt.ylabel('Average Household Expenditures')
plt.title('Flights and Household Expenditures, Corr = '+str(round(corrFlights,2)))
plt.savefig(wdfigs+'flights_expenditures',dpi=700)

x=flightsByAreaM
corrFlights=corr(x,avgExpM)
slope,bInt=np.polyfit(x,avgExpM,1)
yfit=slope*x+bInt
plt.clf()
plt.plot(x,avgExpM,'*b',x,yfit,'-g')
plt.xlabel('Number of Originating Flights by Area')
plt.ylabel('Average Household Expenditures')
plt.title('Flights by Area and Household Expenditures, Corr = '+str(round(corrFlights,2)))
plt.savefig(wdfigs+'flights_area_expenditures',dpi=700)


x=airportsByRegionM
corrFlights=corr(x,avgExpM)
slope,bInt=np.polyfit(x,avgExpM,1)
yfit=slope*x+bInt
plt.clf()
plt.plot(x,avgExpM,'*b',x,yfit,'-g')
plt.xlabel('Number of Airports')
plt.ylabel('Average Household Expenditures')
plt.title('Airports and Household Expenditures, Corr = '+str(round(corrFlights,2)))
plt.savefig(wdfigs+'airports_expenditures',dpi=700)

x=airportsByAreaM
corrFlights=corr(x,avgExpM)
slope,bInt=np.polyfit(x,avgExpM,1)
yfit=slope*x+bInt
plt.clf()
plt.plot(x,avgExpM,'*b',x,yfit,'-g')
plt.xlabel('Number of Airports by Area')
plt.ylabel('Average Household Expenditures')
plt.title('Airports by Area and Household Expenditures, Corr = '+str(round(corrFlights,2)))
plt.savefig(wdfigs+'airports_area_expenditures',dpi=700)







