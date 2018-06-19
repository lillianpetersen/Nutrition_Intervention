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

f=open(wddata+'Africa_conflictss.csv','r')

monthDict={'January':0,'February':1,'March':2,'April':3,'May':4,'June':5,'July':6,'August':7,'September':8,'October':9,'November':10,'December':11}

allDataIDs=np.zeros(shape=(161939),dtype=int)
timePrecision=np.zeros(shape=(161939))
eventTypeDict={}
eventTypeCounter=0
countriesAll=[]
countryCodeAll=[]
latAll=np.zeros(shape=(161939))
lonAll=np.zeros(shape=(161939))
fatalitiesAll=np.zeros(shape=(161939))
fatalitiesOver100=0

i=-2
g=-1
for line in f:
	i+=1
	if i==-1:
		continue
	line=line.replace('\n','')
	tmp=line.split(';')
	if len(tmp)!=27:
		sourceBegins=[]
		sourceEnds=[]
		for k in range(len(tmp)):
			if len(tmp[k])>0:
				if tmp[k][0]=='"' and tmp[k][-1]!='"':
					sourceBegins.append(k)
				if tmp[k][0]!='"' and tmp[k][-1]=='"':
					sourceEnds.append(k)
		for k in range(len(sourceBegins)):
			tmp[sourceBegins[k]]=','.join(tmp[sourceBegins[k]:sourceEnds[k]+1])
		del tmp[sourceBegins[-1]+1:sourceEnds[-1]+1]
		if len(sourceBegins)>1:
			del tmp[sourceBegins[-2]+1:sourceEnds[-2]+1]
		if len(sourceBegins)>2:
			del tmp[sourceBegins[-3]+1:sourceEnds[-3]+1]
	if len(tmp)!=27:
		print 'len(tmp)!=27'
		exit()

	dataID=int(tmp[0])
	allDataIDs[g]=dataID
	if dataID==allDataIDs[g-1]:
		continue
	g+=1

	date=tmp[4]
	date=date.replace('"','')
	date=date.split(' ')
	day=int(date[0])
	month=monthDict[date[1]]
	year=int(date[2])

	timePrecision[g]=tmp[6]
	eventType=tmp[7].replace('"','')

	try:
		eventTypeDict[eventType]
	except:
		eventTypeDict[eventType]=eventTypeCounter
		eventTypeCounter+=1

	country=tmp[13]
	country=country.replace('"','')
	countryCode=tmp[26][:-1]
	countriesAll.append(country)
	countryCodeAll.append(countryCode)
	
	lat=float(tmp[18])
	lon=float(tmp[19])
	latAll[g]=lat
	lonAll[g]=lon
	geoAccuracy=tmp[20] #Probably not important

	notes=tmp[23]

	fatalities=int(tmp[24])
	fatalitiesAll[g]=fatalities

	if fatalities>=88:
		fatalitiesOver100+=1
	if fatalitiesAll[g]>10000:
		fatalitiesAll[g]=3000
	if fatalitiesAll[g]>3000:
		fatalitiesAll[g]=2500

plt.clf()
plt.plot(fatalitiesAll,'*')
plt.title('Fatalities in All African Conflicts since 1997')
plt.savefig(wdfigs+'fatalitiesAll.pdf')

plt.clf()
m = Basemap(llcrnrlon=-20,llcrnrlat=-38,urcrnrlon=54,urcrnrlat=41, projection='lcc',lon_0=23.5,lat_0=3.5)
m.drawcountries()
m.drawcoastlines()
#m.fillcontinents(color='white',lake_color='black',zorder=0)

latAllForPlot=np.zeros(shape=(fatalitiesOver100+1))
lonAllForPlot=np.zeros(shape=(fatalitiesOver100+1))
fatalitiesAllScaled=np.zeros(shape=(fatalitiesAll.shape))
i=-1
for k in range(len(latAll)):
	fatalitiesAllScaled[k]=((fatalitiesAll[k]-np.amin(fatalitiesAll))/(np.amax(fatalitiesAll)-np.amin(fatalitiesAll)))*(300-1)+1
	#if fatalitiesAll[k]>=88:
	#	i+=1
	#	latAllForPlot[i]=latAll[k]
	#	lonAllForPlot[i]=lonAll[k]
x,y = m(lonAll, latAll)
m.scatter(x,y,s=fatalitiesAllScaled,c=fatalitiesAll,cmap=cm.gist_heat_r,alpha=0.7,norm=colors.LogNorm())
plt.title('African Fatalities from Conflicts, 1997-Today')
plt.colorbar()

plt.savefig(wdfigs+'africa_fatalities',dpi=700)











