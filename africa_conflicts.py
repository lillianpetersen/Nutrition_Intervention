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
import shapefile
from shapely.geometry import shape, Point
import matplotlib.patches as patches

def checkIfInPolygon(lon, lat,polygon):
	point = Point(lon, lat)
	return polygon.contains(point)

def scale(var):
	varScaled=(var-np.amin(var))/(np.amax(var)-np.amin(var))
	return varScaled


wddata='/Users/lilllianpetersen/iiasa/data/'
wdfigs='/Users/lilllianpetersen/iiasa/figs/'
wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'

f=open(wddata+'Africa_conflictss.csv','r')

monthDict={'January':0,'February':1,'March':2,'April':3,'May':4,'June':5,'July':6,'August':7,'September':8,'October':9,'November':10,'December':11}

allDataIDs=np.zeros(shape=(161939),dtype=int)
timePrecision=-9999*np.zeros(shape=(119,12,2000))
eventTypeDict={}
eventTypeCounter=0
countriesAll=[]
countryCodeAll=[]
latAll=-9999*np.zeros(shape=(119,12,2000))
lonAll=-9999*np.ones(shape=(119,12,2000))
fatalitiesAll=-9999*np.ones(shape=(119,12,2000))
fatalitiesOver100=0

i=-2
g=-1
d=-1
lastMonth=0
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
	d+=1

	date=tmp[4]
	date=date.replace('"','')
	date=date.split(' ')
	day=int(date[0])
	m=monthDict[date[1]]
	year=int(date[2])-1900
	if m!=lastMonth:
		d=0
	lastMonth=m

	timePrecision[year,m,d]=tmp[6]
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
	latAll[year,m,d]=lat
	lonAll[year,m,d]=lon
	geoAccuracy=tmp[20] #Probably not important

	notes=tmp[23]

	fatalities=int(tmp[24])
	fatalitiesAll[year,m,d]=fatalities

	if fatalities>=88:
		fatalitiesOver100+=1

plt.clf()
plt.plot(np.ma.compressed(fatalitiesAll[fatalitiesAll!=0]),'*')
plt.title('Fatalities in All African Conflicts since 1997')
plt.savefig(wdfigs+'fatalitiesAll.pdf')

fatalitiesAll=np.array(fatalitiesAll)
latAll=np.array(latAll)
lonAll=np.array(lonAll)
timePrecision=np.array(timePrecision)

fatalMask=np.zeros(shape=(fatalitiesAll.shape),dtype=int)
fatalMask[fatalitiesAll==-9999]=1
fatalMask[fatalitiesAll==0]=1
fatalitiesAll=np.ma.masked_array(fatalitiesAll,fatalMask)
latAll=np.ma.masked_array(latAll,fatalMask)
lonAll=np.ma.masked_array(lonAll,fatalMask)
timePrecision=np.ma.masked_array(timePrecision,fatalMask)

################### Read in Ethipia Shape Files ########################
#plt.clf()
#r = shapefile.Reader(wddata+'night_light_ethiopia/2017/2013/Ethiopia_regions_2014/Ethiopia_regions_2014.shp')
#rCountry = shapefile.Reader(wddata+'night_light_ethiopia/2017/2013/Ethiopia_regions_2014/ETH_outline.shp')
## Read in Regions
#for i in range(len(r.shapes())):
#	vars()['shapes'+str(i)] = r.shapes()[i]
#	vars()['shape_points'+str(i)]=vars()['shapes'+str(i)].points
#	vars()['polygon'+str(i)]=patches.Polygon(vars()['shape_points'+str(i)])
#	vars()['points'+str(i)]=np.array(vars()['shape_points'+str(i)])
#
#	plt.plot(vars()['points'+str(i)][:,0],vars()['points'+str(i)][:,1],'.')
#	plt.savefig(wdfigs+'regions',dpi=600)
#
## Entire Country
#shapeC=rCountry.shapes()[0]
#shape_pointsC=shapeC.points
#polygonC=patches.Polygon(shape_pointsC)
#pointsC=np.array(shape_pointsC)
########################################################################


########################## Plot Fatalities #############################
plt.clf()
m = Basemap(llcrnrlon=-20,llcrnrlat=-38,urcrnrlon=54,urcrnrlat=41, projection='lcc',lon_0=23.5,lat_0=3.5)
m.drawcountries()
m.drawcoastlines(linewidth=2)
#m.fillcontinents(color='white',lake_color='blue',zorder=0)

fatalitiesScaled=scale(np.ma.compressed(fatalitiesAll[113:116]))*20
fatalitiesPlot=np.ma.compressed(fatalitiesAll[113:116])

x,y = m(np.ma.compressed(lonAll[113:116]), np.ma.compressed(latAll[113:116]))
m.scatter(x,y,s=fatalitiesScaled,c=fatalitiesPlot,cmap=cm.gist_heat_r,alpha=0.7,norm=colors.LogNorm())
plt.colorbar()
plt.title('African Fatalities from Conflicts, 2013-2015')

plt.savefig(wdfigs+'africa_fatalities',dpi=700)

conflictCoord=np.zeros(shape=(119,12,2000,2))
conflictCoord[:,:,:,0]=latAll
conflictCoord[:,:,:,1]=lonAll

fatalMask3=np.zeros(shape=(119,12,2000,2))
for i in range(2):
	fatalMask3[:,:,:,i]=fatalMask

conflictCoord=np.ma.masked_array(conflictCoord,fatalMask3)

fatalitiesAll.dump(wdvars+'africa_fatalities')
conflictCoord.dump(wdvars+'conflictCoord')


exit()
########################################################################

################ Divide Ethiopia Conflicts into Regions ################

lonLatMask=np.ones(shape=(len(r.shapes()),latAll.shape[0]),dtype=bool)
lonLatMaskSimple=np.ones(shape=(latAll.shape[0]),dtype=bool)
lonLatNums=-1*np.ones(shape=(latAll.shape[0]))
for i in range(len(latAll)):
	if countriesAll[i]=='Ethiopia':
		lonLatMaskSimple[i]=False
		for k in range(len(r.shapes())):
			if checkIfInPolygon(lonAll[i],latAll[i],vars()['polygon'+str(k)])[0]==True:
				lonLatMask[k,i]=False
				lonLatNums[i]=k
				break
		if lonLatNums[i]==-1:
			lonLatMask[:,i]=True
			lonLatMaskSimple[i]=True

#np.save(wdvars+'regionMaskforNightLights',regionMask)
lonAllM=np.ma.compressed(np.ma.masked_array(lonAll,lonLatMaskSimple))
latAllM=np.ma.compressed(np.ma.masked_array(latAll,lonLatMaskSimple))
lonLatNumsM=np.ma.compressed(np.ma.masked_array(lonLatNums,lonLatMaskSimple))
fatalitiesAllScaledM=np.ma.compressed(np.ma.masked_array(fatalitiesAllScaled,lonLatMaskSimple))
plt.clf()
plt.scatter(lonAllM,latAllM,s=fatalitiesAllScaled,c=lonLatNumsM,cmap=cm.nipy_spectral)
for i in range(len(r.shapes())):
	plt.plot(vars()['points'+str(i)][:,0],vars()['points'+str(i)][:,1],linestyle='-',c='k',linewidth=.5)
plt.title('Ethiopian Conflicts by Region')
plt.colorbar()
plt.axes().set_aspect('equal', 'datalim')
plt.savefig(wdfigs+'testMask',dpi=700)






