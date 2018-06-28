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

wddata='/Users/lilllianpetersen/iiasa/data/'
wdfigs='/Users/lilllianpetersen/iiasa/figs/'
wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'

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

################### Read in Ethipia Shape Files ########################
r = shapefile.Reader(wddata+'night_light_ethiopia/2017/2013/Ethiopia_regions_2014/Ethiopia_regions_2014.shp')
rCountry = shapefile.Reader(wddata+'night_light_ethiopia/2017/2013/Ethiopia_regions_2014/ETH_outline.shp')
# Read in Regions
for i in range(len(r.shapes())):
	vars()['shapes'+str(i)] = r.shapes()[i]
	vars()['shape_points'+str(i)]=vars()['shapes'+str(i)].points
	vars()['polygon'+str(i)]=patches.Polygon(vars()['shape_points'+str(i)])
	vars()['points'+str(i)]=np.array(vars()['shape_points'+str(i)])

	plt.plot(vars()['points'+str(i)][:,0],vars()['points'+str(i)][:,1],'.')
	plt.savefig(wdfigs+'regions',dpi=600)

# Entire Country
shapeC=rCountry.shapes()[0]
shape_pointsC=shapeC.points
polygonC=patches.Polygon(shape_pointsC)
pointsC=np.array(shape_pointsC)
########################################################################

########################## Plot Fatalities #############################
#m = Basemap(llcrnrlon=-20,llcrnrlat=-38,urcrnrlon=54,urcrnrlat=41, projection='lcc',lon_0=23.5,lat_0=3.5)
plt.clf()
m = Basemap(llcrnrlon=33,llcrnrlat=2.5,urcrnrlon=48.5,urcrnrlat=15.5, projection='lcc',lon_0=40.49,lat_0=9.145,resolution='i')
m.drawcountries()
m.drawcoastlines(linewidth=2)
#m.fillcontinents(color='white',lake_color='black',zorder=0)

latAllForPlot=np.zeros(shape=(fatalitiesOver100+1))
lonAllForPlot=np.zeros(shape=(fatalitiesOver100+1))
fatalitiesAllScaled=np.zeros(shape=(fatalitiesAll.shape))
i=-1
for k in range(len(latAll)):
	fatalitiesAllScaled[k]=((fatalitiesAll[k]-np.amin(fatalitiesAll))/(np.amax(fatalitiesAll)-np.amin(fatalitiesAll)))*(300-1)+1

x,y = m(lonAll, latAll)
m.scatter(x,y,s=fatalitiesAllScaled,c=fatalitiesAll,cmap=cm.gist_heat_r,alpha=0.7,norm=colors.LogNorm())
plt.colorbar()
for i in range(len(r.shapes())):
	x,y=m(vars()['points'+str(i)][:,0],vars()['points'+str(i)][:,1])
	m.plot(x,y,linestyle='-')
plt.title('African Fatalities from Conflicts, 1997-Today')

#plt.savefig(wdfigs+'africa_fatalities',dpi=700)
plt.savefig(wdfigs+'ethiopia_fatalities',dpi=700)
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






