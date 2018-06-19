import csv
from math import sqrt
from sys import exit
import pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from mpl_toolkits.basemap import Basemap
import os
from scipy.stats import norm
import matplotlib as mpl
from matplotlib.patches import Polygon
import random
import matplotlib.cm as cm
import matplotlib.colors as colors

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
    xOut=math.sqrt(xOut)
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
plt.clf()

OpenStreetMap=False
Rivers=False
Flights=True

################## Ethiopia OpenStreetMap ##################
if OpenStreetMap:
	# Ethiopia lat lons
	map = Basemap(llcrnrlon=33,llcrnrlat=2.5,urcrnrlon=48.5,urcrnrlat=15.5, projection='lcc',lon_0=40.49,lat_0=9.145)
	
	map.readshapefile(wddata+'ethiopia-openstreetmap/gis_osm_roads_free_1', name='roads', drawbounds=True)
	
	ax = plt.gca()
	plt.title('Ethiopian Roads from Open Street Map')
	plt.yticks([])
	plt.xticks([])
	
	plt.savefig(wdfigs+'openstreetmap_ethiopia',dpi=700)

################## Africa Rivers ##################

if Rivers:
	plt.clf()
	map = Basemap(llcrnrlon=-20,llcrnrlat=-38,urcrnrlon=54,urcrnrlat=41, projection='lcc',lon_0=23.5,lat_0=3.5)
	
	map.drawcountries()
	map.drawcoastlines()
	
	map.readshapefile(wddata+'rivers_africa/waterways', name='rivers', drawbounds=True, color='b',linewidth=.7) # rn most accurate
	#map.readshapefile(wddata+'rivers_africa/ne_10m_rivers_lake_centerlines', name='rivers', drawbounds=True, color='b')
	#map.readshapefile(wddata+'rivers_africa/hydrobasins_africa', name='rivers', drawbounds=True, color='b')
	
	ax = plt.gca()
	plt.title('Africa Rivers')
	plt.yticks([])
	plt.xticks([])
	
	plt.savefig(wdfigs+'africa_rivers.pdf')
	#plt.savefig(wdfigs+'africa_rivers_basins',dpi=700)

################## Africa Flights ##################

if Flights:

	f=open(wddata+'Global_Flight_Data_Monthly/country_codes.csv')
	countryCodetoName={}
	for line in f:
		line=line[:-1]
		tmp=line.split(',')
		country=tmp[0]
		countryCode=tmp[1]
		countryCodetoName[countryCode]=country

	############### Get airport lon lats ###############
	f=open(wddata+'Global_Flight_Data_Monthly/Airports_2010.csv')
	
	airportCode=[]
	city=[]
	stateCode=[]
	countryCode=[]
	country=[]
	lon=np.zeros(shape=(3632))
	lat=np.zeros(shape=(3632))
	airportToLonLat={}
	airportCodetoNum={}
	i=-2
	for line in f:
		i+=1
		tmp=line.split(',')
		tmp[5]=tmp[5][:-2]
		
		if i==-1:
			continue

		airportCode.append(tmp[0])
		airportCodetoNum[airportCode[i]]=i
		city.append(tmp[1])
		stateCode.append(tmp[2])
		countryCode.append(tmp[3])
		country.append(countryCodetoName[countryCode[i]])
		lon[i]=float(tmp[4])
		lat[i]=float(tmp[5])
		airportToLonLat[airportCode[i]]=str(lon[i])+','+str(lat[i])

	############### Monthly Airport Use ###############

	f=open(wddata+'Global_Flight_Data_Monthly/Prediction_Monthly.csv')
	
	origin=[]
	dest=[]
	month=np.zeros(shape=(383048),dtype=int)
	numFlights=np.zeros(shape=(383048),dtype=int)
	errors=np.zeros(shape=(383048,2),dtype=int) # lower then upper
	i=-2
	for line in f:
		i+=1
		tmp=line.split(',')

		if i==-1:
			continue

		origin.append(tmp[0])
		dest.append(tmp[1])
		month[i]=int(tmp[2])-1
		numFlights[i]=int(tmp[3])
		errors[i,0]=int(tmp[4])
		errors[i,1]=int(tmp[5])

	flightsByAirport=np.zeros(shape=(len(airportCodetoNum)))
	i=-1
	test=0
	for code in origin:
		i+=1
		flightsByAirport[airportCodetoNum[code]]+=numFlights[i]
			
	plt.clf()
	#m = Basemap(llcrnrlon=-20,llcrnrlat=-38,urcrnrlon=54,urcrnrlat=41, projection='lcc',lon_0=23.5,lat_0=3.5)
	m = Basemap(projection='mill',llcrnrlat=-75,urcrnrlat=82,llcrnrlon=-180,urcrnrlon=180,lon_0=0)
	m.drawcountries()
	m.drawcoastlines()

	flightsScaled=np.zeros(shape=(flightsByAirport.shape))
	for k in range(len(lat)):
		flightsScaled[k]=((flightsByAirport[k]-np.amin(flightsByAirport))/(np.amax(flightsByAirport)-np.amin(flightsByAirport)))*(20-.05)+.05

	x,y = m(lon, lat)
	m.scatter(x,y,s=flightsScaled,c=flightsByAirport,cmap=cm.Blues,alpha=1,norm=colors.LogNorm())
	
	plt.title('World Airports, Scaled by 2010 Flights')
	plt.colorbar()

	plt.savefig(wdfigs+'world_flights',dpi=700)



















