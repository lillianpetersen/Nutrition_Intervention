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
import shapefile
from shapely.geometry import shape, Point
import matplotlib.patches as patches
from math import sqrt

def checkIfInPolygon(lon, lat,polygons):
	i=-1
	done=False
	for polygon in polygons:
		i+=1
		point = Point(lon, lat)
		if polygon.contains(point)[0]==True:
			done=True
			return i
		else:
			continue
	if done==False:
		return -1

wddata='/Users/lilllianpetersen/iiasa/data/'
wdfigs='/Users/lilllianpetersen/iiasa/figs/'
wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'
plt.clf()

OpenStreetMap=False
Rivers=False
Flights=True

################## Ethiopian Regions ##################
plt.clf()
r = shapefile.Reader(wddata+'night_light_ethiopia/2017/2013/Ethiopia_regions_2014/Ethiopia_regions_2014.shp')
rCountry = shapefile.Reader(wddata+'night_light_ethiopia/2017/2013/Ethiopia_regions_2014/ETH_outline.shp')

# Read in Regions
polygons=[]
for i in range(len(r.shapes())):
	vars()['shapes'+str(i)] = r.shapes()[i]
	vars()['shape_points'+str(i)]=vars()['shapes'+str(i)].points
	vars()['polygon'+str(i)]=patches.Polygon(vars()['shape_points'+str(i)])
	polygons.append(vars()['polygon'+str(i)])
	vars()['points'+str(i)]=np.array(vars()['shape_points'+str(i)])

	plt.plot(vars()['points'+str(i)][:,0],vars()['points'+str(i)][:,1],'.')
	plt.savefig(wdfigs+'regions'+str(i),dpi=600)

# Entire Country
shapeC=rCountry.shapes()[0]
shape_pointsC=shapeC.points
polygonC=patches.Polygon(shape_pointsC)
pointsC=np.array(shape_pointsC)

##### Read in Ethiopia Regions Land Masses #######
f = open(wddata+'ethiopia_land_area.csv','r')
regions=[]
area=np.zeros(shape=(11))
perUrban=np.zeros(shape=(11))
i=-2
for line in f:
	i+=1
	if i==-1:
		continue
	tmp=line.split(',')
	regions.append(tmp[0])
	area[i]=float(tmp[2])
	perUrban[i]=float(tmp[6])
np.save(wdvars+'regions',regions)
np.save(wdvars+'percentUrban',perUrban)

################## Ethiopia OpenStreetMap ##################
if OpenStreetMap:
	

	####### Calculate km of roads in each region #######
	rRoads = shapefile.Reader(wddata+'ethiopia-openstreetmap/gis_osm_roads_free_1.shp')
	shapes=rRoads.shapes()

	distanceD=np.load(wdvars+'roadDistanceByRegionDegrees')
	#distanceD=np.zeros(shape=(len(r.shapes())))
	#s=-1
	#for shape in rRoads.shapes():
	#	s+=1
	#	print s
	#	shape_points=shape.points
	#	points=np.array(shape_points)
	#	for i in range(len(points)):
	#		if i==0:
	#			region=checkIfInPolygon(points[i,0],points[i,1],polygons)
	#			if region==-1:
	#				break
	#			continue
	#		if checkIfInPolygon(points[i,0],points[i,1],polygons)==region:
	#			dis=sqrt((points[i,0]-points[i-1,0])**2+(points[i,1]-points[i-1,1])**2)
	#			distanceD[region]+=dis
	#		else:
	#			region=checkIfInPolygon(points[i,0],points[i,1],polygons)
	#			if region==-1:
	#				break
	#			continue
	#np.save(wdvars+'roadDistanceByRegionDegrees',distanceD)
	kilometers=distanceD*111
	roadsByArea=kilometers/area
	np.save(wdvars+'roadsByArea',roadsByArea)
	exit()

	# Ethiopia lat lons
	plt.clf()
	map = Basemap(llcrnrlon=33,llcrnrlat=2.5,urcrnrlon=48.5,urcrnrlat=15.5, projection='lcc',lon_0=40.49,lat_0=9.145,resolution='i')
	map.drawcoastlines(linewidth=1.5)
	map.drawcountries()
	
	map.readshapefile(wddata+'ethiopia-openstreetmap/gis_osm_roads_free_1', name='roads', drawbounds=True)
	
	ax = plt.gca()
	plt.title('Ethiopian Roads from Open Street Map')
	
	plt.savefig(wdfigs+'openstreetmap_ethiopia',dpi=700)
	exit()

################## Africa Rivers ##################

if Rivers:
	plt.clf()
	#map = Basemap(llcrnrlon=-20,llcrnrlat=-38,urcrnrlon=54,urcrnrlat=41, projection='lcc',lon_0=23.5,lat_0=3.5)
	map = Basemap(llcrnrlon=33,llcrnrlat=2.5,urcrnrlon=48.5,urcrnrlat=15.5, projection='lcc',lon_0=40.49,lat_0=9.145)
	
	map.drawcountries()
	map.drawcoastlines()
	
	map.readshapefile(wddata+'rivers_africa/waterways', name='rivers', drawbounds=True, color='b',linewidth=.7) # rn most accurate
	#map.readshapefile(wddata+'rivers_africa/ne_10m_rivers_lake_centerlines', name='rivers', drawbounds=True, color='b')
	#map.readshapefile(wddata+'rivers_africa/hydrobasins_africa', name='rivers', drawbounds=True, color='b')
	
	ax = plt.gca()
	plt.title('Ethiopian Rivers')
	plt.yticks([])
	plt.xticks([])
	
	#plt.savefig(wdfigs+'africa_rivers.pdf')
	plt.savefig(wdfigs+'ethiopian_rivers',dpi=700)

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
	
	eth=[]
	airportCode=[]
	city=[]
	stateCode=[]
	countryCode=[]
	country=[]
	#lon=np.zeros(shape=(3632))
	#lat=np.zeros(shape=(3632))
	lon=np.zeros(shape=(17))
	lat=np.zeros(shape=(17))
	airportToLonLat={}
	airportToCountry={}
	airportCodetoNum={}
	k=-2
	i=-1
	for line in f:
		k+=1
		tmp=line.split(',')
		tmp[5]=tmp[5][:-2]
		if k==-1:
			continue

		countryCodetmp=tmp[3]
		countrytmp=countryCodetoName[countryCodetmp]
		if countrytmp!='Ethiopia':
			continue
		i+=1

		airportCode.append(tmp[0])
		airportCodetoNum[airportCode[i]]=i
		city.append(tmp[1])
		stateCode.append(tmp[2])
		countryCode.append(tmp[3])
		country.append(countryCodetoName[countryCode[i]])
		lon[i]=float(tmp[4])
		lat[i]=float(tmp[5])
		airportToLonLat[airportCode[i]]=str(lon[i])+','+str(lat[i])
		airportToCountry[airportCode[i]]=country[i]
		eth.append(airportCode[i])
	############### Monthly Airport Use ##############

	f=open(wddata+'Global_Flight_Data_Monthly/Prediction_Monthly.csv')
	
	origin=[]
	dest=[]
	month=np.zeros(shape=(383048),dtype=int)
	numFlights=np.zeros(shape=(383048),dtype=int)
	errors=np.zeros(shape=(383048,2),dtype=int) # lower then upper
	k=-2
	i=-1
	for line in f:
		k+=1
		tmp=line.split(',')

		if k==-1:
			continue

		try:
			airportToCountry[tmp[0]]
		except:
			continue
		i+=1
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
	m = Basemap(llcrnrlon=33,llcrnrlat=2.5,urcrnrlon=48.5,urcrnrlat=15.5, projection='lcc',lon_0=40.49,lat_0=9.145,resolution='i')
	#m = Basemap(projection='mill',llcrnrlat=-75,urcrnrlat=82,llcrnrlon=-180,urcrnrlon=180,lon_0=0)
	m.drawcountries()
	m.drawcoastlines()

	flightsScaled=np.zeros(shape=(flightsByAirport.shape))
	for k in range(len(flightsByAirport)):
		flightsScaled[k]=((flightsByAirport[k]-0)/(np.amax(flightsByAirport)-0))*(35-5)+5

	x,y = m(lon, lat)
	m.scatter(x,y,s=flightsScaled,c=flightsByAirport,cmap=cm.Blues,alpha=1,norm=colors.LogNorm())
	
	plt.title('World Airports, Scaled by 2010 Flights')
	plt.colorbar()

	plt.savefig(wdfigs+'world_flights',dpi=700)

	
	airportRegions=np.zeros(shape=(len(flightsByAirport)))
	flightsByRegion=np.zeros(shape=(len(r.shapes())))
	airportsByRegion=np.zeros(shape=(len(r.shapes())))
	for k in range(len(flightsByAirport)):
		region=checkIfInPolygon(lon[k],lat[k],polygons)
		if region==-1:
			break
		airportRegions[k]=region
		flightsByRegion[region]+=flightsByAirport[k]
		airportsByRegion[region]+=1

	flightsByArea=flightsByRegion/area
	airportsByArea=airportsByRegion/area

	np.save(wdvars+'flightsByRegion',flightsByRegion)
	np.save(wdvars+'flightsByArea',flightsByArea)
	np.save(wdvars+'airportsByRegion',airportsByRegion)
	np.save(wdvars+'airportsByArea',airportsByArea)


















