from __future__ import division
import csv
from math import sqrt
from sys import exit
import pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
#from mpl_toolkits.basemap import Basemap
import os
from scipy.stats import norm
import matplotlib as mpl
from matplotlib.patches import Polygon
import random
import matplotlib.cm as cm
import matplotlib.colors as colors
from PIL import Image
from osgeo import gdal
#import ogr
from IPython import embed
import shapefile
#from shapely.geometry import shape, Point
import matplotlib.patches as patches
from math import sin, cos, sqrt, atan2, radians, pi, degrees
# from geopy.geocoders import Nominatim
# geolocator = Nominatim()
# import geopy.distance
from scipy import ndimage
from scipy.signal import convolve2d
from sklearn import linear_model
from scipy.interpolate import RectSphereBivariateSpline
# from libtiff import TIFF
# import googlemaps
import json

try:
	wddata='/Users/lilllianpetersen/iiasa/data/supply_chain/'
	wdfigs='/Users/lilllianpetersen/iiasa/figs/'
	wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'
	f=open(wddata+'population/CAPITALVERSIONcasenumbers.csv','r')
except:
	wddata='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/supply_chain/data/'
	wdfigs='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/supply_chain/figs/'
	wdvars='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/supply_chain/vars/'

##capitalonly
subsaharancountry=[]
subsaharancapital=[]
indexedwasting=np.zeros(shape=43)
indexedSAM=np.zeros(shape=43)
indexedstunting=np.zeros(shape=43)
indexedMAM=np.zeros(shape=43)
f=open(wddata+'population/CAPITALVERSIONcasenumbers.csv','r')
i=-1
for line in f:
    i+=1
    tmp=line.split(',')
    subsaharancountry.append(tmp[0])
    indexedwasting[i]=float(tmp[1])
    indexedSAM[i]=float(tmp[2])*150*11.66
    indexedMAM[i]=float(tmp[3])*50*(365/75)
    indexedstunting[i]=float(tmp[4])
    subsaharancapital.append(tmp[5][:-1])

np.savetxt(wddata + "optiarrays/SAMdemand.csv", indexedSAM, delimiter=",")
np.savetxt(wddata + "optiarrays/MAMdemand.csv", indexedMAM, delimiter=",")

# subsaharancountry=[]
# subsaharancapital=[]
# indexedwasting=np.zeros(shape=50)
# indexedSAM=np.zeros(shape=50)
# indexedstunting=np.zeros(shape=50)
# indexedMAM=np.zeros(shape=50)
# f=open(wddata+'population/nationalcasenumbers.csv','r')
# i=-1
# for line in f:
#     i+=1
#     tmp=line.split(',')
#     subsaharancountry.append(tmp[0])
#     indexedwasting[i]=float(tmp[1])
#     indexedSAM[i]=float(tmp[2])
#     indexedMAM[i]=float(tmp[3])
#     indexedstunting[i]=float(tmp[4])
#     for j in range(len(subsaharancapital1)):
#         if tmp[0]==subsaharancountry1[j]:
#             subsaharancapital.append(subsaharancapital1[j])
       
     	
### cost per country
countrycosted=[]
rutfprice=[]
rusfprice=[]
scplusprice=[]
f=open(wddata+'foodstuffs/pricesCorrected.csv')
code=np.zeros(shape=(247),dtype=int)
i=-1
for line in f:
	i+=1
	tmp=line.split(',')
	countrycosted.append(tmp[0])
	rutfprice.append(tmp[1])
	rusfprice.append(tmp[2])
	scplusprice.append(tmp[3][:-1])

countrycosted[6]="Congo (Republic of the)"
countrycosted[7]="Cote d'Ivoire"

indexedrutf=np.zeros(shape=43)
indexedrusf=np.zeros(shape=43)
indexedscp=np.zeros(shape=43)
conversionlist=[]
convertarray=np.zeros(shape=43)
for i in range(len(countrycosted)):
    for j in range(len(subsaharancountry)):
        if countrycosted[i]==subsaharancountry[j]:
            conversionlist.append(j)
            convertarray[j]=i
            indexedrutf[j]=rutfprice[i]
            indexedrusf[j]=rusfprice[i]
            indexedscp[j]=scplusprice[i][:-1]


#read in transport cost
transportcostArray=np.zeros(shape=(43,43))
f=open(wddata+'travel_time/transportationcostsarray.csv','r')
i=-1
for line in f:
    i+=1
    tmp=line.split(',')
    for j in range(len(tmp)):
        transportcostArray[i,j]=tmp[j]

# import and export costs
importExportCosts=np.zeros(shape=(transportcostArray.shape))
for x in range(len(subsaharancountry)):
	exportCost=-9999
	fx=open(wddata+'trading_across_borders2017.csv','r')
	xCountry=subsaharancountry[x]
	for line in fx:
		tmp=line.split(',')
		if tmp[0]==xCountry:
			exportCost=float(tmp[4])+float(tmp[6])
			break
	print exportCost,xCountry

	for y in range(len(subsaharancountry)):
		importCost=-9999
		yCountry=subsaharancountry[y]
		if xCountry==yCountry:
			continue

		fy=open(wddata+'trading_across_borders2017.csv','r')
		for line in fy:
			tmp=line.split(',')
			if tmp[0]==yCountry:
				importCost=float(tmp[8])+float(tmp[10])
				break

		importExportCosts[x,y]=importCost+exportCost

#cost dabber RUTF ########################################################################
rutfcostarray=np.zeros(shape=(24,43))
for i in range(len(subsaharancountry)):
    if(indexedrutf[i]!=0):
        for j in range(len(indexedSAM)):
			# sums ingredient and transport cost, converts to $/100g delivered
            rutfcostarray[int(convertarray[i]),j]=indexedrutf[i]+100/1000000*transportcostArray[i,j]

        # if(indexedscp[i]*2<indexedrusf[i]):
        #     for j in range(len(indexedSAM)):
        #         costarray[int(convertarray[i]),j]=indexedrutf[i]*indexedSAM[j]*53*11.66+indexedSAM[j]*53*11.66*92/1000000*transportcostArray[i,j]+indexedscp[i]*indexedMAM[j]*50*(365/75)*2+indexedMAM[j]*50*(365/75)*200/1000000*transportcostArray[i,j]
        # else:
        #     for j in range(len(indexedSAM)):
        #         costarray[int(convertarray[i]),j]=indexedrutf[i]*indexedSAM[j]*53*11.66+indexedSAM[j]*53*11.66*92/1000000*transportcostArray[i,j]+indexedrusf[i]*indexedMAM[j]*50*(365/75)+indexedMAM[j]*50*(365/75)*100/1000000*transportcostArray[i,j]
    
np.savetxt(wddata + "foodstuffs/rutfcostarray.csv", rutfcostarray, delimiter=",")

# costarrayfixed=np.zeros(shape=(20,43))
# f=open(wddata+'foodstuffs/rutfcostarray.csv','r')
# i=-1
# for line in f:
#     i+=1
#     tmp=line.split(',')
#     for j in range(len(tmp)):
#         costarrayfixed[i,j]=tmp[j]

rutfdictionary={}
### array to dict
for i in range(len(countrycosted)):
    rutfdictionary[countrycosted[i]]={}
    for j in range(len(subsaharancountry)):
        rutfdictionary[countrycosted[i]][subsaharancountry[j]]=rutfcostarray[i,j]

with open(wddata+'optiarrays/rutfdictionary.json', 'w') as fp:
    json.dump(rutfdictionary, fp, sort_keys=True)

#cost dabber MAM ########################################################################
mamcostarray=np.zeros(shape=(24,43))
for i in range(len(subsaharancountry)):
    if(indexedrusf[i]!=0):
        for j in range(len(indexedSAM)):
                if(indexedscp[i]*2<indexedrusf[i]):
                    for j in range(len(indexedSAM)):
                        mamcostarray[int(convertarray[i]),j]=indexedscp[i]+200/1000000*transportcostArray[i,j]
                        print subsaharancountry[i] +' SC+'
                else:
                    for j in range(len(indexedSAM)):
                        mamcostarray[int(convertarray[i]),j]=indexedrusf[i]+100/1000000*transportcostArray[i,j]

np.savetxt(wddata + "foodstuffs/rutfcostarray.csv", rutfcostarray, delimiter=",")

mamdictionary={}
### array to dict
for i in range(len(countrycosted)):
    mamdictionary[countrycosted[i]]={}
    for j in range(len(subsaharancountry)):
        mamdictionary[countrycosted[i]][subsaharancountry[j]]=mamcostarray[i,j]

with open(wddata+'optiarrays/mamdictionary.json', 'w') as fp:
    json.dump(mamdictionary, fp, sort_keys=True)
