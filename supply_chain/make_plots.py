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

subsaharancountry = np.load(wdvars+'subsaharancountry.npy')
for i in range(len(subsaharancountry)):
	subsaharancountry[i]=subsaharancountry[i].replace(' ','_')

optiLevel = ['AllarOpti','LocalOpti']
loopvar = ['SHIPCOST','impexp','strtup','truckfactor']
Titles = ['Shipping','Import/Export','Startup','Trucking']
titles = ['shipping','importexport','startup','trucking']
factor = np.array([0.2,0.2,0.5,0.2])
maxs = np.array([2.01,4.01,9.51,4.01])

for L in range(len(optiLevel)):
	for V in range(len(loopvar)):
		File =optiLevel[L]+'_'+loopvar[V]
		print File
		cost=np.zeros(shape=(14))
		factoryNum=np.zeros(shape=(14))
		pctTrans=np.zeros(shape=(14))
		pctIngredient=np.zeros(shape=(14))
		avgShipments=np.zeros(shape=(14))
		sizeAvg=np.zeros(shape=(14))
		sizeStdDev=np.zeros(shape=(14))
		sizeMinMax=np.zeros(shape=(14,2))
		factorySizeAll=np.zeros(shape=(14,24))
		factorySizeAllMask=np.ones(shape=(14,24),dtype=bool)
		
		factoryPct=np.zeros(shape=(14,43))

		countrycosted = np.load(wdvars+'countrycosted.npy')
                
                for i in range(len(countrycosted)):
                    countrycosted[i]=countrycosted[i].replace(' ','_')

		i=-1
		for s in np.arange(factor[V],maxs[V],factor[V]):
			s=np.round(s,2)
			i+=1
			countriesWfactories=[]
			factorySize=[]
		
			f = open(wddata+'results/'+str(File)+'/'+str(File)+str(s)+'.csv')
			k=-1
			j=-1
			for line in f:
				k+=1
				tmp=line.split(',')
				if k==0:
					cost[i]=float(tmp[1])
				elif k==1:
					factoryNum[i]=float(tmp[1])
				elif k==2:
					pctTrans[i]=float(tmp[1])
				elif k==3:
					pctIngredient[i]=float(tmp[1])
				elif k==4:
					avgShipments[i]=float(tmp[1])
				else:
					j+=1
					countryi=tmp[0]

					k=np.where(countryi==countrycosted)[0][0]
					country=countrycosted[i]

					#country=country.replace('_',' ')
					countriesWfactories.append(country)
					size=float(tmp[1])
					factorySize.append(size)
					factorySizeAll[i,j]=size
					factorySizeAllMask[i,j]=False
					
					c=np.where(country==subsaharancountry)[0][0]
					factoryPct[i,c]=size
		
			factorySize=np.array(factorySize)
			sizeAvg[i]=np.mean(factorySize)
			sizeStdDev[i]=np.std(factorySize)
			sizeMinMax[i,0]=np.amin(factorySize)
			sizeMinMax[i,1]=np.amax(factorySize)
		factorySizeAll=np.ma.masked_array(factorySizeAll,factorySizeAllMask)
		
		totalCapacity = np.sum(factorySizeAll,axis=1)[0]
		factoryPct = np.swapaxes(factoryPct,0,1)
		factoryPct = 100*factoryPct/totalCapacity
		
		x = np.array([0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5])
		x = x*100
		fig = plt.figure(figsize=(9, 6))
		#### cost ####
		plt.clf()
		plt.plot(x,cost,'*-')
		plt.title('Total Cost by Trucking Cost')
		plt.xlabel('Trucking Cost, % of Today')
		plt.ylabel('Total Procurement Cost for One Year')
		plt.grid(True)
		plt.savefig(wdfigs+'totalCost_vs_truckingCost.pdf')
		
		#### factoryNum ####
		plt.clf()
		plt.plot(x,factoryNum,'*-')
		plt.title('Number of Factories by Trucking Cost')
		plt.xlabel('Trucking Cost, % of Today')
		plt.ylabel('Number of Factories')
		plt.grid(True)
		plt.savefig(wdfigs+'factoryNum_vs_truckingCost.pdf')
		
		#### factorySize ####
		plt.clf()
		plt.plot(x,sizeAvg,'b*-')
		plt.plot(x,sizeMinMax[:,0],'b*--')
		plt.plot(x,sizeMinMax[:,1],'b*--')
		plt.title('Avg Factory Size by Trucking Cost')
		plt.xlabel('Trucking Cost, % of Today')
		plt.ylabel('Avg Factory Size')
		plt.grid(True)
		plt.savefig(wdfigs+'factorySize_vs_truckingCost.pdf')
		
		#### factorySizeAll ####
		plt.clf()
		plt.plot(x,factorySizeAll,'b*')
		plt.title('Factory Size by Trucking Cost')
		plt.xlabel('Trucking Cost, % of Today')
		plt.ylabel('Factory Size')
		plt.grid(True)
		plt.savefig(wdfigs+'factorySize_vs_truckingCost.pdf')
		
		#### PctTrans ####
		plt.clf()
		plt.plot(x,pctTrans*100,'b*-')
		plt.title('% Cost that is Transport by Trucking Cost')
		plt.xlabel('Trucking Cost, % of Today')
		plt.ylabel('% of Total Cost that is Transport')
		plt.grid(True)
		plt.savefig(wdfigs+'pctTrans_vs_truckingCost.pdf')
		
		#### PctIngredient ####
		plt.clf()
		plt.plot(x,pctIngredient*100,'b*-')
		plt.title('% Cost that is Ingredient by Trucking Cost')
		plt.xlabel('Trucking Cost, % of Today')
		plt.ylabel('% of Total Cost that is Ingredient')
		plt.grid(True)
		plt.savefig(wdfigs+'pctIngredient_vs_truckingCost.pdf')
		
		#### FactoryPct ####
		factoryCountries = []
		plt.clf()
		fig = plt.figure(figsize=(15, 10))
		
		for t in range(len(subsaharancountry)):
			country = subsaharancountry[t]
			c=np.where(country==subsaharancountry)[0][0]
			if country=='Congo (Republic of the)':
				country='RepOfCongo'
			elif country=='Congo':
				country='DRC'
			country=country.replace(' ','_')
		
			if np.amax(factoryPct[c])!=0:
				vars()[country+'Pct'] = factoryPct[c]
				factoryCountries.append(country)
		
		countryComparison = np.zeros(shape=(len(factoryCountries)))
		for q in range(len(factoryCountries)):
			country = factoryCountries[q]
			countryComparison[q] = vars()[country+'Pct'][-1]
		sIndex=np.argsort(countryComparison)
		factoryCountries=np.array(factoryCountries)
		factoryCountries=factoryCountries[sIndex][::-1]
		
		
		colors = ['firebrick','navy','m','lawngreen','darkorange','darkgreen','crimson','yellow','mediumseagreen','darkslategray','indianred','goldenrod','blueviolet','lightseagreen','powderblue','mediumpurple','navajowhite','peru','tomato','midnightblue','magenta','forestgreen']
		width = 5 
		pvars=[]
		for i in range(len(factoryCountries)):
			country = factoryCountries[i]
		
			if i==0:
				vars()['p'+str(i)] = plt.bar(x, vars()[country+'Pct'], width, color=colors[i], )
				bottomStuff = vars()[country+'Pct']
			else:
				vars()['p'+str(i)] = plt.bar(x, vars()[country+'Pct'], width, color=colors[i], bottom = bottomStuff)
				bottomStuff+=vars()[country+'Pct']
		
			print vars()[country+'Pct']
			pvars.append(vars()['p'+str(i)])
			
		
		plt.title('Percent Each Facotry Produces by Trucking Cost')
		plt.xlabel('Trucking Cost, % of Today')
		plt.ylabel('% of Total Production')
		plt.legend((pvars[::-1]),(factoryCountries[::-1]))
		plt.savefig(wdfigs+'factoryPct_vs_truckingCost.pdf')
# 
# 
# 
# ##############################################
# # Import/Export Costs
# ##############################################
# cost=np.zeros(shape=(14))
# factoryNum=np.zeros(shape=(14))
# pctTrans=np.zeros(shape=(14))
# pctIngredient=np.zeros(shape=(14))
# sizeAvg=np.zeros(shape=(14))
# sizeStdDev=np.zeros(shape=(14))
# sizeMinMax=np.zeros(shape=(14,2))
# factorySizeAll=np.zeros(shape=(14,24))
# factorySizeAllMask=np.ones(shape=(14,24),dtype=bool)
# 
# factoryPct=np.zeros(shape=(14,43))
# #factoryPctMask=np.ones(shape=(14,43),dtype=bool)
# 
# i=-1
# for s in np.arange(0.2,1.6,0.1):
# 	s=np.round(s,2)
# 	i+=1
# 	countriesWfactories=[]
# 	factorySize=[]
# 
# 	f = open(wddata+'borderfactor_costs/bordercostfactor'+str(s)+'.csv')
# 	k=-1
# 	j=-1
# 	for line in f:
# 		k+=1
# 		tmp=line.split(',')
# 		if k==0:
# 			cost[i]=float(tmp[1])
# 		elif k==1:
# 			factoryNum[i]=float(tmp[1])
# 		elif k==2:
# 			pctTrans[i]=float(tmp[1])
# 		elif k==3:
# 			pctIngredient[i]=float(tmp[1])
# 		else:
# 			j+=1
# 			country=tmp[0]
# 			country=country.replace('_',' ')
# 			countriesWfactories.append(country)
# 			size=float(tmp[1])
# 			factorySize.append(size)
# 			factorySizeAll[i,j]=size
# 			factorySizeAllMask[i,j]=False
# 			
# 			c=np.where(country==subsaharancountry)[0][0]
# 			factoryPct[i,c]=size
# 
# 	factorySize=np.array(factorySize)
# 	sizeAvg[i]=np.mean(factorySize)
# 	sizeStdDev[i]=np.std(factorySize)
# 	sizeMinMax[i,0]=np.amin(factorySize)
# 	sizeMinMax[i,1]=np.amax(factorySize)
# factorySizeAll=np.ma.masked_array(factorySizeAll,factorySizeAllMask)
# 
# totalCapacity = np.sum(factorySizeAll,axis=1)[0]
# factoryPct = np.swapaxes(factoryPct,0,1)
# factoryPct = 100*factoryPct/totalCapacity
# 
# x = np.array([0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5])
# x = x*100
# fig = plt.figure(figsize=(9, 6))
# #### cost ####
# plt.clf()
# plt.plot(x,cost,'*-')
# plt.title('Total Cost by Import/Export Cost')
# plt.xlabel('Import/Export Cost, % of Today')
# plt.ylabel('Total Procurement Cost for One Year')
# plt.grid(True)
# plt.savefig(wdfigs+'totalCost_vs_borderCost.pdf')
# 
# #### factoryNum ####
# plt.clf()
# plt.plot(x,factoryNum,'*-')
# plt.title('Number of Factories by Import/Export Cost')
# plt.xlabel('Import/Export Cost, % of Today')
# plt.ylabel('Number of Factories')
# plt.grid(True)
# plt.savefig(wdfigs+'factoryNum_vs_borderCost.pdf')
# 
# #### factorySize ####
# plt.clf()
# plt.plot(x,sizeAvg,'b*-')
# plt.plot(x,sizeMinMax[:,0],'b*--')
# plt.plot(x,sizeMinMax[:,1],'b*--')
# plt.title('Avg Factory Size by Import/Export Cost')
# plt.xlabel('Import/Export Cost, % of Today')
# plt.ylabel('Avg Factory Size')
# plt.grid(True)
# plt.savefig(wdfigs+'factorySize_vs_borderCost.pdf')
# 
# #### factorySizeAll ####
# plt.clf()
# plt.plot(x,factorySizeAll,'b*')
# plt.title('Factory Size by Import/Export Cost')
# plt.xlabel('Import/Export Cost, % of Today')
# plt.ylabel('Factory Size')
# plt.grid(True)
# plt.savefig(wdfigs+'factorySize_vs_borderCost.pdf')
# 
# #### PctTrans ####
# plt.clf()
# plt.plot(x,pctTrans*100,'b*-')
# plt.title('% Cost that is Transport by Import/Export Cost')
# plt.xlabel('Import/Export Cost, % of Today')
# plt.ylabel('% of Total Cost that is Transport')
# plt.grid(True)
# plt.savefig(wdfigs+'pctTrans_vs_borderCost.pdf')
# 
# #### PctIngredient ####
# plt.clf()
# plt.plot(x,pctIngredient*100,'b*-')
# plt.title('% Cost that is Ingredient by Import/Export Cost')
# plt.xlabel('Import/Export Cost, % of Today')
# plt.ylabel('% of Total Cost that is Ingredient')
# plt.grid(True)
# plt.savefig(wdfigs+'pctIngredient_vs_borderCost.pdf')
# 
# #### FactoryPct ####
# factoryCountries = []
# plt.clf()
# fig = plt.figure(figsize=(15, 10))
# 
# N = 14
# for i in range(len(subsaharancountry)):
# 	country = subsaharancountry[i]
# 	c=np.where(country==subsaharancountry)[0][0]
# 	if country=='Congo (Republic of the)':
# 		country='RepOfCongo'
# 	elif country=='Congo':
# 		country='DRC'
# 	country=country.replace(' ','')
# 
# 	if np.amax(factoryPct[c])!=0:
# 		vars()[country+'Pct'] = factoryPct[c]
# 		factoryCountries.append(country)
# 
# countryComparison = np.zeros(shape=(len(factoryCountries)))
# for i in range(len(factoryCountries)):
# 	country = factoryCountries[i]
# 	countryComparison[i] = vars()[country+'Pct'][-1]
# sIndex=np.argsort(countryComparison)
# factoryCountries=np.array(factoryCountries)
# factoryCountries=factoryCountries[sIndex][::-1]
# 
# 
# colors = ['firebrick','navy','m','lawngreen','darkorange','darkgreen','crimson','yellow','mediumseagreen','darkslategray','indianred','goldenrod','blueviolet','lightseagreen','powderblue','mediumpurple','navajowhite','peru','tomato','midnightblue','magenta','forestgreen']
# width = 5 
# for i in range(len(factoryCountries)):
# 	country = factoryCountries[i]
# 
# 	if i==0:
# 		vars()['p'+str(i)] = plt.bar(x, vars()[country+'Pct'], width, color=colors[i], )
# 		bottomStuff = vars()[country+'Pct']
# 	else:
# 		vars()['p'+str(i)] = plt.bar(x, vars()[country+'Pct'], width, color=colors[i], bottom = bottomStuff)
# 		bottomStuff+=vars()[country+'Pct']
# 
# 	print vars()[country+'Pct']
# 
# plt.title('Percent Each Facotry Produces by Import/Export Cost')
# plt.xlabel('Import/Export Cost, % of Today')
# plt.ylabel('% of Total Production')
# plt.legend(([p0[0],p1[0],p2[0],p3[0],p4[0],p5[0],p6[0],p7[0],p8[0],p9[0],p10[0],p11[0],p12[0],p13[0],p14[0],p15[0],p16[0],p17[0],p18[0],p19[0],p20[0],p21[0]][::-1]),(factoryCountries[::-1]))
# plt.savefig(wdfigs+'factoryPct_vs_borderCost.pdf')
# 
# ##############################################
# # Tariffs 
# ##############################################
# cost=np.zeros(shape=(11))
# factoryNum=np.zeros(shape=(11))
# pctTrans=np.zeros(shape=(11))
# pctIngredient=np.zeros(shape=(11))
# sizeAvg=np.zeros(shape=(11))
# sizeStdDev=np.zeros(shape=(11))
# sizeMinMax=np.zeros(shape=(11,2))
# factorySizeAll=np.zeros(shape=(11,24))
# factorySizeAllMask=np.ones(shape=(11,24),dtype=bool)
# 
# factoryPct=np.zeros(shape=(11,43))
# #factoryPctMask=np.ones(shape=(14,43),dtype=bool)
# 
# i=-1
# for s in np.arange(0.0,0.31,0.03):
# 	s=np.round(s,2)
# 	i+=1
# 	countriesWfactories=[]
# 	factorySize=[]
# 
# 	f = open(wddata+'tariff_multiplier/tarriflevel'+str(s)+'.csv')
# 	k=-1
# 	j=-1
# 	for line in f:
# 		k+=1
# 		tmp=line.split(',')
# 		if k==0:
# 			cost[i]=float(tmp[1])
# 		elif k==1:
# 			factoryNum[i]=float(tmp[1])
# 		elif k==2:
# 			pctTrans[i]=float(tmp[1])
# 		elif k==3:
# 			pctIngredient[i]=float(tmp[1])
# 		else:
# 			j+=1
# 			country=tmp[0]
# 			country=country.replace('_',' ')
# 			countriesWfactories.append(country)
# 			size=float(tmp[1])
# 			factorySize.append(size)
# 			factorySizeAll[i,j]=size
# 			factorySizeAllMask[i,j]=False
# 			
# 			c=np.where(country==subsaharancountry)[0][0]
# 			factoryPct[i,c]=size
# 
# 	factorySize=np.array(factorySize)
# 	sizeAvg[i]=np.mean(factorySize)
# 	sizeStdDev[i]=np.std(factorySize)
# 	sizeMinMax[i,0]=np.amin(factorySize)
# 	sizeMinMax[i,1]=np.amax(factorySize)
# factorySizeAll=np.ma.masked_array(factorySizeAll,factorySizeAllMask)
# 
# totalCapacity = np.sum(factorySizeAll,axis=1)[0]
# factoryPct = np.swapaxes(factoryPct,0,1)
# factoryPct = 100*factoryPct/totalCapacity
# 
# x = np.array([0.0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3])
# x = x*100
# fig = plt.figure(figsize=(9, 6))
# #### cost ####
# plt.clf()
# plt.plot(x,cost,'*-')
# plt.title('Total Cost by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('Total Procurement Cost for One Year')
# plt.grid(True)
# plt.savefig(wdfigs+'totalCost_vs_tariffCost.pdf')
# 
# #### factoryNum ####
# plt.clf()
# plt.plot(x,factoryNum,'*-')
# plt.title('Number of Factories by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('Number of Factories')
# plt.grid(True)
# plt.savefig(wdfigs+'factoryNum_vs_tariffCost.pdf')
# 
# #### factorySize ####
# plt.clf()
# plt.plot(x,sizeAvg,'b*-')
# plt.plot(x,sizeMinMax[:,0],'b*--')
# plt.plot(x,sizeMinMax[:,1],'b*--')
# plt.title('Avg Factory Size by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('Avg Factory Size')
# plt.grid(True)
# plt.savefig(wdfigs+'factorySize_vs_tariffCost.pdf')
# 
# #### factorySizeAll ####
# plt.clf()
# plt.plot(x,factorySizeAll,'b*')
# plt.title('Factory Size by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('Factory Size')
# plt.grid(True)
# plt.savefig(wdfigs+'factorySize_vs_tariffCost.pdf')
# 
# #### PctTrans ####
# plt.clf()
# plt.plot(x,pctTrans*100,'b*-')
# plt.title('% Cost that is Transport by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('% of Total Cost that is Transport')
# plt.grid(True)
# plt.savefig(wdfigs+'pctTrans_vs_tariffCost.pdf')
# 
# #### PctIngredient ####
# plt.clf()
# plt.plot(x,pctIngredient*100,'b*-')
# plt.title('% Cost that is Ingredient by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('% of Total Cost that is Ingredient')
# plt.grid(True)
# plt.savefig(wdfigs+'pctIngredient_vs_tariffCost.pdf')
# 
# #### FactoryPct ####
# factoryCountries = []
# plt.clf()
# fig = plt.figure(figsize=(15, 10))
# 
# for i in range(len(subsaharancountry)):
# 	country = subsaharancountry[i]
# 	c=np.where(country==subsaharancountry)[0][0]
# 	if country=='Congo (Republic of the)':
# 		country='RepOfCongo'
# 	elif country=='Congo':
# 		country='DRC'
# 	country=country.replace(' ','')
# 
# 	if np.amax(factoryPct[c])!=0:
# 		vars()[country+'Pct'] = factoryPct[c]
# 		factoryCountries.append(country)
# 
# countryComparison = np.zeros(shape=(len(factoryCountries)))
# for i in range(len(factoryCountries)):
# 	country = factoryCountries[i]
# 	countryComparison[i] = vars()[country+'Pct'][-1]
# sIndex=np.argsort(countryComparison)
# factoryCountries=np.array(factoryCountries)
# factoryCountries=factoryCountries[sIndex][::-1]
# 
# 
# colors = ['firebrick','navy','m','lawngreen','darkorange','darkgreen','crimson','yellow','mediumseagreen','darkslategray','indianred','goldenrod','blueviolet','lightseagreen','powderblue','mediumpurple','navajowhite','peru','tomato','midnightblue','magenta','forestgreen']
# width = 1.5
# for i in range(len(factoryCountries)):
# 	country = factoryCountries[i]
# 
# 	if i==0:
# 		vars()['p'+str(i)] = plt.bar(x, vars()[country+'Pct'], width, color=colors[i], )
# 		bottomStuff = vars()[country+'Pct']
# 	else:
# 		vars()['p'+str(i)] = plt.bar(x, vars()[country+'Pct'], width, color=colors[i], bottom = bottomStuff)
# 		bottomStuff+=vars()[country+'Pct']
# 
# 	print vars()[country+'Pct']
# 
# plt.title('Percent Each Facotry Produces by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('% of Total Production')
# plt.legend(([p0[0],p1[0],p2[0],p3[0],p4[0],p5[0],p6[0],p7[0],p8[0],p9[0],p10[0],p11[0],p12[0],p13[0],p14[0],p15[0],p16[0],p17[0],p18[0],p19[0]][::-1]),(factoryCountries[::-1]))
# plt.savefig(wdfigs+'factoryPct_vs_tariffCost.pdf')
# 
# 
# ##############################################
# # startup
# ##############################################
# cost=np.zeros(shape=(11))
# factoryNum=np.zeros(shape=(11))
# pctTrans=np.zeros(shape=(11))
# pctIngredient=np.zeros(shape=(11))
# sizeAvg=np.zeros(shape=(11))
# sizeStdDev=np.zeros(shape=(11))
# sizeMinMax=np.zeros(shape=(11,2))
# factorySizeAll=np.zeros(shape=(11,24))
# factorySizeAllMask=np.ones(shape=(11,24),dtype=bool)
# 
# factoryPct=np.zeros(shape=(11,43))
# #factoryPctMask=np.ones(shape=(14,43),dtype=bool)
# 
# i=-1
# for s in np.arange(0.0,0.31,0.03):
# 	s=np.round(s,2)
# 	i+=1
# 	countriesWfactories=[]
# 	factorySize=[]
# 
# 	f = open(wddata+'tariff_multiplier/tarriflevel'+str(s)+'.csv')
# 	k=-1
# 	j=-1
# 	for line in f:
# 		k+=1
# 		tmp=line.split(',')
# 		if k==0:
# 			cost[i]=float(tmp[1])
# 		elif k==1:
# 			factoryNum[i]=float(tmp[1])
# 		elif k==2:
# 			pctTrans[i]=float(tmp[1])
# 		elif k==3:
# 			pctIngredient[i]=float(tmp[1])
# 		else:
# 			j+=1
# 			country=tmp[0]
# 			country=country.replace('_',' ')
# 			countriesWfactories.append(country)
# 			size=float(tmp[1])
# 			factorySize.append(size)
# 			factorySizeAll[i,j]=size
# 			factorySizeAllMask[i,j]=False
# 			
# 			c=np.where(country==subsaharancountry)[0][0]
# 			factoryPct[i,c]=size
# 
# 	factorySize=np.array(factorySize)
# 	sizeAvg[i]=np.mean(factorySize)
# 	sizeStdDev[i]=np.std(factorySize)
# 	sizeMinMax[i,0]=np.amin(factorySize)
# 	sizeMinMax[i,1]=np.amax(factorySize)
# factorySizeAll=np.ma.masked_array(factorySizeAll,factorySizeAllMask)
# 
# totalCapacity = np.sum(factorySizeAll,axis=1)[0]
# factoryPct = np.swapaxes(factoryPct,0,1)
# factoryPct = 100*factoryPct/totalCapacity
# 
# x = np.array([0.0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3])
# x = x*100
# fig = plt.figure(figsize=(9, 6))
# #### cost ####
# plt.clf()
# plt.plot(x,cost,'*-')
# plt.title('Total Cost by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('Total Procurement Cost for One Year')
# plt.grid(True)
# plt.savefig(wdfigs+'totalCost_vs_tariffCost.pdf')
# 
# #### factoryNum ####
# plt.clf()
# plt.plot(x,factoryNum,'*-')
# plt.title('Number of Factories by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('Number of Factories')
# plt.grid(True)
# plt.savefig(wdfigs+'factoryNum_vs_tariffCost.pdf')
# 
# #### factorySize ####
# plt.clf()
# plt.plot(x,sizeAvg,'b*-')
# plt.plot(x,sizeMinMax[:,0],'b*--')
# plt.plot(x,sizeMinMax[:,1],'b*--')
# plt.title('Avg Factory Size by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('Avg Factory Size')
# plt.grid(True)
# plt.savefig(wdfigs+'factorySize_vs_tariffCost.pdf')
# 
# #### factorySizeAll ####
# plt.clf()
# plt.plot(x,factorySizeAll,'b*')
# plt.title('Factory Size by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('Factory Size')
# plt.grid(True)
# plt.savefig(wdfigs+'factorySize_vs_tariffCost.pdf')
# 
# #### PctTrans ####
# plt.clf()
# plt.plot(x,pctTrans*100,'b*-')
# plt.title('% Cost that is Transport by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('% of Total Cost that is Transport')
# plt.grid(True)
# plt.savefig(wdfigs+'pctTrans_vs_tariffCost.pdf')
# 
# #### PctIngredient ####
# plt.clf()
# plt.plot(x,pctIngredient*100,'b*-')
# plt.title('% Cost that is Ingredient by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('% of Total Cost that is Ingredient')
# plt.grid(True)
# plt.savefig(wdfigs+'pctIngredient_vs_tariffCost.pdf')
# 
# #### FactoryPct ####
# factoryCountries = []
# plt.clf()
# fig = plt.figure(figsize=(15, 10))
# 
# for i in range(len(subsaharancountry)):
# 	country = subsaharancountry[i]
# 	c=np.where(country==subsaharancountry)[0][0]
# 	if country=='Congo (Republic of the)':
# 		country='RepOfCongo'
# 	elif country=='Congo':
# 		country='DRC'
# 	country=country.replace(' ','')
# 
# 	if np.amax(factoryPct[c])!=0:
# 		vars()[country+'Pct'] = factoryPct[c]
# 		factoryCountries.append(country)
# 
# countryComparison = np.zeros(shape=(len(factoryCountries)))
# for i in range(len(factoryCountries)):
# 	country = factoryCountries[i]
# 	countryComparison[i] = vars()[country+'Pct'][-1]
# sIndex=np.argsort(countryComparison)
# factoryCountries=np.array(factoryCountries)
# factoryCountries=factoryCountries[sIndex][::-1]
# 
# 
# colors = ['firebrick','navy','m','lawngreen','darkorange','darkgreen','crimson','yellow','mediumseagreen','darkslategray','indianred','goldenrod','blueviolet','lightseagreen','powderblue','mediumpurple','navajowhite','peru','tomato','midnightblue','magenta','forestgreen']
# width = 1.5
# for i in range(len(factoryCountries)):
# 	country = factoryCountries[i]
# 
# 	if i==0:
# 		vars()['p'+str(i)] = plt.bar(x, vars()[country+'Pct'], width, color=colors[i], )
# 		bottomStuff = vars()[country+'Pct']
# 	else:
# 		vars()['p'+str(i)] = plt.bar(x, vars()[country+'Pct'], width, color=colors[i], bottom = bottomStuff)
# 		bottomStuff+=vars()[country+'Pct']
# 
# 	print vars()[country+'Pct']
# 
# plt.title('Percent Each Facotry Produces by Tariff Cost')
# plt.xlabel('Tariff Cost, % of Today')
# plt.ylabel('% of Total Production')
# plt.legend(([p0[0],p1[0],p2[0],p3[0],p4[0],p5[0],p6[0],p7[0],p8[0],p9[0],p10[0],p11[0],p12[0],p13[0],p14[0],p15[0],p16[0],p17[0],p18[0],p19[0]][::-1]),(factoryCountries[::-1]))
# plt.savefig(wdfigs+'factoryPct_vs_tariffCost.pdf')
# 
# 
# 
# 
