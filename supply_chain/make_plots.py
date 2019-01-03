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
geolocator = Nominatim()
import geopy.distance
from scipy import ndimage
from scipy.signal import convolve2d
from sklearn import linear_model
from scipy.interpolate import RectSphereBivariateSpline
# from libtiff import TIFF
# import googlemaps
import json
from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import itertools


try:
	wddata='/Users/lilllianpetersen/iiasa/data/supply_chain/'
	wdfigs='/Users/lilllianpetersen/iiasa/figs/supply_chain/'
	wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'
	f=open(wddata+'population/CAPITALVERSIONcasenumbers.csv','r')
except:
	wddata='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/data/'
	wdfigs='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/figs/'
	wdvars='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/vars/'

MakePlots=False


subsaharancountry = np.load(wdvars+'subsaharancountry.npy')
#for f in range(len(subsaharancountry)):
#	subsaharancountry[f]=subsaharancountry[f].replace(' ','_')
subsaharancountry[subsaharancountry=='Congo']='DRC'
subsaharancountry[subsaharancountry=='Congo (Republic of the)']='Congo'
subsaharancountry[subsaharancountry=="Cote d'Ivoire"]='Ivory Coast'

countrycosted=np.load(wdvars+'countrycosted.npy')
countrycosted[countrycosted=='Congo']='DRC'
countrycosted[countrycosted=='Congo (Republic of the)']='Congo'
countrycosted[countrycosted=="I_Cote d'Ivoire"]='I_Ivory Coast'
countrycosted[countrycosted=="Cote d'Ivoire"]='Ivory Coast'

capitalcosted=np.load(wdvars+'capitalcosted.npy')
subsaharancapital=np.load(wdvars+'subsaharancapital.npy')

try:
	capitalcostedLatLon = np.load(wdvars+'capitalLatLon.npy')
except:
	capitalLatLon=np.zeros(shape=(2,len(capitalcosted)))
	for c in range(len(capitalcosted)):
		if countrycosted[c][:2]=='I_':
			location = geolocator.geocode(capitalcosted[c]+', '+countrycosted[c][2:])
		else:
			location = geolocator.geocode(capitalcosted[c]+', '+countrycosted[c])
		capitalLatLon[0,c] = location.latitude
		capitalLatLon[1,c] = location.longitude
	np.save(wdvars+'capitalLatLon.npy',capitalLatLon)

try:
	SScapitalLatLon = np.load(wdvars+'subsaharancapitalLatLon.npy')
except:
	SScapitalLatLon=np.zeros(shape=(2,len(subsaharancapital)))
	for c in range(len(subsaharancountry)):
		location = geolocator.geocode(subsaharancapital[c]+', '+subsaharancountry[c])
		SScapitalLatLon[0,c] = location.latitude
		SScapitalLatLon[1,c] = location.longitude
	np.save(wdvars+'subsaharancapitalLatLon.npy',SScapitalLatLon)

optiLevel = ['AllarOpti','LocalOpti', 'AllIntl_opti', 'AllIntl']
trfLevel = ['','_trf']
loopvar = ['shipcost','impexp','strtup','truckfactor']

LTitles = ['All Optimized','Local Optimized','Optimized Intl','Current Intl']
TTitles = [' (No Tariff)',' (Tariff)']
VTitles = ['Shipping','Import/Export','Startup','Trucking']

Ltitles = ['AllOpti','LocalOpti', 'AllIntlOpti', 'AllIntl']
Ttitles = ['_noTrf','_trf']
Vtitles = ['shipping','importexport','startup','trucking']

factor = np.array([0.2,0.2,0.5,0.2])
maxs = np.array([2.01,4.01,9.51,4.01])
#LTitles = ['All Optimized','Local Optimized','Local Producing Optimized International with Tariff','Local Producing International with Tariff','Local Optimized with Tariff','All Optimized with Tariff','Local Producing Optimized International']
#optiLevel = ['AllarOpti','LocalOpti', 'AllIntl_opti_trf', 'AllIntl_trf', 'LocalOpti_trf', 'AllarOpti_trf', 'AllIntl_opti']

#[AllOpti,LocalOpti,AllIntl,AllIntlOpti] , [Trf,NoTrf] , [Shipcost,imexp,strtup,truckfactor] , [scenarios]
cost=np.zeros(shape=(4,2,4,20)) 
Mask=np.ones(shape=(4,2,4,20),dtype=bool) 
factoryNum=np.zeros(shape=(4,2,4,20))
pctLocal=np.zeros(shape=(4,2,4,20,2))

for L in range(len(optiLevel)):
	for T in range(len(trfLevel)):
		for V in range(len(loopvar)):
			File =optiLevel[L]+trfLevel[T]+'_'+loopvar[V]
			print optiLevel[L],trfLevel[T],loopvar[V]
	
			shp=len(np.arange(factor[V],maxs[V],factor[V]))
			pctTrans=np.zeros(shape=(shp))
			pctIngredient=np.zeros(shape=(shp))
			avgShipments=np.zeros(shape=(shp))
			sizeAvg=np.zeros(shape=(2,shp))
			sizeStdDev=np.zeros(shape=(2,shp))
			sizeMinMax=np.zeros(shape=(2,shp,2))
			factorySizeAll=np.zeros(shape=(2,shp,33))
			factorySizeAllMask=np.ones(shape=(2,shp,33),dtype=bool)
			
			factoryPct=np.zeros(shape=(2,shp,len(countrycosted)))
			factoryPctOne=np.zeros(shape=(2,len(countrycosted)))
	
			for f in range(len(countrycosted)):
				countrycosted[f]=countrycosted[f].replace(' ','_')
	
			i=-1
			for s in np.arange(factor[V],maxs[V],factor[V]):
				s=np.round(s,2)
				#print '\n',s
				i+=1
				countriesWfactories=[]
				factorySizeS=[]
				factorySizeM=[]
				capacity=np.zeros(shape=2)
			
				try:
					f = open(wddata+'results/'+str(File)+'/'+str(File)+str(s)+'.csv')
				except:
					continue
				k=-1
				j=-1
				for line in f:
					k+=1
					tmp=line.split(',')
					if k==0:
						cost[L,T,V,i]=float(tmp[1])
						Mask[L,T,V,i]=0
						if s==1:
							costOne=float(tmp[1])
						#print 'cost',cost[i]
					elif k==1:
						factoryNum[L,T,V,i]=float(tmp[1])
						#print 'factoryNum',factoryNum[i]
					elif k==2:
						pctTrans[i]=float(tmp[1])
						#print 'pctTrans',pctTrans[i]
					elif k==3:
						pctIngredient[i]=float(tmp[1])
						#print 'pctIngredient',pctIngredient[i]
					elif k==4:
						avgShipments[i]=float(tmp[1])
						#print 'avgShipments',avgShipments[i]
					else:
						j+=1
						country=tmp[0]
						if country=='Congo': country='DRC'
						if country=='Congo_(Republic_of_the)': country='Congo'
						if country=="I_Cote_d'Ivoire": country='I_Ivory_Coast'
						if country=="Cote_d'Ivoire": country='Ivory_Coast'

						c=np.where(country==countrycosted)[0][0]
	
						#country=country.replace('_',' ')
						countriesWfactories.append(country)
						factorySizeS.append(float(tmp[1]))
						factorySizeM.append(float(tmp[2]))
						factorySizeAll[0,i,c]=float(tmp[1])
						factorySizeAll[1,i,c]=float(tmp[2])
						factorySizeAllMask[:,i,c]=False
						capacity[0]+=float(tmp[1])
						capacity[1]+=float(tmp[2])
						
						factoryPct[0,i,c]=float(tmp[1])
						factoryPct[1,i,c]=float(tmp[2])
						if s==1:
							factoryPctOne[0,c]=float(tmp[1])
							factoryPctOne[1,c]=float(tmp[2])
						if country[:2]!='I_':
							pctLocal[L,T,V,i,0]+=float(tmp[1])
							pctLocal[L,T,V,i,1]+=float(tmp[2])

				pctLocal[L,T,V,i,:]=pctLocal[L,T,V,i,:]/capacity[:]
					
				factorySizeS=np.array(factorySizeS)
				factorySizeM=np.array(factorySizeM)
				sizeAvg[0,i]=np.mean(factorySizeS)
				sizeAvg[1,i]=np.mean(factorySizeM)
				sizeStdDev[0,i]=np.std(factorySizeS)
				sizeStdDev[1,i]=np.std(factorySizeM)
				sizeMinMax[0,i,0]=np.amin(factorySizeS)
				sizeMinMax[1,i,0]=np.amin(factorySizeM)
				sizeMinMax[0,i,1]=np.amax(factorySizeS)
				sizeMinMax[1,i,1]=np.amax(factorySizeM)
			factorySizeAll=np.ma.masked_array(factorySizeAll,factorySizeAllMask)
			
			totalCapacity = np.zeros(shape=(2,len(factorySizeAll[0])))
			totalCapacity[0] = np.sum(factorySizeAll[0],axis=1)
			totalCapacity[1] = np.sum(factorySizeAll[1],axis=1)
			factoryPct = np.swapaxes(factoryPct,1,2)
			for p in range(33):
				factoryPct[0,p] = 100*factoryPct[0,p]/totalCapacity[0]
				factoryPct[1,p] = 100*factoryPct[1,p]/totalCapacity[1]
	
			if not os.path.exists(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Vtitles[V]):
					os.makedirs(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Vtitles[V])

			x = np.arange(factor[V],maxs[V],factor[V])
			x = x*100
			if MakePlots:
				fig = plt.figure(figsize=(9, 6))
				#### cost ####
				ydata=np.ma.compressed(np.ma.masked_array(cost[L,T,V,:],Mask[L,T,V,:]))
				plt.clf()
				plt.plot(x,ydata,'b*-')
				plt.title(LTitles[L]+TTitles[T]+': Effect of '+VTitles[V]+' Cost on Total Cost')
				plt.xlabel(VTitles[V]+' Cost, % of Today')
				plt.ylabel('Total Procurement Cost for One Year')
				plt.grid(True)
				plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Vtitles[V]+'/'+Ltitles[L]+Ttitles[T]+'__totalCost_vs_'+Vtitles[V]+'.pdf')
			
				#### factoryNum ####
				ydata=np.ma.compressed(np.ma.masked_array(factoryNum[L,T,V,:],Mask[L,T,V,:]))
				plt.clf()
				plt.plot(x,ydata,'b*-')
				plt.title(LTitles[L]+TTitles[T]+': Effect of '+VTitles[V]+' Number of Factories')
				plt.xlabel(VTitles[V]+' Cost, % of Today')
				plt.ylabel('Number of Factories')
				plt.grid(True)
				plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Vtitles[V]+'/'+Ltitles[L]+Ttitles[T]+'__factoryNum_vs_'+Vtitles[V]+'.pdf')
				
				#### factorySize ####
				plt.clf()
				plt.plot(x,sizeAvg[0],'b*-')
				plt.plot(x,sizeAvg[0]-sizeStdDev[0],'b*--')
				plt.plot(x,sizeAvg[0]+sizeStdDev[0],'b*--')
				plt.plot(x,sizeAvg[1],'g*-')
				plt.plot(x,sizeAvg[1]-sizeStdDev[1],'g*--')
				plt.plot(x,sizeAvg[1]+sizeStdDev[1],'g*--')
				plt.title(LTitles[L]+': Avg Factory Size by '+VTitles[V]+' Cost')
				plt.xlabel(VTitles[V]+' Cost, % of Today')
				plt.ylabel('Avg Factory Size')
				plt.grid(True)
				plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Vtitles[V]+'/'+Ltitles[L]+'factorySize_vs_'+Vtitles[V]+'.pdf')
				
				#### factorySizeAll ####
				for g in range(2):
					plt.clf()
					plt.plot(x,factorySizeAll[g],'b*')
					plt.title(LTitles[L]+': Factory Size by '+VTitles[V]+' Cost')
					plt.xlabel(''+VTitles[V]+' Cost, % of Today')
					plt.ylabel('Factory Size')
					plt.grid(True)
					plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Vtitles[V]+'/'+Ltitles[L]+'factorySize'+str(g)+'_vs_'+Vtitles[V]+'Cost.pdf')
				
				#### PctTrans ####
				plt.clf()
				plt.plot(x,pctTrans*100,'b*-')
				plt.title(LTitles[L]+': % Cost that is Transport by '+VTitles[V]+' Cost')
				plt.xlabel(''+VTitles[V]+' Cost, % of Today')
				plt.ylabel('% of Total Cost that is Transport')
				plt.grid(True)
				plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Vtitles[V]+'/'+Ltitles[L]+'pctTrans_vs_'+Vtitles[V]+'Cost.pdf')
				
				#### PctIngredient ####
				plt.clf()
				plt.plot(x,pctIngredient*100,'b*-')
				plt.title(LTitles[L]+': % Cost that is Ingredient by '+VTitles[V]+' Cost')
				plt.xlabel(VTitles[V]+' Cost, % of Today')
				plt.ylabel('% of Total Cost that is Ingredient')
				plt.grid(True)
				plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Vtitles[V]+'/'+Ltitles[L]+'pctIngredient_vs_'+Vtitles[V]+'Cost.pdf')
			
			#### FactoryPct ####
			SMtitles=['SAM','MAM']
			for g in range(2):
				factoryCountries = []
				plt.clf()
				fig = plt.figure(figsize=(18, 7))
				ax = plt.subplot(1,2,1)
				
				for t in range(len(countrycosted)):
					country = countrycosted[t]
					c=np.where(country==countrycosted)[0][0]
					if country=='Congo (Republic of the)':
						country='RepOfCongo'
					elif country=='Congo':
						country='DRC'
					country=country.replace(' ','_')
				
					if np.amax(factoryPct[g,c])!=0:
						vars()[country+'Pct'] = factoryPct[g,c]
						factoryCountries.append(country)
				
				countryComparison = np.zeros(shape=(len(factoryCountries)))
				for q in range(len(factoryCountries)):
					country = factoryCountries[q]
					countryComparison[q] = vars()[country+'Pct'][-1]
				sIndex=np.argsort(countryComparison)
				factoryCountries=np.array(factoryCountries)
				factoryCountries=factoryCountries[sIndex][::-1]
				
				OfactoryCountries = []
				Otitles = []
				for country in factoryCountries:
					if country[:2]!='I_':
						OfactoryCountries.append(country)
						Otitles.append(country+' Factory')
				for country in factoryCountries:
					if country[:2]=='I_':
						OfactoryCountries.append(country)
						countryTitle = 'Intl: '+country[2:]+' Port'
						Otitles.append(countryTitle)
				
				if MakePlots:
					Dcolors = ['firebrick','m','darkorange','crimson','yellow','indianred','goldenrod','mediumpurple','navajowhite','peru','tomato','magenta','deeppink','lightcoral','lemonchiffon','sandybrown','r','gold','moccasin','peachpuff','orangered','orange','rosybrown','papayawhip']
					Icolors = ['navy','lawngreen','darkgreen','deepskyblue','darkslategray','mediumseagreen','lightseagreen','powderblue','midnightblue','forestgreen']
					width = 7 
					pvars=[]
					inter=-1
					domes=-1
					for l in range(len(OfactoryCountries)):
						country = OfactoryCountries[l]
						if country[:2]=='I_':
							inter+=1
							clr=Icolors[inter]
						else:
							domes+=1
							clr=Dcolors[domes]
					
						if l==0:
							vars()['p'+str(l)] = ax.bar(x, vars()[country+'Pct'], width, color=clr, )
							bottomStuff = vars()[country+'Pct']
						else:
							vars()['p'+str(l)] = ax.bar(x, vars()[country+'Pct'], width, color=clr, bottom = bottomStuff)
							bottomStuff+=vars()[country+'Pct']
					
						pvars.append(vars()['p'+str(l)])
						
					
					fontP = FontProperties()
					fontP.set_size('small')
		
					plt.title(LTitles[L]+':'+ SMtitles[g]+' Treatment: Factory Production by '+VTitles[V]+' Cost')
					plt.xlabel(''+VTitles[V]+' Cost, % of Today')
					plt.ylabel('% of Total Production')
					ax.legend((pvars[::-1]),(Otitles[::-1]),bbox_to_anchor=(0.97, 0.6),prop=fontP)
					plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Vtitles[V]+'/'+Ltitles[L]+'factoryPct'+str(g)+'_vs_'+Vtitles[V]+'Cost.pdf')

			##################################################################
			# MAPS
			##################################################################

			##### Export #####
			countrycosted=np.load(wdvars+'countrycosted.npy')
			countrycosted[countrycosted=='Congo']='DRC'
			countrycosted[countrycosted=='Congo (Republic of the)']='Congo'
			countrycosted[countrycosted=="I_Cote d'Ivoire"]='I_Ivory Coast'
			countrycosted[countrycosted=="Cote d'Ivoire"]='Ivory Coast'

			shapename = 'admin_0_countries'
			countries_shp = shpreader.natural_earth(resolution='110m',
				category='cultural', name=shapename)
			
			plt.clf()
			cmapArray=plt.cm.hot_r(np.arange(256))
			cmin=0
			cmax=np.amax(factoryPctOne[0,:]) #*0.9
			y1=0
			y2=255
			
			fig = plt.figure(figsize=(10, 8))
			MinMaxArray=np.ones(shape=(3,2))
			subPlot1 = plt.axes([0.61, 0.07, 0.2, 0.8])
			MinMaxArray[0,0]=cmin
			MinMaxArray[1,0]=cmax
			plt.imshow(MinMaxArray,cmap='hot_r')
			plt.colorbar()
			plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Ltitles[L]+'_export_map.pdf')
			
			ax = plt.axes([0.05,0.05,0.8,0.85],projection=ccrs.PlateCarree())
			ax.set_extent([-19, 53, -37, 39], ccrs.PlateCarree())
			ax.coastlines()

			plt.plot(capitalLatLon[1,8], capitalLatLon[0,8], marker='*', markersize=7, color='g',label='Factories')
			plt.plot(capitalLatLon[1,8], capitalLatLon[0,8], marker='*', markersize=5, color='darkred', label='Possible Factories (Not Producing)')
			plt.plot(capitalLatLon[1,30], capitalLatLon[0,30], marker='o', markersize=7, color='g', label = 'Intl Shipment Port')
			plt.plot(capitalLatLon[1,30], capitalLatLon[0,30], marker='o', markersize=5, color='darkred', label = 'Intl Shipment Port (No Shipments)')

			factoryNumOne=0
			IntlNumOne=0
			
			for country in shpreader.Reader(countries_shp).records():
				cName=country.attributes['NAME_LONG']
				if cName[-6:]=='Ivoire':
					cName="Ivory Coast"
				if cName=='Democratic Republic of the Congo':
					cName='DRC'
				if cName=='Republic of the Congo':
					cName='Congo'
				if cName=='eSwatini':
					cName='Swaziland'
				if cName=='The Gambia':
					cName='Gambia'
				if np.amax(cName==subsaharancountry)==0:
					continue
				if np.amax(cName==countrycosted)==0:
					x=0
					y=y1+(y2-y1)/(cmax-cmin)*(x-cmin)
					icmap=min(255,int(round(y,1)))
					icmap=max(0,int(round(icmap,1)))
					ax.add_geometries(country.geometry, ccrs.PlateCarree(), edgecolor='black',
						facecolor=[cmapArray[icmap,0],cmapArray[icmap,1],cmapArray[icmap,2]],label=cName)
				else:
					c=np.where(cName==countrycosted)[0][0]
					x=factoryPctOne[0,c]
					y=y1+(y2-y1)/(cmax-cmin)*(x-cmin)
					icmap=min(255,int(round(y,1)))
					icmap=max(0,int(round(icmap,1)))
					ax.add_geometries(country.geometry, ccrs.PlateCarree(), edgecolor='black', facecolor=[cmapArray[icmap,0],cmapArray[icmap,1],cmapArray[icmap,2]],label=cName)

					if x!=0:
						size = 10*(1+factoryPctOne[0,c]/cmax)
						plt.plot(capitalLatLon[1,c], capitalLatLon[0,c], marker='*', markersize=size, color='g')
						factoryNumOne+=1
					if x==0:
						plt.plot(capitalLatLon[1,c], capitalLatLon[0,c], marker='*', markersize=7, color='darkred')

			
			for icoast in range(24,len(countrycosted)):
				x=factoryPctOne[0,icoast]
				if x!=0:
					size = 10*(1+factoryPctOne[0,icoast]/cmax)
					plt.plot(capitalLatLon[1,icoast], capitalLatLon[0,icoast], marker='o', markersize=size, color='g')
					IntlNumOne+=1
				if x==0:
					plt.plot(capitalLatLon[1,icoast], capitalLatLon[0,icoast], marker='o', markersize=7, color='darkred')

			local = str(int(np.round(100*np.sum(factoryPctOne[0,:24])/np.sum(factoryPctOne[0,:]),0)))
			intl = str(np.round(100*np.sum(factoryPctOne[0,24:])/np.sum(factoryPctOne[0,:]),0))
			costOne = str(int(round(costOne/1000000.,0)))

			plt.title('Exports of RUTF by Factory, Packets \n' + LTitles[L] + TTitles[T])
			plt.legend(loc = 'lower left')
			plt.text(-15,-10,str(factoryNumOne)+' Factories Open\n'+str(IntlNumOne)+' Ports Open\n'+local+'% Produced Locally\nTotal Cost = $'+costOne+' Million', bbox=dict(fc="none", boxstyle="round"), size = 10)
			plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Ltitles[L]+'_export_map.pdf')

			##### Skeleton #####
			shapename = 'admin_0_countries'
			countries_shp = shpreader.natural_earth(resolution='110m',
				category='cultural', name=shapename)
			
			plt.clf()
			ax = plt.axes([0.05,0.05,0.8,0.85],projection=ccrs.PlateCarree())
			ax.set_extent([-19, 53, -37, 39], ccrs.PlateCarree())
			ax.coastlines()

			plt.plot(capitalLatLon[1,8], capitalLatLon[0,8], marker='*', markersize=8, color='slateblue', label='Possible Factories')
			plt.plot(capitalLatLon[1,30], capitalLatLon[0,30], marker='o', markersize=8, color='c', label = 'Possible Intl Shipment Ports')
			
			for country in shpreader.Reader(countries_shp).records():
				cName=country.attributes['NAME_LONG']
				if cName[-6:]=='Ivoire':
					cName="Ivory Coast"
				if cName=='Democratic Republic of the Congo':
					cName='DRC'
				if cName=='Republic of the Congo':
					cName='Congo'
				if cName=='eSwatini':
					cName='Swaziland'
				if cName=='The Gambia':
					cName='Gambia'
				if np.amax(cName==subsaharancountry)==0:
					continue
				plt.plot(
				if np.amax(cName==countrycosted)==0:
					ax.add_geometries(country.geometry, ccrs.PlateCarree(), edgecolor='black',
						facecolor='lightgray')
				else:
					c=np.where(cName==countrycosted)[0][0]
					ax.add_geometries(country.geometry, ccrs.PlateCarree(), edgecolor='black', facecolor='lightgreen',label=cName)

					plt.plot(capitalLatLon[1,c], capitalLatLon[0,c], marker='*', markersize=8, color='slateblue')

			
			for icoast in range(24,len(countrycosted)):
				x=factoryPctOne[0,icoast]
				plt.plot(capitalLatLon[1,icoast], capitalLatLon[0,icoast], marker='o', markersize=8, color='c')

			#for l in range(len(

			plt.title('Multi-Product Capacitated Facility Location Model\nPossible Factory Locations')
			plt.legend(loc = 'lower left')
			plt.text(-15,-10,'24 Possible Factories\n9 Possible Ports', bbox=dict(fc="none", boxstyle="round"), size = 10)
			plt.savefig(wdfigs+'skeleton_map.pdf')
			exit()

exit()


fig = plt.figure(figsize=(9, 6))
LTitles = ['All Optimized','Local Optimized','Optimized Intl','None Optimized']

for T in range(len(trfLevel)):
	for V in range(len(loopvar)):
		x = np.arange(factor[V],maxs[V],factor[V])
		x = x*100
		plt.clf()
		plt.plot(x,np.ma.compressed(np.ma.masked_array(cost[0,T,V,:],Mask[0,T,V,:])),'g*-',label=LTitles[0])
		plt.plot(x,np.ma.compressed(np.ma.masked_array(cost[1,T,V,:],Mask[1,T,V,:])),'c*-',label=LTitles[1])
		plt.plot(x,np.ma.compressed(np.ma.masked_array(cost[2,T,V,:],Mask[2,T,V,:])),'b*-',label=LTitles[2])
		try:
			plt.plot(x,np.ma.compressed(np.ma.masked_array(cost[3,T,V,:],Mask[3,T,V,:])),'r*-',label=LTitles[3])
		except:
			print T,V

		plt.title('Effect of '+VTitles[V]+' on Total Cost'+TTitles[T])
		plt.xlabel(VTitles[V]+' Cost, % of Today')
		plt.ylabel('Total Procurement Cost for One Year')
		#plt.ylim([0,2e9])
		plt.grid(True)
		plt.legend()
		plt.savefig(wdfigs+'totalCost_vs_'+Vtitles[V]+'_'+Ttitles[T]+'.pdf')

for T in range(len(trfLevel)):
	for V in range(len(loopvar)):
		x = np.arange(factor[V],maxs[V],factor[V])
		x = x*100
		plt.clf()
		plt.plot(x,np.ma.compressed(np.ma.masked_array(factoryNum[0,T,V,:],Mask[0,T,V,:])),'g*-',label=LTitles[0])
		plt.plot(x,np.ma.compressed(np.ma.masked_array(factoryNum[1,T,V,:],Mask[1,T,V,:])),'c*-',label=LTitles[1])
		plt.plot(x,np.ma.compressed(np.ma.masked_array(factoryNum[2,T,V,:],Mask[2,T,V,:])),'b*-',label=LTitles[2])
		plt.plot(x,np.ma.compressed(np.ma.masked_array(factoryNum[3,T,V,:],Mask[3,T,V,:])),'r*-',label=LTitles[3])

		plt.title('Effect of '+VTitles[V]+' on Number of Factories'+TTitles[T])
		plt.xlabel(VTitles[V]+' Cost, % of Today')
		plt.ylabel('Number of Factories')
		#plt.ylim([0,2e9])
		plt.grid(True)
		plt.legend()
		plt.savefig(wdfigs+'factoryNum_vs_'+Vtitles[V]+'_'+Ttitles[T]+'.pdf')

for T in range(len(trfLevel)):
	for V in range(len(loopvar)):
		x = np.arange(factor[V],maxs[V],factor[V])
		x = x*100
		plt.clf()
		plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[0,T,V,:,0],Mask[0,T,V,:])),'g*-',label=LTitles[0])
		plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[1,T,V,:,0],Mask[1,T,V,:])),'c*-',label=LTitles[1])
		plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[2,T,V,:,0],Mask[2,T,V,:])),'b*-',label=LTitles[2])
		try:
			plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[3,T,V,:,0],Mask[3,T,V,:])),'r*-',label=LTitles[3])
		except:
			print T,V
		plt.title('Effect of '+VTitles[V]+' on % RUTF Produced Locally'+TTitles[T])
		plt.xlabel(VTitles[V]+' Cost, % of Today')
		plt.ylabel('Percent of RUTF Produced Locally')
		plt.ylim([0,101])
		plt.grid(True)
		plt.legend()
		plt.savefig(wdfigs+'0pctLocal_vs_'+Vtitles[V]+'_'+Ttitles[T]+'.pdf')

for T in range(len(trfLevel)):
	for V in range(len(loopvar)):
		x = np.arange(factor[V],maxs[V],factor[V])
		x = x*100
		plt.clf()
		plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[0,T,V,:,1],Mask[0,T,V,:])),'g*-',label=LTitles[0])
		plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[1,T,V,:,1],Mask[1,T,V,:])),'c*-',label=LTitles[1])
		plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[2,T,V,:,1],Mask[2,T,V,:])),'b*-',label=LTitles[2])
		try:
			plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[3,T,V,:,1],Mask[3,T,V,:])),'r*-',label=LTitles[3])
		except:
			print T,V
		plt.title('Effect of '+VTitles[V]+' on % RUSF Produced Locally'+TTitles[T])
		plt.xlabel(VTitles[V]+' Cost, % of Today')
		plt.ylabel('Percent of RUSF Produced Locally')
		plt.ylim([0,101])
		plt.grid(True)
		plt.legend()
		plt.savefig(wdfigs+'1pctLocal_vs_'+Vtitles[V]+'_'+Ttitles[T]+'.pdf')



#######################
# One One One array
####################### 
rusfarray = np.load(wddata+'results/example/Rrusfarray.npy')
rutfarray = np.load(wddata+'results/example/Rrutfarray.npy')
Rcountry = np.load(wddata+'results/example/Rcountry.npy')
Rcountry[Rcountry=='Congo']='DRC'
Rcountry[Rcountry=='Congo_(Republic_of_the)']='Congo'
Rcountry[Rcountry=='Cote_d\'Ivoire']='Ivory Coast'

Rcountry2=[]
for i in range(len(Rcountry)):
	country=Rcountry[i]
	if country[:2]=='I_':
		if country=="I_Cote_d'Ivoire":
			country='I_Ivory_Coast'
		Rcountry2.append('Intl ('+country[2:]+')')
	else:
		Rcountry2.append(country)
Rcountry2=np.array(Rcountry2)
for i in range(len(Rcountry2)):
	Rcountry2[i]=Rcountry2[i].replace('_',' ')

plt.clf()
fig = plt.figure(figsize=(13, 8))
plt.imshow(rusfarray,cmap=cm.terrain_r,vmax=2.5e8)
plt.colorbar()
plt.xticks(np.arange(43), subsaharancountry, rotation='vertical')
plt.yticks(np.arange(len(Rcountry2)), Rcountry2)
plt.title('Imports and Exports (Packets, Ordered by Lon)')
plt.savefig(wdfigs+'impexp_array.pdf')

