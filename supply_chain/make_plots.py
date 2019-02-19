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
from matplotlib.font_manager import FontProperties
# from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import itertools

################################################3
# Functions
################################################3
def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap
################################################3

try:
    wddata='/Users/lilllianpetersen/iiasa/data/supply_chain/'
    wdfigs='/Users/lilllianpetersen/iiasa/figs/supply_chain/'
    wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'
    f=open(wddata+'population/CAPITALVERSIONcasenumbers.csv','r')
except:
    wddata='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/supply_chain/data/'
    wdfigs='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/supply_chain/figs/'
    wdvars='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/supply_chain/vars/'

MakePlots=True

subsaharancountry = np.load(wdvars+'subsaharancountry.npy')


#for f in range(len(subsaharancountry)):
#    subsaharancountry[f]=subsaharancountry[f].replace(' ','_')
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
    capitalLatLon = np.load(wdvars+'capitalLatLon.npy')
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

#  'AllIntl_opti',
optiLevel = ['AllarOpti','LocalOpti', 'AllIntl']
loopvar = ['shipcost','impexp','strtup','truckfactor','tariff']

LTitles = ['All Optimized','Local Optimized','Optimized Intl','Current Intl']
TTitles = [' (No Tariff)',' (Tariff)']
VTitles = ['Shipping','Import/Export','Startup','Trucking', 'Tariff']

Ltitles = ['AllOpti','LocalOpti', 'AllIntlOpti', 'AllIntl']
Ttitles = ['_noTrf','_trf']
Vtitles = ['shipping','importexport','startup','trucking', 'tariff']

mins= np.array([0.2,0.2,0.5,0.2,0])
factor = np.array([0.2,0.2,0.5,0.2,0.2])
maxs = np.array([2.01,4.01,9.51,4.01,2.6])
#LTitles = ['All Optimized','Local Optimized','Local Producing Optimized International with Tariff','Local Producing International with Tariff','Local Optimized with Tariff','All Optimized with Tariff','Local Producing Optimized International']
#optiLevel = ['AllarOpti','LocalOpti', 'AllIntl_opti_trf', 'AllIntl_trf', 'LocalOpti_trf', 'AllarOpti_trf', 'AllIntl_opti']

#[AllOpti,LocalOpti,AllIntl,AllIntlOpti] , [Trf,NoTrf] , [Shipcost,imexp,strtup,truckfactor] , [scenarios]
cost=np.zeros(shape=(4,1,5,20)) 
costOneAll=np.zeros(shape=(4,1)) 
Mask=np.ones(shape=(4,1,5,20),dtype=bool) 
factoryNum=np.zeros(shape=(4,1,5,20))
factoryNumOneAll=np.zeros(shape=(4,1))
portNumOneAll=np.zeros(shape=(4,1))
pctLocal=np.zeros(shape=(4,1,5,20,2))
pctLocalOneAll=np.zeros(shape=(4,1,2)) #optiLevel, trf, SAM-MAM

for L in range(len(optiLevel)):
    for T in range(1):
        for V in range(len(loopvar)):
            File =optiLevel[L]+'_'+loopvar[V]
            print optiLevel[L],loopvar[V]
    
            shp=len(np.arange(mins[V],maxs[V],factor[V]))
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
            capacityOne=np.zeros(shape=2)
    
            for f in range(len(countrycosted)):
                countrycosted[f]=countrycosted[f].replace(' ','_')
    
            i=-1
            for s in np.arange(mins[V],maxs[V],factor[V]):
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
                            costOneAll[L,T]=float(tmp[1])
                        #print 'cost',cost[i]
                    elif k==1:
                        factoryNum[L,T,V,i]=float(tmp[1])
                        if s==1:
                            factoryNumOneAll[L,T]=float(tmp[1])
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
                            capacityOne[0]+=float(tmp[1])
                            capacityOne[1]+=float(tmp[2])
                            if country[:2]!='I_':
                                pctLocalOneAll[L,T,0]+=float(tmp[1])
                                pctLocalOneAll[L,T,1]+=float(tmp[2])
                            if country[:2]=='I_' and V==0:
                                portNumOneAll[L,T]+=1
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
            pctLocalOneAll[L,T,:]=pctLocalOneAll[L,T,:]/capacityOne[:]
            
            totalCapacity = np.zeros(shape=(2,len(factorySizeAll[0])))
            totalCapacity[0] = np.sum(factorySizeAll[0],axis=1)
            totalCapacity[1] = np.sum(factorySizeAll[1],axis=1)
            factoryPct = np.swapaxes(factoryPct,1,2)
            for p in range(33):
                factoryPct[0,p] = 100*factoryPct[0,p]/totalCapacity[0]
                factoryPct[1,p] = 100*factoryPct[1,p]/totalCapacity[1]
    
            if not os.path.exists(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Vtitles[V]):
                os.makedirs(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Vtitles[V])
            if not os.path.exists(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/geographical'):
                os.makedirs(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/geographical')
            if not os.path.exists(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/exports_by_country/'):
                os.makedirs(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/exports_by_country/')

            x = np.arange(mins[V],maxs[V],factor[V])
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
                
                    if np.amax(factoryPct[:,c])!=0:
                        vars()[country+'Pct'] = factoryPct[g,c]
                        factoryCountries.append(country)
                
                if g==0:
                    countryComparison = np.zeros(shape=(len(factoryCountries)))
                    for q in range(len(factoryCountries)):
                        country = factoryCountries[q]
                        countryComparison[q] = vars()[country+'Pct'][5]
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
        
                    plt.title(SMtitles[g]+' Treatment: Factory Production by '+VTitles[V]+' Cost ('+LTitles[L]+')')
                    plt.xlabel(''+VTitles[V]+' Cost, % of Today')
                    plt.ylabel('% of Total Production')
                    if g==1:
                        ax.legend((pvars[::-1]),(Otitles[::-1]),bbox_to_anchor=(1, 0.98),prop=fontP)
                    plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Vtitles[V]+'/'+Ltitles[L]+'factoryPct'+str(g)+'_vs_'+Vtitles[V]+'Cost.pdf')

        ##################################################################
        # MAPS
        ##################################################################

        countrycosted=np.load(wdvars+'countrycosted.npy')
        countrycosted[countrycosted=='Congo']='DRC'
        countrycosted[countrycosted=='Congo (Republic of the)']='Congo'
        countrycosted[countrycosted=="I_Cote d'Ivoire"]='I_Ivory Coast'
        countrycosted[countrycosted=="Cote d'Ivoire"]='Ivory Coast'

        ###########################
        # Export
        ###########################
        if MakePlots:
            colors = [(255,255,255),(152, 240, 152), (97, 218, 97), (65, 196, 65), (42, 175, 42), (28, 162, 28), (17, 149, 17), (7, 135, 7), (0, 118, 0)]
            my_cmap = make_cmap(colors,bit=True)
            shapename = 'admin_0_countries'
            countries_shp = shpreader.natural_earth(resolution='110m',
                category='cultural', name=shapename)
            
            for g in range(2): # SAM MAM treatment
                plt.clf()
                cmapArray=my_cmap(np.arange(256))
                cmin=0
                cmax=np.amax(factoryPctOne[g,:]) #*0.9
                y1=0
                y2=255
                
                fig = plt.figure(figsize=(10, 8))
                MinMaxArray=np.ones(shape=(3,2))
                subPlot1 = plt.axes([0.61, 0.07, 0.2, 0.8])
                MinMaxArray[0,0]=cmin
                MinMaxArray[1,0]=cmax
                plt.imshow(MinMaxArray,cmap=my_cmap)
                plt.colorbar()
                plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Ltitles[L]+'_export_map.pdf')
                
                ax = plt.axes([0.05,0.05,0.8,0.85],projection=ccrs.PlateCarree())
                ax.set_extent([-19, 53, -37, 39], ccrs.PlateCarree())
                ax.coastlines()
        
                plt.plot(capitalLatLon[1,8], capitalLatLon[0,8], marker='*', markersize=12, color=[97/255., 218/255., 97/255.], markeredgewidth=1.5, markeredgecolor='k',label='Factories')
                plt.plot(capitalLatLon[1,8], capitalLatLon[0,8], marker='*', markersize=7, color='darkred', label='Possible Factories (Not Producing)')
                plt.plot(capitalLatLon[1,30], capitalLatLon[0,30], marker='o', markersize=12, color=[97/255., 218/255., 97/255.], markeredgewidth=1.5, markeredgecolor='k', label = 'Intl Shipment Port')
                plt.plot(capitalLatLon[1,30], capitalLatLon[0,30], marker='o', markersize=7, color='darkred', label = 'Intl Shipment Port (No Shipments)')
        
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
                    if cName=='Somaliland':
                        cName='Somalia'
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
                        x=factoryPctOne[g,c]
                        y=y1+(y2-y1)/(cmax-cmin)*(x-cmin)
                        icmap=min(255,int(round(y,1)))
                        icmap=max(0,int(round(icmap,1)))
                        ax.add_geometries(country.geometry, ccrs.PlateCarree(), edgecolor='black', facecolor=[cmapArray[icmap,0],cmapArray[icmap,1],cmapArray[icmap,2]],label=cName)
        
                        if x!=0:
                            size = 10*(1+factoryPctOne[g,c]/cmax)
                            plt.plot(capitalLatLon[1,c], capitalLatLon[0,c], marker='*', markersize=size, color=[cmapArray[icmap,0],cmapArray[icmap,1],cmapArray[icmap,2]], markeredgewidth=1.5, markeredgecolor='k')
                            factoryNumOne+=1
                        if x==0:
                            plt.plot(capitalLatLon[1,c], capitalLatLon[0,c], marker='*', markersize=7, color='darkred')
        
                
                for icoast in range(24,len(countrycosted)):
                    x=factoryPctOne[g,icoast]
                    y=y1+(y2-y1)/(cmax-cmin)*(x-cmin)
                    icmap=min(255,int(round(y,1)))
                    icmap=max(0,int(round(icmap,1)))
                    if x!=0:
                        size = 10*(1+factoryPctOne[g,icoast]/cmax)
                        plt.plot(capitalLatLon[1,icoast], capitalLatLon[0,icoast], marker='o', markersize=size, color=[cmapArray[icmap,0],cmapArray[icmap,1],cmapArray[icmap,2]], markeredgewidth=1.5, markeredgecolor='k')
                        IntlNumOne+=1
                    if x==0:
                        plt.plot(capitalLatLon[1,icoast], capitalLatLon[0,icoast], marker='o', markersize=7, color='darkred')
        
                local = str(int(np.round(100*np.sum(factoryPctOne[g,:24])/np.sum(factoryPctOne[g,:]),0)))
                intl = str(np.round(100*np.sum(factoryPctOne[g,24:])/np.sum(factoryPctOne[g,:]),0))
                #costOne = str(int(round(costOne/1000000.,0)))
        
                plt.title('Production of '+SMtitles[g]+' Treatment by Factory and Port\n' + LTitles[L] + TTitles[T])
                plt.legend(loc = 'lower left')
                #plt.text(-15,-10,str(factoryNumOne)+' Factories Open\n'+str(IntlNumOne)+' Ports Open\n'+local+'% Produced Locally\nTotal Cost = $'+costOne+' Million', bbox=dict(fc="none", boxstyle="round"), size = 10)
                plt.text(-15,-10,str(factoryNumOne)+' Factories Open\n'+str(IntlNumOne)+' Ports Open\n'+local+'% Produced Locally', bbox=dict(fc="none", boxstyle="round"), size = 10)
                
                plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/geographical/'+SMtitles[g]+'_export_map.pdf')

        ###########################
        # Skeleton
        ###########################
        if MakePlots:
            shapename = 'admin_0_countries'
            countries_shp = shpreader.natural_earth(resolution='110m',
                category='cultural', name=shapename)
            
            plt.clf()
            ax = plt.axes([0.05,0.05,0.8,0.85],projection=ccrs.PlateCarree())
            ax.set_extent([-19, 53, -37, 39], ccrs.PlateCarree())
            ax.coastlines()

            plt.plot(capitalLatLon[1,8], capitalLatLon[0,8], marker='*', markersize=9, color='orangered', label='Possible Factories')
            plt.plot(capitalLatLon[1,30], capitalLatLon[0,30], marker='o', markersize=8, color='dodgerblue', label = 'Possible Intl Shipment Ports')
            plt.plot(SScapitalLatLon[1,9],SScapitalLatLon[0,9], marker='^', markersize=8, color='mediumpurple', label = 'Recieves Treatment')

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
                if cName=='Somaliland':
                    cName='Somalia'
                if np.amax(cName==subsaharancountry)==0:
                    continue
                if np.amax(cName==countrycosted)!=0:
                    c=np.where(cName==countrycosted)[0][0]
                    lon1=capitalLatLon[1,c]
                    lat1=capitalLatLon[0,c]
                    for iSS in range(len(subsaharancountry)):
                        lat2=SScapitalLatLon[0,iSS]
                        lon2=SScapitalLatLon[1,iSS]
                        dist=np.sqrt((lat2-lat1)**2+(lon2-lon1)**2)
                        if dist<15:
                            plt.plot([lon1,lon2] , [lat1,lat2], color='gray', linestyle='--', linewidth = 0.5, transform=ccrs.PlateCarree() )
            
            for icoast in range(24,len(countrycosted)):
                lon1=capitalLatLon[1,icoast]
                lat1=capitalLatLon[0,icoast]
                for iSS in range(len(subsaharancountry)):
                    lat2=SScapitalLatLon[0,iSS]
                    lon2=SScapitalLatLon[1,iSS]
                    dist=np.sqrt((lat2-lat1)**2+(lon2-lon1)**2)
                    if dist<17:
                        plt.plot([lon1,lon2] , [lat1,lat2], color='gray', linestyle='--', linewidth = 0.5, transform=ccrs.PlateCarree() )

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
                if cName=='Somaliland':
                    cName='Somalia'
                if np.amax(cName==subsaharancountry)==0:
                    continue
                if np.amax(cName==countrycosted)==0:
                    ax.add_geometries(country.geometry, ccrs.PlateCarree(), edgecolor='black',
                        facecolor='lightgray')
                    c=np.where(cName==subsaharancountry)[0][0]
                    plt.plot(SScapitalLatLon[1,c],SScapitalLatLon[0,c], marker='^', markersize=8, color='mediumpurple')
                else:
                    c=np.where(cName==countrycosted)[0][0]
                    ax.add_geometries(country.geometry, ccrs.PlateCarree(), edgecolor='black', facecolor='lightgreen',label=cName)
                    plt.plot(capitalLatLon[1,c], capitalLatLon[0,c], marker='*', markersize=9, color='orangered')

            for icoast in range(24,len(countrycosted)):
                plt.plot(capitalLatLon[1,icoast], capitalLatLon[0,icoast], marker='o', markersize=8, color='dodgerblue')

            plt.title('Supply Chain Optimization\nPossible Factory and Port Locations')
            plt.legend(loc = 'lower left')
            plt.text(-15,-10,'24 Possible Factories\n9 Possible Ports', bbox=dict(fc="none", boxstyle="round"), size = 10)
            plt.savefig(wdfigs+'skeleton_map.pdf')

        ###########################
        # Supply Zones
        ###########################
        if MakePlots:
            ruftitles=['rutf','rusf']
            for g in range(2):
                productarray = np.load(wddata+'results/example/'+optiLevel[L]+'/RN'+ruftitles[g]+'array.npy')
                Rcountrycosted1=np.load(wddata+'results/example/'+optiLevel[L]+'/RNcountry.npy')
                Rsubsaharancountry1=np.load(wddata+'results/example/'+optiLevel[L]+'/RNsubsaharancountry.npy')
                Rcountrycosted=[]
                for i in range(len(Rcountrycosted1)):
                    country=Rcountrycosted1[i]
                    if country[:2]=='I_':
                        countrytmp=country[2:].replace('_',' ')
                        Rcountrycosted.append('I_'+countrytmp)
                    else:
                        countrytmp=country.replace('_',' ')
                        Rcountrycosted.append(countrytmp)
                Rcountrycosted=np.array(Rcountrycosted)
                Rsubsaharancountry=[]
                for i in range(len(Rsubsaharancountry1)):
                    country=Rsubsaharancountry1[i]
                    if country[:2]=='I_':
                        countrytmp=country[2:].replace('_',' ')
                        Rsubsaharancountry.append('I_'+countrytmp)
                    else:
                        countrytmp=country.replace('_',' ')
                        Rsubsaharancountry.append(countrytmp)
                Rsubsaharancountry=np.array(Rsubsaharancountry)

                Rsubsaharancountry[Rsubsaharancountry=='Congo']='DRC'
                Rsubsaharancountry[Rsubsaharancountry=='Congo (Republic of the)']='Congo'
                Rsubsaharancountry[Rsubsaharancountry=="Cote d'Ivoire"]='Ivory Coast'
                Rcountrycosted[Rcountrycosted=='Congo']='DRC'
                Rcountrycosted[Rcountrycosted=='Congo (Republic of the)']='Congo'
                Rcountrycosted[Rcountrycosted=="I_Cote d'Ivoire"]='I_Ivory Coast'
                Rcountrycosted[Rcountrycosted=="Cote d'Ivoire"]='Ivory Coast'
                
                shapename = 'admin_0_countries'
                countries_shp = shpreader.natural_earth(resolution='110m', category='cultural', name=shapename)
                colors = [(240,59,32),(252,146,114),(254,178,76),(255,237,160),(35,132,67),(49,163,84),(229,245,224),(0,0,139),(49,130,189),(158,202,225),(136,86,167),(158,188,218)]                
                        # colors = [(240,59,32),(252,146,114),(254,178,76),(255,237,160),(49,163,84),(161,217,155),(229,245,224),(49,130,189),(158,202,225),(136,86,167),(158,188,218)]
                # colors = [(128,0,0),(170,110,40),(128,128,0),(0,128,128),(0,0,128),(0,0,128),(0,0,0),(230,25,75),(245,130,48),(255,225,25),(210,245,60),(60,180,75),(70,240,240),(0,130,200),(145,30,180),(240,50,230),(128,128,128),(250,190,190),(255,215,180),(255,250,200),(170,255,195),(230,190,255),(255,255,255)]
                colors=colors[:len(Rcountrycosted)+1]
                my_cmap = make_cmap(colors,bit=True)
                
                plt.clf()
                cmapArray=my_cmap(np.arange(256))
                cmin=0
                cmax=len(Rcountrycosted)
                y1=0
                y2=255
                
                fig = plt.figure(figsize=(10, 8))
                MinMaxArray=np.ones(shape=(3,2))
                subPlot1 = plt.axes([0.61, 0.07, 0.2, 0.8])
                MinMaxArray[0,0]=cmin
                MinMaxArray[1,0]=cmax
                plt.imshow(MinMaxArray,cmap=my_cmap)
                plt.colorbar()
                
                ax = plt.axes([0.05,0.05,0.8,0.85],projection=ccrs.PlateCarree())
                ax.set_extent([-19, 53, -37, 39], ccrs.PlateCarree())
                ax.coastlines()

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
                    if cName=='Somaliland':
                        cName='Somalia'
                    if np.amax(cName==Rsubsaharancountry)==0:
                        continue
                    else:
                        poz=np.where(cName==Rsubsaharancountry)[0][0]
                        c=np.where(productarray[:,poz]==np.amax(productarray[:,poz]))[0][0]
                        y=y1+(y2-y1)/(cmax-cmin)*(c-cmin)
                        icmap=min(255,int(round(y,1)))
                        icmap=max(0,int(round(icmap,1)))
                        ax.add_geometries(country.geometry, ccrs.PlateCarree(), edgecolor='black', facecolor=[cmapArray[icmap,0],cmapArray[icmap,1],cmapArray[icmap,2]],label=cName)

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
                    if cName=='Somaliland':
                        cName='Somalia'
                    if np.amax(cName==Rsubsaharancountry)==0:
                        continue
                    else:
                        poz=np.where(cName==Rsubsaharancountry)[0][0]
                        c=np.where(productarray[:,poz]==np.amax(productarray[:,poz]))[0][0]
                        width=0.2+(productarray[c,poz]/np.amax(productarray))

                        p = np.where(cName==subsaharancountry)[0][0]
                        lat2=SScapitalLatLon[0,p]
                        lon2=SScapitalLatLon[1,p]

                        supplier=Rcountrycosted[c]
                        p2=np.where(supplier==countrycosted)[0][0]
                        lat1=capitalLatLon[0,p2]
                        lon1=capitalLatLon[1,p2]

                        dlat=lat2-lat1
                        dlon=lon2-lon1
                        if dlat!=0:
                            plt.arrow(lon1, lat1, dlon, dlat, color='k', linestyle='-', width=width, head_width=2.5*width, head_length=width*2, length_includes_head=True, transform=ccrs.PlateCarree() )

                for f in range(len(Rcountrycosted)):
                    factory=Rcountrycosted[f]
                    if factory[:2]!='I_':
                        continue

                    y=y1+(y2-y1)/(cmax-cmin)*(f-cmin)
                    icmap=min(255,int(round(y,1)))
                    icmap=max(0,int(round(icmap,1)))

                    p=np.where(factory==countrycosted)[0][0]

                    if factoryPctOne[g,p]!=0:
                        size = 10*(1+factoryPctOne[g,p]/np.amax(factoryPctOne[g,:]))
                        plt.plot(capitalLatLon[1,p], capitalLatLon[0,p], marker='o', markersize=size, markerfacecolor=[cmapArray[icmap,0],cmapArray[icmap,1],cmapArray[icmap,2]], markeredgewidth=1.5, markeredgecolor='k', label=factory[2:]+' Port',linestyle = 'None')
                        IntlNumOne+=1

                for f in range(len(Rcountrycosted)):
                    factory=Rcountrycosted[f]
                    if factory[:2]=='I_':
                        continue
                    y=y1+(y2-y1)/(cmax-cmin)*(f-cmin)
                    icmap=min(255,int(round(y,1)))
                    icmap=max(0,int(round(icmap,1)))

                    p=np.where(factory==countrycosted)[0][0]

                    if factoryPctOne[g,p]!=0:
                        size = 10*(1+factoryPctOne[g,p]/np.amax(factoryPctOne[g,:]))
                        plt.plot(capitalLatLon[1,p], capitalLatLon[0,p], marker='*', markersize=size, markerfacecolor=[cmapArray[icmap,0],cmapArray[icmap,1],cmapArray[icmap,2]], markeredgewidth=1.5, markeredgecolor='k',label=factory,linestyle = 'None')
                        factoryNumOne+=1
                    #if x==0:
                    #    plt.plot(capitalLatLon[1,p], capitalLatLon[0,p], marker='*', markersize=7, color='darkred')
                

                local = str(int(np.round(100*np.sum(factoryPctOne[g,:24])/np.sum(factoryPctOne[g,:]),0)))
                intl = str(np.round(100*np.sum(factoryPctOne[g,24:])/np.sum(factoryPctOne[g,:]),0))
                # costOne = str(int(round(costOne/1000000.,0)))

                plt.legend(loc = 'lower left')

                plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/exports_by_country/'+Ltitles[L]+'_'+factory+'_exports.pdf')

                plt.title(SMtitles[g]+' Treatment: Primary Supplier by Country\n' + LTitles[L] + TTitles[T])
                plt.legend(loc = 'lower left',ncol=1,numpoints=1)
                plt.text(10,25,str(factoryNumOne)+' Factories Open\n'+str(IntlNumOne)+' Ports Open\n'+local+'% Produced Locally', bbox=dict(fc="none", boxstyle="round"), size = 10)
                plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/geographical/'+SMtitles[g]+'supplyzone_map.pdf')

        ###########################
        # By factory import/export
        ###########################
        if MakePlots:
            for g in range(2):
                productarray = np.load(wddata+'results/example/'+optiLevel[L]+'/RN'+ruftitles[g]+'array.npy')
                Rcountrycosted1=np.load(wddata+'results/example/'+optiLevel[L]+'/RNcountry.npy')
                Rsubsaharancountry1=np.load(wddata+'results/example/'+optiLevel[L]+'/Rsubsaharancountry.npy')
                Rcountrycosted=[]
                for i in range(len(Rcountrycosted1)):
                    country=Rcountrycosted1[i]
                    if country[:2]=='I_':
                        countrytmp=country[2:].replace('_',' ')
                        Rcountrycosted.append('I_'+countrytmp)
                    else:
                        countrytmp=country.replace('_',' ')
                        Rcountrycosted.append(countrytmp)
                Rcountrycosted=np.array(Rcountrycosted)
                Rsubsaharancountry=[]
                for i in range(len(Rsubsaharancountry1)):
                    country=Rsubsaharancountry1[i]
                    if country[:2]=='I_':
                        countrytmp=country[2:].replace('_',' ')
                        Rsubsaharancountry.append('I_'+countrytmp)
                    else:
                        countrytmp=country.replace('_',' ')
                        Rsubsaharancountry.append(countrytmp)
                Rsubsaharancountry=np.array(Rsubsaharancountry)

                Rsubsaharancountry[Rsubsaharancountry=='Congo']='DRC'
                Rsubsaharancountry[Rsubsaharancountry=='Congo (Republic of the)']='Congo'
                Rsubsaharancountry[Rsubsaharancountry=="Cote d'Ivoire"]='Ivory Coast'
                Rcountrycosted[Rcountrycosted=='Congo']='DRC'
                Rcountrycosted[Rcountrycosted=='Congo (Republic of the)']='Congo'
                Rcountrycosted[Rcountrycosted=="I_Cote d'Ivoire"]='I_Ivory Coast'
                Rcountrycosted[Rcountrycosted=="Cote d'Ivoire"]='Ivory Coast'
                
                colors = [(255,255,255), (203,208,255), (160,169,255), (121,133,255), (79, 95, 255), (43, 62, 255), (0, 23, 255)]
                my_cmap = make_cmap(colors,bit=True)
                shapename = 'admin_0_countries'
                countries_shp = shpreader.natural_earth(resolution='110m',
                    category='cultural', name=shapename)
                
                for f in range(len(productarray)):

                    factory = Rcountrycosted[f]

                    plt.clf()
                    cmapArray=my_cmap(np.arange(256))
                    cmin=0
                    cmax=np.amax(productarray[f,:]) #*0.9
                    if cmax==0:
                        continue
                    y1=0
                    y2=255
                    
                    fig = plt.figure(figsize=(10, 8))
                    MinMaxArray=np.ones(shape=(3,2))
                    subPlot1 = plt.axes([0.61, 0.07, 0.2, 0.8])
                    MinMaxArray[0,0]=cmin
                    MinMaxArray[1,0]=cmax
                    plt.imshow(MinMaxArray,cmap=my_cmap)
                    plt.colorbar()
                    plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/'+Ltitles[L]+'_'+factory+'_exports.pdf')
                    
                    ax = plt.axes([0.05,0.05,0.8,0.85],projection=ccrs.PlateCarree())
                    ax.set_extent([-19, 53, -37, 39], ccrs.PlateCarree())
                    ax.coastlines()
            
                    plt.plot(-16.1, -34.7, marker='*', markersize=9, color='limegreen', label='Factory',linestyle = 'None')
                    plt.plot(-16.1, -34.7, marker='o', markersize=8, color='limegreen', label = 'Intl Shipment Port',linestyle = 'None')
                    plt.plot(-16.1, -34.7, marker='^', markersize=8, color='mediumpurple', label = 'Recieves Treatment',linestyle = 'None')
                    impCountries=[]
                    impPct=[]
                        
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
                        impc=np.where(cName==subsaharancountry)[0][0]
                        x=productarray[f,impc]
                        y=y1+(y2-y1)/(cmax-cmin)*(x-cmin)
                        icmap=min(255,int(round(y,1)))
                        icmap=max(0,int(round(icmap,1)))
                        ax.add_geometries(country.geometry, ccrs.PlateCarree(), edgecolor='black', facecolor=[cmapArray[icmap,0],cmapArray[icmap,1],cmapArray[icmap,2]],label=cName)
            
                        if x!=0:
                            impCountries.append(cName)
                            impPct.append(x)
                            size = 10*(1+x/cmax)
                            plt.plot(SScapitalLatLon[1,impc], SScapitalLatLon[0,impc], marker='^', markersize=8, color='mediumpurple')
                            facc=np.where(factory==countrycosted)[0][0]
                            if factory[:2]=='I_':
                                plt.plot(capitalLatLon[1,facc], capitalLatLon[0,facc], marker='o', markersize=12, color='limegreen')
                            else:
                                plt.plot(capitalLatLon[1,facc], capitalLatLon[0,facc], marker='*', markersize=13, color='limegreen')

                            factoryNumOne+=1
            
                    
                    #for icoast in range(24,len(countrycosted)):
                    #    x=factoryPctOne[g,icoast]
                    #    if x!=0:
                    #        size = 10*(1+factoryPctOne[g,icoast]/cmax)
                    #        plt.plot(capitalLatLon[1,icoast], capitalLatLon[0,icoast], marker='o', markersize=size, color='g')
                    #        IntlNumOne+=1
                    #    if x==0:
                    #        plt.plot(capitalLatLon[1,icoast], capitalLatLon[0,icoast], marker='o', markersize=7, color='darkred')
            
                    totalshipments = np.sum(productarray[f])
                    impPct=np.array(impPct)
                    impPct=100*impPct/totalshipments
                    order=np.argsort(impPct)
                    impPct=impPct[order][::-1]
                    impCountries=np.array(impCountries)[order][::-1]
                    totalshipments = str(int(round(np.sum(productarray[f])/1000000.)))
                    #local = str(int(np.round(100*np.sum(factoryPctOne[g,:24])/np.sum(factoryPctOne[g,:]),0)))
                    #intl = str(np.round(100*np.sum(factoryPctOne[g,24:])/np.sum(factoryPctOne[g,:]),0))
                    #costOne = str(int(round(costOne/1000000.,0)))
            
                    plt.title('Exports of RUTF for '+factory+', Packets \n' + LTitles[L] + TTitles[T])
                    if factory[:2]=='I_':
                        plt.title(SMtitles[g]+' Treatment Supplied by '+factory[2:]+' Port\n' + LTitles[L] + TTitles[T])
                        plt.text(-15,-8,factory[2:]+' Port\n'+totalshipments+' Million Packets Procured', size = 10)
                        for r in range(len(impCountries)):
                            plt.text(-13,-10.4-1.7*r,'- '+impCountries[r]+', '+str(int(round(impPct[r])))+'%',size=10)
                    else:
                        plt.title(SMtitles[g]+' Treatment Supplied by '+factory+' Factory\n' + LTitles[L] + TTitles[T])
                        plt.text(-15,-8,factory+' Factory\n'+totalshipments+' Million Packets Procured', size = 10)
                        for r in range(len(impCountries)):
                            plt.text(-13,-10.4-1.7*r,'- '+impCountries[r]+', '+str(int(round(impPct[r])))+'%',size=10)

                    plt.legend(loc = 'lower left')
                    plt.savefig(wdfigs+Ltitles[L]+'/'+Ttitles[T]+'/exports_by_country/'+Ltitles[L]+'_'+factory+'_exports.pdf')
    
## cost barchart ##
fig = plt.figure(figsize=(6, 5))
plt.clf()
x=np.array([1,2,3])
ydata = (np.array([costOneAll[0,1],costOneAll[1,1],costOneAll[3,1]])/1e9)[::-1]
colors=['g','b','r'][::-1]
plt.bar(x,ydata,color=colors,tick_label=['Current','Local Optimized','All Optimized'])
plt.ylabel('Total Cost of Procurement for 1 Year (Billion USD)')
plt.title('Total Modeled Cost')
plt.savefig(wdfigs+'cost_optimized_vs_current_barchart.pdf')

## % local barchart ##
fig = plt.figure(figsize=(6, 5))
plt.clf()
x=np.array([1,2,3])
pctLocalOneAll1 = np.mean(pctLocalOneAll[:,1,:],axis=1)
ydata = (np.array([pctLocalOneAll1[0],pctLocalOneAll1[1],pctLocalOneAll1[3]])*100)[::-1]
colors=['g','b','r'][::-1]
plt.bar(x,ydata,color=colors,tick_label=['Current','Local Optimized','All Optimized'])
plt.ylabel('% Treatment Produced Locally')
plt.title('Percent Produced Locally')
plt.savefig(wdfigs+'pctLocal_optimized_vs_current_barchart.pdf')

## factoryNum barchart ##
fig = plt.figure(figsize=(6, 5))
plt.clf()
x1=np.array([0.79,1.79,2.79])
x2=np.array([1.21,2.21,3.21])

bar_width=0.4
ydata = (np.array([factoryNumOneAll[0,1],factoryNumOneAll[1,1],factoryNumOneAll[3,1]]))[::-1]
colors=['g','b','r'][::-1]
plt.bar(x1,ydata,color=colors,width=bar_width) #,tick_label=['Factories','Factories','Factories'])

ydata = (np.array([portNumOneAll[0,1],portNumOneAll[1,1],portNumOneAll[3,1]]))[::-1]
plt.bar(x2,ydata,color=colors,width=bar_width) #,tick_label=['Current','Local Optimized','All Optimized'])

plt.yticks([0,2,4,6,8,10,12,14,16,18,20])
plt.xticks([1,2,3],['Current','Local Optimized','All Optimized'])
plt.text(0.62,0.4,'Factories',size=8)
plt.text(1.62,0.4,'Factories',size=8)
plt.text(2.62,0.4,'Factories',size=8)
plt.text(1.12,0.4,'Ports',size=8)
plt.text(2.12,0.4,'Ports',size=8)
plt.text(3.12,0.4,'Ports',size=8)
plt.ylabel('Number of Factories or Ports')
plt.title('Number of Factories and Ports')
plt.savefig(wdfigs+'factoryNum_optimized_vs_current_barchart.pdf')

fig = plt.figure(figsize=(7, 4))
LTitles = ['All Optimized','Local Optimized','Optimized Intl','Current']

T=1
cost1=cost/1e9
for V in range(len(loopvar)):
    x = np.arange(mins[V],maxs[V],factor[V])
    x = x*100
    plt.clf()
    #plt.plot(x,np.ma.compressed(np.ma.masked_array(cost[2,T,V,:],Mask[2,T,V,:])),'b*-',label=LTitles[2])
    plt.plot(x,np.ma.compressed(np.ma.masked_array(cost1[3,T,V,:],Mask[3,T,V,:])),'r*-',label=LTitles[3])
    plt.plot(x,np.ma.compressed(np.ma.masked_array(cost1[1,T,V,:],Mask[1,T,V,:])),'b*-',label=LTitles[1])
    plt.plot(x,np.ma.compressed(np.ma.masked_array(cost1[0,T,V,:],Mask[0,T,V,:])),'g*-',label=LTitles[0])

    plt.title('Effect of '+VTitles[V]+' on Total Cost')
    plt.xlabel(VTitles[V]+' Cost, % of Today')
    plt.ylabel('Total Procurement Cost for One Year (Billion USD)')
    plt.ylim([0.5,1.35])
    plt.grid(True)
    plt.legend(loc='lower right')
    plt.savefig(wdfigs+'totalCost_vs_'+Vtitles[V]+'_'+Ttitles[T]+'.pdf')

for T in range(len(1)):
    for V in range(len(loopvar)):
        x = np.arange(mins[V],maxs[V],factor[V])
        x = x*100
        plt.clf()
        plt.plot(x,np.ma.compressed(np.ma.masked_array(factoryNum[0,T,V,:],Mask[0,T,V,:])),'g*-',label=LTitles[0])
        plt.plot(x,np.ma.compressed(np.ma.masked_array(factoryNum[1,T,V,:],Mask[1,T,V,:])),'c*-',label=LTitles[1])
        #plt.plot(x,np.ma.compressed(np.ma.masked_array(factoryNum[2,T,V,:],Mask[2,T,V,:])),'b*-',label=LTitles[2])
        plt.plot(x,np.ma.compressed(np.ma.masked_array(factoryNum[3,T,V,:],Mask[3,T,V,:])),'r*-',label=LTitles[3])

        plt.title('Effect of '+VTitles[V]+' on Number of Factories'+TTitles[T])
        plt.xlabel(VTitles[V]+' Cost, % of Today')
        plt.ylabel('Number of Factories')
        #plt.ylim([0,2e9])
        plt.grid(True)
        plt.legend()
        plt.savefig(wdfigs+'factoryNum_vs_'+Vtitles[V]+'_'+Ttitles[T]+'.pdf')

for T in range(len(1)):
    for V in range(len(loopvar)):
        x = np.arange(mins[V],maxs[V],factor[V])
        x = x*100
        plt.clf()
        plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[0,T,V,:,0],Mask[0,T,V,:])),'g*-',label=LTitles[0])
        plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[1,T,V,:,0],Mask[1,T,V,:])),'c*-',label=LTitles[1])
        #plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[2,T,V,:,0],Mask[2,T,V,:])),'b*-',label=LTitles[2])
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

for T in range(len(1)):
    for V in range(len(loopvar)):
        x = np.arange(mins[V],maxs[V],factor[V])
        x = x*100
        plt.clf()
        plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[0,T,V,:,1],Mask[0,T,V,:])),'g*-',label=LTitles[0])
        plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[1,T,V,:,1],Mask[1,T,V,:])),'c*-',label=LTitles[1])
        #plt.plot(x,100*np.ma.compressed(np.ma.masked_array(pctLocal[2,T,V,:,1],Mask[2,T,V,:])),'b*-',label=LTitles[2])
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
exit()
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

#Rcountry2=[]
#for i in range(len(Rcountry)):
#    country=Rcountry[i]
#    if country[:2]=='I_':
#        countrytmp=country[:2].replace('_',' ')
#        Rcountry2.append('I_'+countrytmp)
#    else:
#        countrytmp=country.replace('_',' ')
#        Rcountry2.append(country)
#Rcountry2=np.array(Rcountry2)

plt.clf()
fig = plt.figure(figsize=(13, 8))
plt.imshow(rusfarray,cmap=cm.terrain_r,vmax=2.5e8)
plt.colorbar()
plt.xticks(np.arange(43), subsaharancountry, rotation='vertical')
plt.yticks(np.arange(len(Rcountry2)), Rcountry2)
plt.title('Imports and Exports (Packets, Ordered by Lon)')
plt.savefig(wdfigs+'impexp_array.pdf')

