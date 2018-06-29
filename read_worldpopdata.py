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
from PIL import Image
from libtiff import TIFF

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

Poverty=True

if Poverty:
	countriesT=['Malawi','Nigeria','Uganda','Tanzania']
	countries=[]
	for c in countriesT:
		countries.append(c.lower())
	countryAbb=['mwi','nga','uga','tza']
	years=['11','10','10','10']
	
	i=-1
	for country in countries:
		i+=1
	
		tif125 = TIFF.open(wddata+'poverty/'+country+'/'+countryAbb[i]+years[i]+'povcons125.tif',mode='r')
		tif200 = TIFF.open(wddata+'poverty/'+country+'/'+countryAbb[i]+years[i]+'povcons200.tif',mode='r')
		tifuncert200 = TIFF.open(wddata+'poverty/'+country+'/'+countryAbb[i]+years[i]+'povcons200uncert.tif',mode='r')
		tifuncert125 = TIFF.open(wddata+'poverty/'+country+'/'+countryAbb[i]+years[i]+'povcons125uncert.tif',mode='r')
	
		pov125 = tif125.read_image()*100
		pov200 = tif200.read_image()*100
		uncert125 = tifuncert125.read_image()*100
		uncert200 = tifuncert200.read_image()*100
		
		imageMask=np.ones(shape=(pov125.shape))
		for x in range(pov125.shape[0]):
			for y in range(pov125.shape[1]):
				if pov125[x,y]>0:
					imageMask[x,y]=0
		
		pov125=np.ma.masked_array(pov125,imageMask)
		pov125=pov125[~imageMask.all(axis=1)]
		pov125=pov125[:,~imageMask.all(axis=0)]
		
		pov200=np.ma.masked_array(pov200,imageMask)
		pov200=pov200[~imageMask.all(axis=1)]
		pov200=pov200[:,~imageMask.all(axis=0)]
		
		uncert125=np.ma.masked_array(uncert125,imageMask)
		uncert125=uncert125[~imageMask.all(axis=1)]
		uncert125=uncert125[:,~imageMask.all(axis=0)]
		
		uncert200=np.ma.masked_array(uncert200,imageMask)
		uncert200=uncert200[~imageMask.all(axis=1)]
		uncert200=uncert200[:,~imageMask.all(axis=0)]
		
		plt.clf()
		plt.imshow(pov125,cmap=cm.hot_r,vmin=0,vmax=100)
		plt.yticks([])
		plt.xticks([])
		plt.title('Percent living under 1.25/day '+countriesT[i]+' 20'+years[i])
		plt.colorbar()
		plt.savefig(wdfigs+country+'_poverty_125',dpi=700)
		
		plt.clf()
		plt.imshow(pov200,cmap=cm.hot_r,vmin=0,vmax=100)
		plt.yticks([])
		plt.xticks([])
		plt.title('Percent living under 2USD/day '+countriesT[i]+' 20'+years[i])
		plt.colorbar()
		plt.savefig(wdfigs+country+'_poverty_200',dpi=700)
		
		plt.clf()
		plt.imshow(uncert200,cmap=cm.winter_r,vmin=0,vmax=100)
		plt.yticks([])
		plt.xticks([])
		plt.title(countriesT[i]+' 95 Percent Confidence Interval for Pov2USD')
		plt.colorbar()
		plt.savefig(wdfigs+country+'_uncertainty200',dpi=700)
		
		plt.clf()
		plt.imshow(uncert125,cmap=cm.winter_r,vmin=0,vmax=100)
		plt.yticks([])
		plt.xticks([])
		plt.title(countriesT[i]+' 95 Percent Confidence Interval for Pov1.25USD')
		plt.colorbar()
		plt.savefig(wdfigs+country+'_uncertainty125',dpi=700)

######################################################################
# Age Structures
######################################################################
for y in range(0,66,5):
	for g in ('M','F'):
		if y<10:
			print 'A0'+str(y)+'0'+str(y+4)+g
			# load tiff
			vars()['tifA0'+str(y)+'0'+str(y+4)+g]=TIFF.open(wddata+'africa_1km_age_structures/AFR_PPP_A0'+str(y)+'0'+str(y+4)+'_'+g+'_2015_adj_v5.tif',mode='r')

			# convert tiff to array
			vars()['A0'+str(y)+'0'+str(y+4)+g]=vars()['tifA0'+str(y)+'0'+str(y+4)+g].read_image() 
			# Mask
			imageMask=vars()['A0'+str(y)+'0'+str(y+4)+g]<0
			vars()['A0'+str(y)+'0'+str(y+4)+g]=np.ma.masked_array(vars()['A0'+str(y)+'0'+str(y+4)+g],imageMask)
			vars()['A0'+str(y)+'0'+str(y+4)+g]=vars()['A0'+str(y)+'0'+str(y+4)+g][~imageMask.all(axis=1)]
			vars()['A0'+str(y)+'0'+str(y+4)+g]=vars()['A0'+str(y)+'0'+str(y+4)+g][:,~imageMask.all(axis=0)]
			
		elif y==65:
			print 'A'+str(y)+'PL'+g
			# load tiff
			vars()['tifA'+str(y)+'PL'+g]=TIFF.open(wddata+'africa_1km_age_structures/AFR_PPP_A'+str(y)+'PL'+'_'+g+'_2015_adj_v5.tif',mode='r')

			# convert tiff to array
			vars()['A'+str(y)+'PL'+g]=vars()['tifA'+str(y)+'PL'+g].read_image()

			imageMask=vars()['A'+str(y)+'PL'+g]<0
			vars()['A'+str(y)+'PL'+g]=np.ma.masked_array(vars()['A'+str(y)+'PL'+g],imageMask)
			vars()['A'+str(y)+'PL'+g]=vars()['A'+str(y)+'PL'+g][~imageMask.all(axis=1)]
			vars()['A'+str(y)+'PL'+g]=vars()['A'+str(y)+'PL'+g][:,~imageMask.all(axis=0)]
		else:
k			print 'A'+str(y)+str(y+4)+g
			# load tiff
			vars()['tifA'+str(y)+str(y+4)+g]=TIFF.open(wddata+'africa_1km_age_structures/AFR_PPP_A'+str(y)+str(y+4)+'_'+g+'_2015_adj_v5.tif',mode='r')

			# convert tiff to array
			vars()['A'+str(y)+str(y+4)+g]=vars()['tifA'+str(y)+str(y+4)+g].read_image()

			imageMask=vars()['A'+str(y)+str(y+4)+g]<0
			vars()['A'+str(y)+str(y+4)+g]=np.ma.masked_array(vars()['A'+str(y)+str(y+4)+g],imageMask)
			vars()['A'+str(y)+str(y+4)+g]=vars()['A'+str(y)+str(y+4)+g][~imageMask.all(axis=1)]
			vars()['A'+str(y)+str(y+4)+g]=vars()['A'+str(y)+str(y+4)+g][:,~imageMask.all(axis=0)]

ageStructures=np.zeros(shape=(8,2))
gender=['M','F']
xc,yc=3738,1064
i=-1
for y in range(0,36,5):
	i+=1
	for g in range(2):
		if y<10:
			ageStructures[i,g]=vars()['A0'+str(y)+'0'+str(y+4)+gender[g]][xc,yc]
		else:
			ageStructures[i,g]=vars()['A'+str(y)+str(y+4)+gender[g]][xc,yc]

X = np.arange(0,36,5)


plt.clf()
plt.barh(X, ageStructures[:,1],height=4, color = 'm',alpha=.7,linewidth=2,edgecolor='k')
plt.barh(X, -ageStructures[:,0],height=4, color = 'b',alpha=.7,linewidth=2,edgecolor='k')
plt.axvline(x=0,color='k',linewidth=4)
plt.title('Age Pyramid in a 1km Grid Cell')
plt.ylabel('Ages')
plt.xlabel('Number')
plt.savefig(wdfigs+'sample_age_structure',dpi=700)

exit()


tifA0004F=TIFF.open(wddata+'africa_1km_age_structures/AFR_PPP_A0004_F_2015_adj_v5.tif',mode='r')
A0004F=tifA0004F.read_image()
tifA0004M=TIFF.open(wddata+'africa_1km_age_structures/AFR_PPP_A0004_M_2015_adj_v5.tif',mode='r')
A0004M=tifA0004F.read_image()
tifWOCBA=TIFF.open(wddata+'africa_1km_age_structures/AFR_PPP_WOCBA_2015_adj_v5.tif',mode='r')
WOCBA=tifWOCBA.read_image()

imageMask=A0004F<0

A0004F=np.ma.masked_array(A0004F,imageMask)
A0004M=np.ma.masked_array(A0004M,imageMask)
WOCBA=np.ma.masked_array(WOCBA,imageMask)

A0004=A0004F+A0004M

plt.clf()
plt.imshow(A0004,cmap=cm.jet,vmin=0,vmax=40)
plt.yticks([])
plt.xticks([])
plt.title('Number of Children Age 0-4 2015')
plt.colorbar()
plt.savefig(wdfigs+'A0004',dpi=700)

plt.clf()
plt.imshow(WOCBA,cmap=cm.jet,vmin=0,vmax=50)
plt.yticks([])
plt.xticks([])
plt.title('Women of Child Bearing Age 2015')
plt.colorbar()
plt.savefig(wdfigs+'pop_WOCBA',dpi=700)

A04perWOCBA=A0004/WOCBA
A04perWOCBA=A04perWOCBA[~imageMask.all(axis=1)]
A04perWOCBA=A04perWOCBA[:,~imageMask.all(axis=0)]

plt.clf()
plt.imshow(A04perWOCBA,cmap=cm.jet)
plt.yticks([])
plt.xticks([])
plt.title('Children 0-4 per Women of Child Bearing Age 2015')
plt.colorbar()
plt.savefig(wdfigs+'pop_A04perWOCBA',dpi=700)
