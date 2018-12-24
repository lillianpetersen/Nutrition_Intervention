import sys
import os
import matplotlib.pyplot as plt
#import descarteslabs as dl
import numpy as np
import math
from sys import exit
import sklearn
import time
from collections import Counter
from sklearn.preprocessing import StandardScaler
from operator import and_
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from datetime import datetime

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

def get_second_highest(a):
	hi = mid = 0
	for x in a:
		if x > hi:
			mid = hi
			hi = x
		elif x < hi and x > mid:
			lo = mid
			mid = x
	return mid

def get_third_highest(a):
	hi = mid = third = 0
	for x in a:
		if x > hi:
			third = mid
			mid = hi
			hi = x
		elif x < hi and x > mid:
			third = mid
			mid = x
		elif x < mid and x > third:
			third = x
	return third

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

colors = [(.4,0,.6), (0,0,.7), (0,.6,1), (.9,.9,1), (1,.8,.8), (1,1,0), (.8,1,.5), (.1,.7,.1), (.1,.3,.1)]
my_cmap = make_cmap(colors)
my_cmap_r=make_cmap(colors[::-1])
###############################################

try:
    wddata='/Users/lilllianpetersen/iiasa/data/'
    wddata2='/Users/lilllianpetersen/data/'
    wdfigs='/Users/lilllianpetersen/iiasa/figs/'
    wdvars='/Users/lilllianpetersen/iiasa/saved_vars/'
except:
    wddata='C:/Users/garyk/Documents/python_code/riskAssessmentFromPovertyEstimations/data/'
    wdfigs='C:/Users/garyk/Documents/python_code/riskAssessmentFromPovertyEstimations/figs/'
    wdvars='C:/Users/garyk/Documents/python_code/riskAssessmentFromPovertyEstimations/vars/'


nyears=38
ncountries=47

countryList=[]
indexAll=[]
cropAll=[]
pop=np.zeros(shape=(ncountries,nyears))
countryCounter=-1

for icountry in range(47):
	countryNum=str(icountry+1)
	
	cropDataDir=wddata2+'africa_crop/'
	
	f=open(wddata2+'africa_latlons.csv')
	for line in f:
		tmp=line.split(',')
		if tmp[0]==countryNum:
			country=tmp[1]
			countryl=country.lower()
			countryl=countryl.replace(' ','_')
			countryl=countryl.replace('-','_')
			harvestMonth=float(tmp[7])
			seasons=tmp[8][:-1].split('/')
			twoSeasons=seasons[0]
			if twoSeasons!='no':
				corrMonth1=int(seasons[1])
				corrMonth2=int(seasons[2])
				corrMonth=corrMonth1
			break

	####################################################################

	print '\n',country,countryNum

	files = [filename for filename in os.listdir(cropDataDir) if filename.startswith(countryl+'-')]
	crop=[]
	for n in files:
		tmp=n.split('-')
		croptmp=tmp[1].title()
		crop.append(tmp[1].title())
	
	cropYield=np.zeros(shape=(len(files),38))
	for cp in range(len(files)):
		fCropData=open(wddata2+'africa_crop/'+countryl+'-'+crop[cp]+'-production.csv')
		j=-1
		k=-1
		for line in fCropData:
			j+=1
			if j==0:
				continue
			tmp=line.split(',')
			year=int(tmp[0].replace('"',''))
			if year<1980 or year>=2000+nyears:
				continue
			k+=1
			cropYield[cp,k]=float(tmp[1].replace('"',''))

	for cp in range(len(crop)):
		if crop[cp]=='Centrifugal_Sugar':
			crop[cp]='Sugar'
		elif crop[cp]=='Milled_Rice':
			crop[cp]='Rice'
		elif crop[cp]=='Green_Coffee' or crop[cp]=='Green_Coffe':
			crop[cp]='Coffee'
	
	if np.amax(cropYield)==np.amin(cropYield):
		continue
	print crop

	#fpop = open(wddata+'country_indices/country_population.csv')
	#for line in fpop:
	#	line=line.replace('"','')
	#	tmp = line.split(',')
	#	if tmp[0]==country:
	#		i=3
	#		y=-1
	#		for year in range(1960,2018):
	#			i+=1
	#			if year<2000 or year>2000+nyears-1:
	#				continue
	#			y+=1
	#			pop[icountry,y]=tmp[i]
	#		break

	#plt.clf()
	#plt.plot(range(2000,2016),pop[0,:])
	#plt.title('Ethiopia Populaiton')
	#plt.grid(True)
	#plt.savefig(wdfigs+'ethiopia_pop.pdf')

	cropYieldD=np.zeros(shape=(cropYield.shape))
	x = np.arange(1980,1980+nyears)
	for cp in range(len(crop)):
		ydata=np.ma.compressed(cropYield[cp,:])

		ydataAvg=np.mean(ydata)
		A,B,C=np.polyfit(x,ydata,2)
		yfit2=A*(x**2)+B*x+C
		m,b=np.polyfit(x,ydata,1)
		yfit1=m*x+b

		plt.clf()
		plt.plot(x,cropYield[cp,:],'*-')
		plt.plot(x,yfit1,'-')
		plt.plot(x,yfit2,'g-')
		plt.title('Ethiopia Corn Production')
		plt.ylabel('Corn Yield')
		plt.grid(True)
		plt.savefig(wdfigs+'ethipia_yield.pdf')
	
		cropYieldD[cp,:]=cropYield[cp,:]-yfit2
	
	### Find Crop Yield Anomaly ###
	cropYieldAnom=np.zeros(shape=(cropYieldD.shape))
	meanYield=np.zeros(shape=(len(crop)))
	stdYield=np.zeros(shape=(len(crop)))
	for cp in range(len(crop)):
		meanYield[cp]=np.mean(cropYieldD[cp,:])
		stdYield[cp]=np.std(cropYieldD[cp,:])
		cropYieldAnom[cp,:]=cropYieldD[cp,:]/stdYield[cp]
	
	if np.ma.amax(cropYieldAnom)==0. and np.amin(cropYieldAnom)==0:
		cropYieldAnom[:]=1

	plt.clf()
	plt.plot(x,cropYieldAnom[0,:],'*-')
	plt.plot([1900,2026],[0,0],'k')
	plt.xlim([1979.9,2017.5])
	plt.title('Ethiopia Corn Anomaly')
	plt.ylabel('Corn Yield, Std Dev from Avg')
	plt.grid(True)
	plt.savefig(wdfigs+'yield_anom.pdf')
	
	countryList.append(country)
	countryCounter+=1
	exit()

