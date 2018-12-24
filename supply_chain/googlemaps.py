import csv
import numpy as np

listof50k = open('C:/Users/garyk/Documents/python_code/riskAssessmentFromPovertyEstimations/data/population/MSH_50K_TX_3/MSH_50K_TX.csv','r')

citynames=["" for x in range(629)]
countrynames=["" for x in range(629)]
previouscity=1
firstline=True
index = 0
for line in listof50k:
    if firstline:
        firstline = False
        continue
    tmp=line.split(',')
    marketname=tmp[0]
    countrycode=tmp[2]
    if(np.amax([marketname==citynames[i] for i in range(len(citynames))])==0):
        citynames[index]=marketname
        countrynames[index]=countrycode
        index=index+1
    previouscity=tmp[0]