from pulp import *
import math
import json
import numpy as np
import re

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

scaleaverage=np.zeros(shape=43)
f=open(wddata+'travel_time/averagetkmcost.csv','r')
i=-1
for line in f:
    i+=1
    scaleaverage[i]=line
transportcostArray=np.zeros(shape=(43,43))
f=open(wddata+'travel_time/capitaldistanceArray.csv','r')
i=-1
for line in f:
    i+=1
    tmp=line.split(',')
    for j in range(len(tmp)):
        avg=(scaleaverage[i]+scaleaverage[j])/2
        transportcostArray[i,j]=float(tmp[j])*avg

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
for z in range(16):
    #cost dabber RUTF ########################################################################
    rutfcostarray=np.zeros(shape=(24,43))
    for i in range(len(subsaharancountry)):
        if(indexedrutf[i]!=0):
            for j in range(len(indexedSAM)):
 			# sums ingredient and transport cost, converts to $/100g delivered
                rutfcostarray[int(convertarray[i]),j]=indexedrutf[i]+100/1000000*transportcostArray[i,j]+z/10*100/1000000*importExportCosts[i,j]/15
    
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
                            mamcostarray[int(convertarray[i]),j]=indexedscp[i]+200/1000000*transportcostArray[i,j]+z/10*200/1000000*importExportCosts[i,j]/15
    
                            print subsaharancountry[i] +' SC+'
                    else:
                        for j in range(len(indexedSAM)):
                            mamcostarray[int(convertarray[i]),j]=indexedrusf[i]+100/1000000*transportcostArray[i,j]+z/10*100/1000000*importExportCosts[i,j]/15
    
    
    np.savetxt(wddata + "foodstuffs/mamcostarray.csv", mamcostarray, delimiter=",")
    
    wddata='C:/Users/garyk/Documents/code/riskAssessmentFromPovertyEstimations/supply_chain/data/'
    
    facility = ['Angola',
    'Benin',
    'Burkina Faso',
    'Burundi',
    'Cameroon',
    'Chad',
    'Congo (Republic of the)',
    "Cote d'Ivoire",
    'Ethiopia',
    'Gambia',
    'Ghana',
    'Kenya',
    'Malawi',
    'Mali',
    'Mozambique',
    'Nigeria',
    'Rwanda',
    'Senegal',
    'South Africa',
    'Tanzania',
    'Togo',
    'Uganda',
    'Zambia',
    'Zimbabwe']
    
    location = ['Angola',
    'Benin',
    'Botswana',
    'Burkina Faso',
    'Burundi',
    'Cameroon',
    'Central African Republic',
    'Chad',
    'Congo (Republic of the)',
    'Congo',
    "Cote d'Ivoire",
    'Djibouti',
    'Equatorial Guinea',
    'Eritrea',
    'Ethiopia',
    'Gabon',
    'Gambia',
    'Ghana',
    'Guinea',
    'Guinea-Bissau',
    'Kenya',
    'Lesotho',
    'Liberia',
    'Malawi',
    'Mali',
    'Mauritania',
    'Mozambique',
    'Namibia',
    'Niger',
    'Nigeria',
    'Rwanda',
    'Senegal',
    'Sierra Leone',
    'Somalia',
    'South Africa',
    'Sudan',
    'South Sudan',
    'Swaziland',
    'Tanzania',
    'Togo',
    'Uganda',
    'Zambia',
    'Zimbabwe']
    
    #fixed costs for facility
    startupcost= 100000.0
    # startupcost = dict(zip(facility, [1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000]))
    upgradecost= 20000.0
    
    #machinery costs for production, small machinery
    costM1 = 6000.0
    #machinery costs for production, large machinery
    costM2 = 10000.0
    #fixed capacity per piece of machinery, small machinery
    Capacity1 = 780000
    #fixed capacity per piece of machinery, large machinery
    Capacity2 = (780000*2.0)
    #demand by location
    rutfdemandarray = np.genfromtxt(wddata+'optiarrays/SAMdemand.csv', delimiter=',')
    DemandRUTF = dict(zip(location, rutfdemandarray))
    
    mamdemandarray = np.genfromtxt(wddata+'optiarrays/MAMdemand.csv', delimiter=',')
    DemandMAM = dict(zip(location, mamdemandarray))
    
    # cost to location
    with open(wddata+'optiarrays/rutfdictionary.json', 'r') as fp:
        data = json.load(fp)
    CostRUTF = data
    
    with open(wddata+'optiarrays/mamdictionary.json', 'r') as fp:
        data = json.load(fp)
    CostMAM = data
    
    QuantityRUTF = LpVariable.dicts ('Supply of RUTF %s%s',(facility, location),cat = 'Continuous',lowBound = 0,upBound = None)
    
    QuantityMAM = LpVariable.dicts ('Supply of MAM Treatment %s%s',(facility, location),cat = 'Continuous',lowBound = 0,upBound = None)
    
    #number of small machines
    Machine1 = LpVariable.dicts('Machine 1 %s', facility,
                        lowBound = 0,
                        cat='Integer')
    #number of large machines
    Machine2 = LpVariable.dicts('Machine 2 %s', facility,
                        lowBound = 0,
                        cat='Integer')
    # factory open? Y/N
    Open = LpVariable.dicts('Factory Status %s', facility,
                        cat='Binary')
    #factory size upgrades
    Factorysize = LpVariable.dicts('Factory Size %s', facility,
                        lowBound = 0,
                        cat='Integer')
    
    prob = LpProblem('Fixed Charge', LpMinimize)
    tmp1 = sum(costM1 * Machine1[i] for i in facility)
    tmp2 = sum(costM2 * Machine2[i] for i in facility)
    tmp3 = sum(startupcost * Open[i] for i in facility)
    tmp4 = sum(upgradecost * Factorysize[i] for i in facility)
    tmp5 = sum(sum(CostRUTF[i][j] * QuantityRUTF[i][j] for j in location) for i in facility)
    tmp6 = sum(sum(CostMAM[i][j] * QuantityMAM[i][j] for j in location) for i in facility)
    prob+=tmp1+tmp2+tmp3+tmp4+tmp5+tmp6
    
    #must be less than small machinery
    for i in facility:
        prob += sum(QuantityRUTF[i][j]+QuantityMAM[i][j] for j in location) <= (Capacity1*Machine1[i] + Capacity2*Machine2[i])
    #must be less than maximum amount a factory can produce (this ensures fixed cost is used)
    for i in facility:
        prob += sum(QuantityRUTF[i][j]+QuantityMAM[i][j] for j in location) <= 1000000000000*Open[i]
    for i in facility:
        prob += sum(Machine2[i]+Machine1[i]) <= 10*Factorysize[i]
    for j in location:
        prob += sum(QuantityRUTF[i][j] for i in facility) >= DemandRUTF[j]
    for j in location:
        prob += sum(QuantityMAM[i][j] for i in facility) >= DemandMAM[j]
        
    print(prob)
    
    prob.solve()
    print("Status:")
    print(LpStatus[prob.status])
    print("Objective:")
    print(value(prob.objective))
    for v in prob.variables():
        if(v.varValue>0):
            print (v.name, "=", v.varValue)
    
    varsdict = {}
    for v in prob.variables():
        if(v.varValue>0):
            varsdict[v.name] = v.varValue
    
    with open(wddata+'optiarrays/temp_results.json', 'w') as fp:
        json.dump(varsdict, fp, sort_keys=True)
    
    cost=0
    factorynum=0
    averageshipments=0
    countries=[]
    cost=value(prob.objective)
    sizes=[]
    for v in prob.variables():
        if (v.varValue>0):
            if (v.name[0:10]=="Factory_St"):
                factorynum+=1
                # print v.name[15:]
                countries.append(v.name[15:])
            if (v.name[0:9]=="Machine_2"):
                sizes.append(value(v))
            if (v.name[0:10]=="Supply_of_"):
                averageshipments+=1
    averageshipments=averageshipments/(factorynum*2)
    
    sizes=sizes*780000
    
    f = open(wddata + 'bordercostfactor'+str(z)+'.csv','w')
    f.write('cost'+','+str(cost)+'\n')
    f.write('num_factories'+','+str(factorynum)+'\n')
    for i in range(len(countries)):
        f.write(str(countries[i])+','+str(sizes[i])+'\n')
    f.close()
    
            # if v.name[0,5]="
        
    # resultsarray=np.zeros(shape=(86,2))
    # i=-1
    # for v in varsdict():
    #     if(v.varValue>0 and v.name[:4] == 'Supp'):
    #         i+=1
    #         # resultsarray[i,0]=v.name
    #         resultsarray[i,1]=v.value
    #         # if(v.name[10:14]=="RUTF"):
    #         #     resultsarray[i,0]= 'RUTF'
    #         # else:
    #         #     resultsarray[i,0]= 'RUSF'
    #         # tmp = re.findall('[A-Z][^A-Z]*', v.name[15:])
    #         # resultsarray[i,1]=float(tmp[0])
    #         # resultsarray[i,2]=v.varValue
    # 
    # np.savetxt(wddata + "optiarrays/results/results.csv", resultsarray, delimiter=",")