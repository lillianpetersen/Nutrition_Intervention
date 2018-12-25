from pulp import *
import math
import json
import numpy as np
import re

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
    prob += sum(QuantityRUTF[i][j]+QuantityMAM[i][j] for j in location) <= 1000000000*Open[i]
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