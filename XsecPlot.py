import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import csv

sigmaB = []
sigmaNLO = []
sigmaRes = []
sigmaExp = []
sigmaResHS = []
sigmaExpHS = []
sigmaResHSConvert = []
sigmaExpHSConvert = []
sigmaNLL=[]
RatioNLO=[]
RatioExp =[]
RatioRes=[]
RatioExpHS =[]
RatioResHS=[]
RatioExpHSConvert =[]
RatioResHSConvert=[]
RatioNLL=[]
s = []

with open('NLOXsec.txt','r') as NLO:
    plotNLO = csv.reader(NLO, delimiter=',')
    for row in plotNLO:
        sigmaNLO.append(float(row[0]))    
        s.append(float(row[1]))
        
with open('LOXsec.txt','r') as LO:
    plotLO = csv.reader(LO, delimiter=',')
    for row in plotLO:
        sigmaB.append(float(row[0]))
        
with open('ResumXsec.txt','r') as Resum:
    plotRes = csv.reader(Resum, delimiter=',')
    for row in plotRes:
        sigmaRes.append(float(row[0]))
        
with open('ExpandedXsec.txt','r') as Exp:
    plotExp = csv.reader(Exp, delimiter=',')
    for row in plotExp:
        sigmaExp.append(float(row[0]))

# with open('ResumXsecHS.txt','r') as ResumHS:
#     plotResHS = csv.reader(ResumHS, delimiter=',')
#     for row in plotResHS:
#         sigmaResHS.append(float(row[0]))
        
# with open('ExpandedXsecHS.txt','r') as ExpHS:
#     plotExpHS = csv.reader(ExpHS, delimiter=',')
#     for row in plotExpHS:
#         sigmaExpHS.append(float(row[0]))

# with open('ResumXsecHSConvert.txt','r') as ResumHSConvert:
#     plotResHSConvert = csv.reader(ResumHSConvert, delimiter=',')
#     for row in plotResHSConvert:
#         sigmaResHSConvert.append(float(row[0]))
        
# with open('ExpandedXsecHSConvert.txt','r') as ExpHSConvert:
#     plotExpHSConvert= csv.reader(ExpHSConvert, delimiter=',')
#     for row in plotExpHSConvert:
#         sigmaExpHSConvert.append(float(row[0]))

for k in range(len(sigmaB)):
    sigmaNLL.append(sigmaNLO[k]+sigmaRes[k]-sigmaExp[k])
    RatioExp.append(sigmaExp[k]/sigmaB[k])
    RatioNLO.append(sigmaNLO[k]/sigmaB[k])
    RatioRes.append(sigmaRes[k]/sigmaB[k])
    # RatioExpHS.append(sigmaExpHS[k]/sigmaB[k])
    # RatioResHS.append(sigmaResHS[k]/sigmaB[k])
    # RatioExpHSConvert.append(sigmaExpHSConvert[k]/sigmaB[k])
    # RatioResHSConvert.append(sigmaResHSConvert[k]/sigmaB[k])
    RatioNLL.append(sigmaNLL[k]/sigmaB[k])

LO1, =plt.plot(s,sigmaB, label='LO cross section')
NLO1, =plt.plot(s,sigmaNLO, label='NLO cross section')
Exp1, = plt.plot(s,sigmaExp, label=r'$NLL_{NLO}$ cross section')
Res1, = plt.plot(s,sigmaRes, label='NLL cross section')
plt.xlabel(r'$\sqrt{s}$ in GeV')
plt.ylabel(r'$\sigma$ in pb')
plt.title(r'LO and NLO cross section distribution for q q > z process')
plt.legend(handles=[LO1, NLO1, Exp1])#plt.legend(handles=[LO1, NLO1, Exp1, Res1])
plt.savefig('XSec.png')
plt.legend([])
plt.close()

LO2, =plt.plot(s,sigmaB, label='LO cross section')
NLO2, =plt.plot(s,sigmaNLO, label='NLO cross section')
Exp2, = plt.plot(s,sigmaExp, label=r'$NLL_{NLO}$ cross section')
Res2, = plt.plot(s,sigmaRes, label='NLL cross section')
plt.xscale("log")
plt.xlabel(r'$\sqrt{s}$ in GeV')
plt.ylabel(r'$\sigma$ in pb')
plt.title(r'LO and NLO cross section distribution for q q > z process')
plt.legend(handles=[LO2, NLO2, Exp2, Res2])
#plt.legend(handles=[LO2, NLO2, Exp2])
plt.savefig('XSecLog.png')
plt.legend([])
plt.close()

LO2, =plt.plot(s,sigmaB, label='LO cross section')
NLO2, =plt.plot(s,sigmaNLO, label='NLO cross section')
Exp2, = plt.plot(s,sigmaExp, label=r'$NLL_{NLO}$  cross section')
Res2, = plt.plot(s,sigmaRes, label='NLL cross section')
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r'$\sqrt{s}$ in GeV')
plt.ylabel(r'$\sigma$ in pb')
plt.title(r'LO and NLO cross section distribution for q q > z process')
plt.legend(handles=[LO2, NLO2, Exp2, Res2])
#plt.legend(handles=[LO2, NLO2, Exp2])
plt.savefig('XSecLogLog.png')
plt.legend([])
plt.close()

NLOR, =plt.plot(s,RatioNLO, label='NLO / LO ')
ExpR, = plt.plot(s,RatioExp, label=r'$NLL_{NLO}$ / LO')
ResR, = plt.plot(s,RatioRes, label='NLL / LO')
NLL, = plt.plot(s,RatioNLL, label=r'NLL+NLO-$NLL_{NLO}$ / LO')
plt.xscale("log")
plt.xlabel(r'$\sqrt{s}$ in GeV')
#plt.ylabel(r'$\dfrac{\sigma}{\sigma_0}$')
plt.title(r'Distribution normalized to LO for q q > z process')
plt.legend(handles=[NLOR,ExpR,ResR,NLL])
#plt.legend(handles=[NLOR, ExpR])
plt.savefig('XSecRatiosLog.png')
plt.legend([])
plt.close()


# NLOR2, =plt.plot(s,RatioNLO, label='NLO / LO ')
# ResR2, = plt.plot(s,RatioRes, label='Resummed / LO')
# #ResRHS, = plt.plot(s,RatioResHS, label='$Resummed_{HS}$ / LO')
# ResRHSConvert, = plt.plot(s,RatioResHSConvert, label='$Resummed_{Converted}$ / LO')
# plt.xscale("log")
# plt.xlabel(r'$\sqrt{s}$ in GeV')
# #plt.ylabel(r'$\dfrac{\sigma}{\sigma_0}$')
# plt.title(r'Distribution normalized to LO for q q > z process')
# plt.legend(handles=[NLOR2,ResR2,ResRHSConvert])
# #plt.legend(handles=[NLOR, ExpR])
# plt.savefig('XSecRatiosLogComparisonHSRes.png')
# plt.legend([])
# plt.close()

# NLOR3, =plt.plot(s,RatioNLO, label='NLO / LO ')
# ExpR3, = plt.plot(s,RatioExp, label=r'$NLL_{NLO}$ / LO')
# #ExpRHS, = plt.plot(s,RatioExpHS, label=r'$NLL_{NLO, \, HS}$ / LO')
# ExpRHSConvert, = plt.plot(s,RatioExpHSConvert, label=r'$NLL_{NLO, \, Converted}$ / LO')
# plt.xscale("log")
# plt.xlabel(r'$\sqrt{s}$ in GeV')
# #plt.ylabel(r'$\dfrac{\sigma}{\sigma_0}$')
# plt.title(r'Distribution normalized to LO for q q > z process')
# plt.legend(handles=[NLOR3,ExpR3,ExpRHSConvert])
# #plt.legend(handles=[NLOR, ExpR])
# plt.savefig('XSecRatiosLogComparisonHSExp.png')
# plt.legend([])
# plt.close()
