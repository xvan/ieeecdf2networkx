from ieeecdf2networkx import IEEECDFParser
import pandas as pd
import numpy as np
import itertools

x=IEEECDFParser('Testcdf.txt')
Gn=x.generate_networkx_graph()

nodes=Gn.node
dnodes=Gn.nodes(data=True)

dim=len(nodes)
Y=np.zeros((dim,dim),dtype=np.complex)
dfY=pd.DataFrame(Y,columns=pd.Index(nodes,name="bus"),index=pd.Index(nodes,name="bus"))
#Uso pandas dataframe para indexar con los labes de los buses

for node in nodes:
    for edge in Gn.edges([node],data=True):
        dfY.loc[node,edge[1]] -= np.complex(edge[2]["attr_dict"]["r"],edge[2]["attr_dict"]["x"])**-1

for tnode in dnodes:
    dfY.loc[tnode[0],tnode[0]]=np.complex(tnode[1]["attr_dict"]["shunt_g"],tnode[1]["attr_dict"]["shunt_b"])-dfY.loc[tnode[0],:].sum()


dfYmod=dfY.apply(np.abs)
dfYang=dfY.apply(np.angle)

def P(n,dfV):
    return   (dfV.loc[("V",n)] * dfV.loc["V"] * dfYmod.loc[n,:] * (dfV.loc["d"] + dfYang.loc[n,:] - dfV.loc[("d",n)]).apply(np.cos) ).sum()

def Q(n,dfV):
    return - (dfV.loc[("V",n)] * dfV.loc["V"] * dfYmod.loc[n,:] * (dfV.loc["d"] + dfYang.loc[n,:] - dfV.loc[("d",n)]).apply(np.sin) ).sum()

def dPdd(n,m,dfV):
    if n == m:
        return ( dfV.loc[("V",n)] * dfV.loc["V"] * dfYmod.loc[n,:] * ( dfV.loc["d"] + dfYang.loc[n,:] - dfV.loc[("d",n)]).apply(np.sin) ).drop(n).sum()
    else:
        return - dfV.loc[("V",n)] * dfV.loc[("V",m)] * dfYmod.loc[n,m] * np.sin(dfV.loc[("d",m)] + dfYang.loc[n,m] - dfV.loc[("d",n)])

def dQdd(n,m,dfV):
    if n == m:
        return ( dfV.loc[("V",n)] * dfV.loc["V"] * dfYmod.loc[n,:] * ( dfV.loc["d"] + dfYang.loc[n,:] - dfV.loc[("d",n)]).apply(np.cos) ).drop(n).sum()
    else:
        return - dfV.loc[("V",n)] * dfV.loc[("V",m)] * dfYmod.loc[n,m] * np.cos(dfV.loc[("d",m)] + dfYang.loc[n,m] - dfV.loc[("d",n)])

######## Hasta aca llegue

def dPdV(n,m,dfV):
    if n == m:
        return ( dfV.loc["V"] * dfYmod.loc[n,:] * ( dfV.loc["d"] + dfYang.loc[n,:] - dfV.loc[("d",n)]).apply(np.cos) ).sum() + dfV.loc[("V",n)]* dfYmod.loc[n,n] * np.cos(dfYang.loc[n,n])
    else:
        return dfV.loc[("V",n)] * dfYmod.loc[n,m] * np.cos( dfV.loc[("d",m)] + dfYang.loc[n,m] - dfV.loc[("d",n)])

def dQdV(n,m,dfV):
    if n == m:
       return  - ( dfV.loc["V"] * dfYmod.loc[n,:] * ( dfV.loc["d"] + dfYang.loc[n,:] - dfV.loc[("d",n)]).apply(np.sin) ).sum() - dfV.loc[("V",n)]* dfYmod.loc[n,n] * np.sin(dfYang.loc[n,n])
    else:
       return - dfV.loc[("V",n)] * dfYmod.loc[n,m] * np.sin( dfV.loc[("d",m)] + dfYang.loc[n,m] - dfV.loc[("d",n)])





#TODO: Cambia con condiciones de contorno?
dfV=pd.Series(np.concatenate([np.ones(dim),np.zeros(dim)],0),index=pd.MultiIndex.from_product([ ["V","d"], nodes ],names=["vCoord","bus"]))
for tnode in dnodes:
    if tnode[1]["attr_dict"]["bus_type"] in [2,3]:
        dfV[("V",tnode[0])] = tnode[1]["attr_dict"]["desired_volts"]

#agrupo buses por tipo, (PV,PQ,swing)
bus_group=[[] for _ in range(4)]
for tnode in dnodes:
    bus_group[ tnode[1]["attr_dict"]["bus_type"] ].append(tnode)

if len(bus_group[3])!=1:
    raise ValueError('Only one swing bus admitted')


#dfV.loc["V"].as_matrix()
#dimension
#
# |DP|   |J1 J2| |Dd|
# |  | = |     | |  |
# |DQ|   |J3 J4| |DV|
#

#TODO: Agregar condiciones de contorno
PQ_buses=bus_group[0]+bus_group[1]
PV_buses=bus_group[2]



J_index={"cols": { "d": dict(PQ_buses+PV_buses) , "V": dict(PQ_buses) } ,
         "rows": { "P": dict(PQ_buses+PV_buses) , "Q": dict(PQ_buses) }
        }

J_colindex=pd.MultiIndex.from_product([ ["d"],J_index["cols"]["d"].keys() ],names=["vCoord","bus"]).union(pd.MultiIndex.from_product([ ["V"],J_index["cols"]["V"].keys() ],names=["vCoord","bus"]))
J_rowindex=pd.MultiIndex.from_product([ ["P"],J_index["rows"]["P"].keys() ],names=["sCoord","bus"]).union(pd.MultiIndex.from_product([ ["Q"],J_index["rows"]["Q"].keys() ],names=["sCoord","bus"]))
J_matrix=np.matrix(np.zeros(( len(J_rowindex),len(J_colindex)) ))

def Jac(dfV):
    dfJ1=pd.DataFrame(J_matrix,columns=J_colindex,index=J_rowindex)
    for i in dfJ1.index:
        for j in dfJ1.columns:
            if   i[0]=="P" and j[0]=="d":
                dfJ1.loc[i,j] = dPdd(i[1],j[1],dfV)
            elif i[0]=="P" and j[0]=="V":
                dfJ1.loc[i,j] = dPdV(i[1],j[1],dfV)
            elif i[0]=="Q" and j[0]=="d":
                dfJ1.loc[i,j] = dQdd(i[1],j[1],dfV)
            elif i[0]=="Q" and j[0]=="V":
                dfJ1.loc[i,j] = dQdd(i[1],j[1],dfV)
    return dfJ1




SSched = pd.Series(np.zeros(len(J_rowindex)),index=J_rowindex)
for i in SSched.index:
    if i[0]=="P":
      SSched.loc[i]=dnodes[i[1]]["attr_dict"]["gen_p"]-dnodes[i[1]]["attr_dict"]["load_p"]
    elif i[0]=="Q":
      SSched.loc[i]=dnodes[i[1]]["attr_dict"]["gen_q"]-dnodes[i[1]]["attr_dict"]["load_q"]

#Calculo de residuos
def SRes(dfV):
    dfDP = SSched.copy()
    for i in dfDP.index:
        if i[0]=="P":
            dfDP.loc[i]+=P(i[1],dfV)
        elif i[0]=="Q":
            dfDP.loc[i]+=P(i[1],dfV)
    return dfDP




dfRes=SRes(dfV)
while np.abs(dfRes.max()) > 1E-4:
    dfJ=Jac(dfV) #Jacobiano
    dfIJ=pd.DataFrame(np.linalg.inv(dfJ.values),dfJ.columns,dfJ.index) #Inversa
    dfNewV = dfIJ.dot(dfRes) # Producto por residuos
    dfV.update(dfNewV) #Actualizo tensiones
    dfRes=SRes(dfV)  #Recalculo residuos
    print dfV
