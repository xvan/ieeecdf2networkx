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
dfY=pd.DataFrame(Y,columns=pd.Index(nodes,name="vCoord"),index=pd.Index(nodes,name="vCoord"))
#Uso pandas dataframe para indexar con los labes de los buses

for node in nodes:
    for edge in Gn.edges([node],data=True):
        dfY.loc[node,edge[1]] -= np.complex(edge[2]["attr_dict"]["r"],edge[2]["attr_dict"]["x"])**-1

for tnode in dnodes:
    dfY.loc[tnode[0],tnode[0]]=np.complex(tnode[1]["attr_dict"]["shunt_g"],tnode[1]["attr_dict"]["shunt_b"])-df2.loc[tnode[0],:].sum()


dfYmod=dfY.apply(np.abs)
dfYang=dfY.apply(np.angle)

dfV=pd.DataFrame(np.stack([np.ones(dim),np.zeros(dim)],1),index=pd.Index(nodes,name="bus"),columns=pd.Index(["V","d"],name="vCoord"))

def P(n,dfV):
    return (dfV.loc[n,"V"] * dfV.loc[:,"V"] * dfYmod.loc[n,:] * (dfV.loc[:,"d"] + dfYang.loc[n,:] - dfV.loc[n,"d"]).apply(np.cos) ).sum()


def Q(n,V,d):
    return - (dfV.loc[n,"V"] * dfV.loc[:,"V"] * dfYmod.loc[n,:] * (dfV.loc[:,"d"] + dfYang.loc[n,:] - dfV.loc[n,"d"]).apply(np.sin) ).sum()

def dPdd(n,m,dfV):
    if n == m:
        return ( dfV.loc[n,"V"] * dfV.loc[:,"V"] * dfYmod.loc[n,:] * ( dfV.loc[:,"d"] + dfYang.loc[n,:] - dfV.loc[n,"d"]).apply(np.sin) ).drop(n).sum()
    else:
        return dfV.loc[n,"V"] * dfV.loc[m,"V"] * dfYmod.loc[n,m] * np.sin(dfV.loc[m,"d"] + dfYang.loc[n,m] - dfV.loc[n,"d"])


######## Hasta aca llegué


def dQdd(n,m,V,d):
    if n == m:
       return np.sum(np.delete( V[n] * np.multiply( np.multiply(V,Ymod[n,:]) , np.cos(d + Yang[n,:] - d[n]) ),n))
    else:
       return -V[n] * V[m] * Ymod[n,m]* np.cos(d[m] + Yang[n,m] - d[n])

def dPdV(n,m,V,d):
    if n == m:
        return np.sum ( np.multiply( np.multiply(V,Ymod[n,:]) , np.cos(d + Yang[n,:] - d[n]) )) + V[n]*Ymod[n,n]*np.cos(Yang[n,n])
    else:
        return V[n] * Ymod[n,m] * np.cos(d[m] + Yang[n,m] - d[n])

def dQdV(n,m,V,d):
    if n == m:
       return - np.sum( np.multiply( np.multiply(V,Ymod[n,:]) , np.sin(d + Yang[n,:] - d[n]) )) - V[n]*Ymod[n,n]*np.sin(Yang[n,n])
    else:
       return - V[n] * Ymod[n,m] * np.sin(d[m] + Yang[n,m] - d[n])


#agrupo buses por tipo, (PV,PQ,swing)
bus_group=[[] for _ in range(4)]
for node in dnodes:
    bus_group[ node[1]["attr_dict"]["bus_type"] ].append(node)

if len(bus_group[3])!=1:
    raise ValueError('Only one swing bus admitted')


dfV.loc["V"].as_matrix()
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

def Jac(V,d):
    dfJ1=pd.DataFrame(J_matrix,columns=J_colindex,index=J_rowindex)
    for i in dfJ1.index:
        for j in dfJ1.columns:
            if i[0]=="P" and j[0]=="d":
                dfJ1.loc[i,j] = dPdd(i[1],j[1],V,d)
            else if i[0]=="P" and j[0]=="V":
                dfJ1.loc[i,j] = dPdV(i[1],j[1],V,d)
            else if i[0]=="Q" and j[0]=="d":
                dfJ1.loc[i,j] = dQdd(i[1],j[1],V,d)
            else if i[0]=="Q" and j[0]=="V":
                dfJ1.loc[i,j] = dQdd(i[1],j[1],V,d)
