from ieeecdf2networkx import IEEECDFParser
import pandas as pd
import numpy as np
import itertools

x=IEEECDFParser('ieee14cdf.txt')
Gn=x.generate_networkx_graph()

nodes=Gn.node
dnodes=Gn.nodes(data=True)

dim=len(nodes)
Y=np.zeros((dim,dim),dtype=np.complex)
df2=pd.DataFrame(Y,columns=nodes,index=nodes)
#Uso pandas dataframe para indexar con los labes de los buses

for node in nodes:
    for edge in Gn.edges([node],data=True):
        df2.loc[node,edge[1]] -= np.complex(edge[2]["attr_dict"]["r"],edge[2]["attr_dict"]["x"])**-1

for tnode in dnodes:
    df2.loc[tnode[0],tnode[0]]=np.complex(tnode[1]["attr_dict"]["shunt_g"],tnode[1]["attr_dict"]["shunt_b"])-df2.loc[tnode[0],:].sum()


Ymod=df2.apply(np.abs)
Yang=df2.apply(np.angle)


def P(n,V,d):
    return np.sum ( V[n] * np.multiply( np.multiply(V,Ymod[n,:]) , np.cos(d + Yang[n,:] - d[n]) ))

def Q(n,V,d):
    return - np.sum ( V[n] * np.multiply( np.multiply(V,Ymod[n,:]) , np.sin(d + Yang[n,:] - d[n]) ))

def dPdd(n,m,V,d):
    if n == m:
       return np.sum( V[n] * np.multiply( np.multiply(V,np.delete(Ymod[n,:],n) , np.sin(d + np.delete(Yang[n,:],n)) - d[n]) ))
    else:
       return -V[n] * V[m] * Ymod[n,m]* np.sin(d[m] + Yang[n,m] - d[n])

def dQdd(n,m,V,d):
    if n == m:
       return np.sum( V[n] * np.multiply( np.multiply(V,np.delete(Ymod[n,:],n) , np.cos(d + np.delete(Yang[n,:],n)) - d[n]) ))
    else:
       return -V[n] * V[m] * Ymod[n,m]* np.cos(d[m] + Yang[n,m] - d[n])

def dPdV(n,m,V,d):
    if n == m:
        return np.sum ( np.multiply( np.multiply(V,Ymod[n,:]) , np.cos(d + Yang[n,:] - d[n]) )) + V[n]*Ymod[n,n]*np.cos(Yang[n,n])
    else:
        return V[n] * Ymod[n,m] * np.cos(d[m] + Yang[n,m] - d[n])

def dQdV(n,m,V,d):
    if n == m:
       return - np.sum( np.multiply( np.multiply(V,Ymod[n,:]) , np.sin(d + Yang[n,:] - d[n]) )) - V[n]*Ymod[n,n]*np.cos(Yang[n,n])
    else:
       return - V[n] * Ymod[n,m] * np.sin(d[m] + Yang[n,m] - d[n])



V=np.ones(dim)
d=np.zeros(dim)

#agrupo buses por tipo, (PV,PQ,swing)
bus_group=[[] for _ in range(4)]
for node in dnodes:
    bus_group[ node[1]["attr_dict"]["bus_type"] ].append(node)

if len(bus_group[3])!=1:
    raise ValueError('Only one swing bus admitted')
#dimension
#
# |DP|   |J1 J2| |Dd|
# |  | = |     | |  |
# |DQ|   |J3 J4| |DV|
#

#TODO: Agregar problemas de contorno
PQ_buses=bus_group[0]+bus_group[1]
PV_buses=bus_group[2]

dj1=len(PQ_buses)+len(PV_buses)
J1=np.matrix(np.zeros((dj1,dj1)))

for i in range(dj1):
    for j in range(dj1):
        J1[i,j]=dPdd(i,j,V,d)    

def J1(V,d):
    dj1=
