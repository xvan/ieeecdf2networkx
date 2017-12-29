import numpy as np


Y=np.zeros((3,3),dtype=np.complex)


Y[0,1]=Y[1,0]=-np.complex(0.02,0.04)**-1
Y[0,2]=Y[2,0]=-np.complex(0.01,0.03)**-1
Y[1,2]=Y[2,1]=-np.complex(0.0125,0.025)**-1

Y[0,0]=-np.sum(Y[0,:])
Y[1,1]=-np.sum(Y[1,:])
Y[2,2]=-np.sum(Y[2,:])

Ymod=np.abs(Y)
Yang=np.angle(Y)


V=np.ones(3)
d=np.zeros(3)


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


def J(V,d):
    return np.matrix([
                [dPdd(1,1,V,d),dPdd(1,2,V,d), dPdV(1,1,V,d)],
                [dPdd(2,1,V,d),dPdd(2,2,V,d), dPdV(2,1,V,d)],
                [dQdd(1,1,V,d),dQdd(1,2,V,d), dQdd(2,2,V,d)]
                ])


S1sch = np.complex (-4,-2.5)
P2sch=2
V[0]=1.05
V[2]=1.04
d[0]=0


#Residuos (RP1, RP2, RQ1)
def R(V,d):
 return np.array(
            np.real(S1sch) - P(1,V,d),
            P2sch - P(2,V,d),
            np.imag(S1sch) - Q(1,V,d)
        )


 print J(V,d)
