from decimal import *
import numpy as np
def derivativen(f,x,n,epsilon=0.001):
    M=[[k**j for k in range(n+1)] for j in range (n+1)]
    B=[[0] for k in range (n)]+[[1]]
    res=resolution_systeme(M,B)
    ftot=Decimal(0)
    for k in range (n+1):
        fk=Decimal(f(x+k*epsilon))
        ftot+=Decimal(res[k][0])*fk
    ftot=ftot*Decimal(factorielle(n))/(Decimal(epsilon)**n)
    return ftot

def serie_entiere(f,n):
    P=[derivativen(f,0,k) for k in range(n)]

#getcontext.prec()=200
## Inversion de matrices (pour inverser Vandermonde):
def permutation_lignes(A,i,j):
    p=np.shape(A)[1]
    for k in range(p):
        A[i][k],A[j][k]=A[j][k],A[i][k]

def resolution_systeme(A0,B0):
    A=np.copy(A0)
    Y=np.copy(B0)
    n,p=np.shape(A)
    nb_1,nb_c=np.shape(Y)          
    for j in range(n-1):
        m=abs(A[j,j])
        i_pivot=j
        for i in range(j+1,n):
            if abs(A[i,j])>m:
                m=abs(A[i,j])
                i_pivot=i  
        permutation_lignes(A,i_pivot,j)
        permutation_lignes(Y,i_pivot,j)
        if A[j,j]!=0:
            for k in range(j+1,n):
                Y[k,:]-=A[k,j]/A[j,j]*Y[j,:]
                A[k,:]-=A[k,j]/A[j,j]*A[j,:]
    X=np.zeros((n,nb_c))
    for i in range(n-1,-1,-1):
        for k in range(i+1,n):
            Y[i,:]-=A[i,k]*X[k,:]
        X[i,:]=Y[i,:]/A[i,i]
    return X


