## Fonctions de base sur les polynômes (addition, multiplication, dérivation, etc...):

import numpy as np
import random as rd
from random import uniform
import matplotlib.pyplot as plt
from scipy import misc
# on envisage ici les polynômes en tant que liste de coefficients
#P+Q: max(deg(P),deg(Q)) opérations
def somme_polynomes(P,Q):
    somme=[]
    n=len(P)
    m=len(Q)
    if n<m:
        P+=(m-n)*[0]
    elif n<m:
        Q+=(n-m)*[0]
    deg_max=max(n,m)
    for k in range(0,deg_max):
        somme+=[P[k]+Q[k]]
    return somme
#aP: deg(P) opérations
def polynome_x_constante(P,a):
    n=len(P)
    for k in range (0,n):
        P[k]*=a
    return P
#PQ: (deg(P)+deg(Q))**2 opérations
def produit_polynomes(P,Q):
    n=len(P)
    m=len(Q)
    P+=m*[0]
    Q+=n*[0]
    p=n+m
    #donc deg(PQ)=n+m
    liste=[]
    for k in range (0,p-1):
        a_k=0
        for i in range (0,p-1):
            a_k+=P[i]*Q[k-i]
        liste+=[a_k]
    return(liste)
#P(x): deg(P) opérations
def image_polynome(P,x):
    n=len(P)
    sigma=0
    for k in range(0,n):
        sigma+=P[k]*x**k
    return sigma
#int(P): 3deg(P) opérations  
def integrale_polynome(P,a,b):
    n=len(P)
    Q=[]
    for k in range(0,n):
        Q+=[P[k]/(k+1)]
    Q=[0]+Q
    I=image_polynome(Q,b)-image_polynome(Q,a)
    return I
#P: deg(P) opérations
def derivee_polynome(P):
    n=len(P)
    dP=[]
    for k in range(1,n):
        dP+=[P[k]*(k)]
    return dP

def trace_polynomes(P,a,b):
    coordsx=np.linspace(a,b,1000)
    coordsy=[image_polynome(P,x) for x in coordsx]
    plt.clf
    plt.plot(coordsx,coordsy,"r-",linewidth=2)

#deux outils usuels pour le calcul de séries entières: la factorielle et les coefficients binomiaux

#complexité log(n)
def factorielle(n):
    if n==0:
        return 1
    else:
        return n*factorielle(n-1)
#complexité 3*log(n)
def coeff_binomiaux(k,n):
    u=factorielle(n)/(factorielle(k)*factorielle(n-k))
    return u