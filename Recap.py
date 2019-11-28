## Fonctions de base sur les polynômes (addition, multiplication, dérivation, etc...):

import numpy as np
import random as rd
from random import uniform
import matplotlib.pyplot as plt
import scipy as sp
from scipy import misc
# on envisage ici les polynômes en tant que liste de coefficients
#P+Q
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
#aP
def polynome_x_constante(P,a):
    n=len(P)
    for k in range (0,n):
        P[k]*=a
    return P
#PQ
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
#P(x)
def image_polynome(P,x):
    n=len(P)
    sigma=0
    for k in range(0,n):
        sigma+=P[k]*x**k
    return sigma
#int(P)    
def integrale_polynome(P,a,b):
    n=len(P)
    Q=[]
    for k in range(0,n):
        Q+=[P[k]/(k+1)]
    Q=[0]+Q
    I=image_polynome(Q,b)-image_polynome(Q,a)
    return I
#P'
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
def factorielle(n):
    u=1
    for k in range(n):
        u=u*(k+1)
    return u

def coeff_binomiaux(k,n):
    u=factorielle(n)/(factorielle(k)*factorielle(n-k))
    return u
    
## Interpolation de Lagrange, utilisation des abscisses de Tchebychev, formule de Newton-Côtes:
#Calcul de Li(dans la formule de l'analyse-synthèse du problème de Lagrange)
def polynome_Lagrange_i(xcoords,i):
    # prend en entrée les abscisses d'interpolation, renvoie les coefficients du ième polynôme d'interpolation
    li=[1]
    n=len(xcoords)
    for k in range(0,n):
        if k!=i:
            Pj=[-xcoords[k]/(xcoords[i]-xcoords[k]),1/(xcoords[i]-xcoords[k])]
            li=produit_polynomes(li,Pj)
    return li
#Calcul du polynôme complet de Lagrange (formule de l'analyse synthèse)
def polynome_Lagrange(xcoords,ycoords):
    L=[0]
    n=len(xcoords)
    for i in range (0,n):
        Li=polynome_Lagrange_i(xcoords,i)
        Li=polynome_x_constante(Li,ycoords[i])
        L=somme_polynomes(L,Li)
    return L
#Tracé de ce polynôme
def trace_Lagrange(xcoords,ycoords):
    L=polynome_Lagrange(xcoords,ycoords)
    a,b=min(xcoords),max(xcoords)
    coordsx=np.linspace(a,b,1000)
    coordsy=[image_polynome(L,x) for x in coordsx]
    plt.plot(coordsx,coordsy,"r-")
    plt.show
#Tracé du polynôme de Lagrange pour des points choisis d'une fonction
def trace_Lagrange_fonction(xcoords,f):
    ycoords=[f(x) for x in xcoords]
    L=polynome_Lagrange(xcoords,ycoords)
    a,b=min(xcoords),max(xcoords)
    coordsx=np.linspace(a,b,1000)
    coords=[f(x) for x in coordsx]
    coordsy=[image_polynome(L,x) for x in coordsx]
    plt.clf
    plt.plot(coordsx,coords,"g-",label="fonction à interpoler")
    plt.plot(coordsx,coordsy,"r--",label="polynôme interpolateur de Lagrange")
    plt.legend(loc="best")
    plt.show
#Calcul du polynôme de Lagrange de DEGRE n (n+1 points d'interpolation) par les abscisses de Tchebychev
def Lagrange_Tchebychev(f,a,b,n):
    n=n+1
    racines=[(a+b)/2+(a-b)/2*np.cos((2*k+1)*(np.pi)/(2*n)) for k in range(0,n)]
    ycoords=[f(x) for x in racines]
    L=polynome_Lagrange(racines,ycoords)
#Tracé de ce polynôme
def trace_Lagrange_Tchebychev(f,a,b,n):
    n=n+1
    racines=[(a+b)/2+(a-b)/2*np.cos((2*k+1)*(np.pi)/(2*n)) for k in range(0,n)]
    ycoords=[f(x) for x in racines]
    L=polynome_Lagrange(racines,ycoords)
    coordsx=np.linspace(a,b,1000)
    coords=[f(x) for x in coordsx]
    coordsy=[image_polynome(L,x) for x in coordsx]
    plt.clf
    plt.plot(coordsx,coords,"g-",label="fonction à interpoler")
    plt.plot(coordsx,coordsy,"b--",label="polynôme interpolateur de Tchebychev")
    plt.legend(loc="best")
    plt.show()
#Comparaison du tracé d'un polynôme de Lagrange de degré n avec des points hasardeux avec celui des points de Tchebychev:
def random_vs_Tn(f,n,a,b):
    liste=[]
    for k in range(0,n+1):
        liste+=[uniform(a,b)]
    ycoords=[f(x) for x in liste]
    L=polynome_Lagrange(liste,ycoords)
    coordsx=np.linspace(a,b,1000)
    coords=[f(x) for x in coordsx]
    coordsy=[image_polynome(L,x) for x in coordsx]
    plt.clf
    plt.plot(coordsx,coordsy,"r--",label="polynôme interpolateur de Lagrange")
    tracé_Lagrange_Tchebychev(f,a,b,n)
    
#Formule de Newton-Côtes utilisant la base des polynômes de Lagrange pour une approximation intégrale
def Newton_Cotes(f,n,a,b):
    xcoords=[a+k*(b-a)/n for k in range(0,n)]
    sigma=0
    for i in range (0,n):
        li=polynome_Lagrange_i(xcoords,i)
        sigma+=f(xcoords[i])*integrale_polynome(li,a,b)
    return sigma
    #méthode divergente: reste correct pour n<25 en général. Méthode tout de même assez instable
    #on peut retrouver l'intégrale d'une fonction non primitivable usuellement, commme l'intégrale de Gauss, mais il ne faut pas prendre un intervalle trop grand. En effet, pour la primitive de Gauss, le fonction devient constante sur les bords de l'intervalle et le polynôme devient trop irrégulier pour être constant de la sorte tout en étant de degré 10,20...

#Approximation intégrale en traçant le polynôme de Tchebychev d'ordre n et en calculant son intégrale comme approximation de celle de f:
def integrale_interp(f,a,b,n):
    n=n+1
    racines=[(a+b)/2+(a-b)/2*np.cos((2*k+1)*(np.pi)/(2*n)) for k in range(0,n)]
    ycoords=[f(x) for x in racines]
    L=polynome_Lagrange(racines,ycoords)
    return integrale_polynome(L,a,b)

##Interpolation d'Hermite
#Calcul du polynôme d'interpolation d'Hermite pour une liste de points et leurs dérivées associées données
def polynome_lhermite(xcoords,ycoords,dycoords):
    P=[]
    n=len(xcoords)
    for i in range(n):
        Li=polynome_Lagrange_i(xcoords,i)
        dLi=derivee_polynôme(Li)
        dLi_xi=image_polynome(dLi,xcoords[i])
        hi=produit_polynomes(somme_polynomes([1],polynome_x_constante([2*xcoords[i],-2],dLi_xi)),produit_polynomes(Li,Li))
        ki=produit_polynomes([-xcoords[i],1],produit_polynômes(Li,Li))
        P=somme_polynomes(P,somme_polynomes(polynome_x_constante(hi,ycoords[i]),polynome_x_constante(ki,dycoords[i])))
    return P
#Tracé de ce polynôme
def trace_lhermite(xcoords,ycoords,dycoords):
    P=polynome_lhermite(xcoords,ycoords,dycoords)
    a,b=min(xcoords),max(xcoords)
    coordsx=np.linspace(a,b,1000)
    coordsy=[image_polynôme(P,x) for x in coordsx]
    plt.plot(coordsx,coordsy,"r-")
#Utilisation des modules pythons pour déterminer les dérivées en les points choisis pour une fonction donnée et interpolation à l'aide de ces listes de points:
def trace_lhermite_fonction(xcoords,f):
    ycoords=[f(x) for x in xcoords]
    dycoords=[misc.derivative(f,x) for x in xcoords]
    P=polynome_lhermite(xcoords,ycoords,dycoords)
    a,b=min(xcoords),max(xcoords)
    coordsx=np.linspace(a,b,1000)
    coordsy=[image_polynome(P,x) for x in coordsx]
    fcoords=[f(x) for x in coordsx]
    plt.clf
    plt.plot(coordsx,fcoords,"b-",label="fonction à interpoler")
    plt.plot(coordsx,coordsy,"r-",label="polynôme interpolateur d'Hermite")
    plt.legend(loc="best")
#Comparaison de l'efficacité d'interpolation d'Hermite par les abscisses de Tchebychev et par des abscisses normales: à approfondir (utiliser majoration de l'erreur)
def trace_lhermite_subdivise(f,a,b,n):
    xcoords=[a+k*(b-a)/n for k in range(0,n)]
    trace_lhermite_fonction(xcoords,f)

def trace_lhermite_Tchebychev(f,a,b,n):
    racines=[(a+b)/2+(a-b)/2*np.cos((2*k+1)*(np.pi)/(2*n)) for k in range(0,n)]
    trace_lhermite_fonction(racines,f)

##Splines cubiques
#Lors de la résolution de ce problème, résolution d'un sytème matriciel à 4 inconnues. Utilisation d'une résolution "à la main" et d'une résolution en utilisant les outils de Python:
def splines_cubiques(f,a,b,n):
    xcoords=[a+k*(b-a)/n for k in range(n+1)]
    ycoords=[f(k) for k in xcoords]
    dycoords=[misc.derivative(f,k) for k in xcoords]
    for k in range(n):
        a0=ycoords[k+1]/(xcoords[k+1]-xcoords[k])
        a1=ycoords[k]/(xcoords[k]-xcoords[k+1])
        a2=(dycoords[k]-a0-a1)/((xcoords[k]-xcoords[k+1])**2)
        a3=(dycoords[k+1]-a0-a1)/((xcoords[k+1]-xcoords[k])**2)
        nu0=-a0*xcoords[k]-a1*xcoords[k+1]-a2*xcoords[k]*xcoords[k+1]**2-a3*xcoords[k]**2*xcoords[k+1]
        nu1=a0+a1+a2*(xcoords[k+1]**2+2*xcoords[k+1]*xcoords[k])+a3*(xcoords[k]**2+2*xcoords[k]*xcoords[k+1])
        nu2=-a2*(xcoords[k]+2*xcoords[k+1])-a3*(xcoords[k+1]+2*xcoords[k])
        nu3=a2+a3
        P=[nu0,nu1,nu2,nu3]
        trace_polynomes(P,xcoords[k],xcoords[k+1])
    coordsx=np.linspace(a,b,1000)
    coordsy=[f(x) for x in coordsx]
    plt.plot(coordsx,coordsy,"b-",label="fonction originale")
    plt.legend(loc="best")

def splines_cubiques2(f,a,b,n):
    xcoords=[a+k*(b-a)/n for k in range(n+1)]
    ycoords=[f(k) for k in xcoords]
    dycoords=[misc.derivative(f,k) for k in xcoords]
    for k in range(n):
        M=np.array([(1,xcoords[k],xcoords[k]**2,xcoords[k]**3),(1,xcoords[k+1],xcoords[k+1]**2,xcoords[k+1]**3),(0,1,2*xcoords[k],3*xcoords[k]**2),(0,1,2*xcoords[k+1],3*xcoords[k+1]**2)])
        N=np.array([(ycoords[k]),(ycoords[k+1]),(dycoords[k]),(dycoords[k+1])])
        O=np.dot(np.linalg.inv(M),N)
        a0,a1,a2,a3=O[0],O[1],O[2],O[3]
        P=[a0,a1,a2,a3]
        trace_polynomes(P,xcoords[k],xcoords[k+1])
    coordsx=np.linspace(a,b,1000)
    coordsy=[f(x) for x in coordsx]
    plt.plot(coordsx,coordsy,"b-",label="fonction originale")
    plt.legend(loc="best")
#Approximation intégrale par les splines
def interp_splines(f,a,b,n):
    I=0
    xcoords=[a+k*(b-a)/n for k in range(n+1)]
    ycoords=[f(k) for k in xcoords]
    dycoords=[misc.derivative(f,k) for k in xcoords]
    for k in range(n):
        a0=ycoords[k+1]/(xcoords[k+1]-xcoords[k])
        a1=ycoords[k]/(xcoords[k]-xcoords[k+1])
        a2=(dycoords[k]-a0-a1)/((xcoords[k]-xcoords[k+1])**2)
        a3=(dycoords[k+1]-a0-a1)/((xcoords[k+1]-xcoords[k])**2)
        nu0=-a0*xcoords[k]-a1*xcoords[k+1]-a2*xcoords[k]*xcoords[k+1]**2-a3*xcoords[k]**2*xcoords[k+1]
        nu1=a0+a1+a2*(xcoords[k+1]**2+2*xcoords[k+1]*xcoords[k])+a3*(xcoords[k]**2+2*xcoords[k]*xcoords[k+1])
        nu2=-a2*(xcoords[k]+2*xcoords[k+1])-a3*(xcoords[k+1]+2*xcoords[k])
        nu3=a2+a3
        P=[nu0,nu1,nu2,nu3]
        I+=integrale_polynome(P,xcoords[k],xcoords[k+1])
    return I
#comparaison de l'efficacité (complexité temporelle)des deux méthodes sur des fonctions usuelles:
def test_splines_cubiques(a,b,n):
    liste_f_ref=[np.sin,np.cos,np.tan,np.exp,np.cosh,np.sinh,np.tanh]
    time1,time2=[],[]
    p=len(liste_f_ref)
    for i in range(p):
        #méthode1
        f=liste_f_ref[i]
        début1=time.time()
        xcoords=[a+k*(b-a)/n for k in range(n+1)]
        ycoords=[f(k) for k in xcoords]
        dycoords=[misc.derivative(f,k) for k in xcoords]
        for k in range(n):
            a0=ycoords[k+1]/(xcoords[k+1]-xcoords[k])
            a1=ycoords[k]/(xcoords[k]-xcoords[k+1])
            a2=(dycoords[k]-a0-a1)/((xcoords[k]-xcoords[k+1])**2)
            a3=(dycoords[k+1]-a0-a1)/((xcoords[k+1]-xcoords[k])**2)
            nu0=-a0*xcoords[k]-a1*xcoords[k+1]-a2*xcoords[k]*xcoords[k+1]**2-a3*xcoords[k]**2*xcoords[k+1]
            nu1=a0+a1+a2*(xcoords[k+1]**2+2*xcoords[k+1]*xcoords[k])+a3*(xcoords[k]**2+2*xcoords[k]*xcoords[k+1])
            nu2=-a2*(xcoords[k]+2*xcoords[k+1])-a3*(xcoords[k+1]+2*xcoords[k])
            nu3=a2+a3
            P=[nu0,nu1,nu2,nu3]
        fin1=time.time()
        time1+=[fin1-début1]
        #méthode2
        début2=time.time()
        xcoords=[a+k*(b-a)/n for k in range(n+1)]
        ycoords=[f(k) for k in xcoords]
        dycoords=[misc.derivative(f,k) for k in xcoords]
        for k in range(n):
            M=np.array([(1,xcoords[k],xcoords[k]**2,xcoords[k]**3),(1,xcoords[k+1],xcoords[k+1]**2,xcoords[k+1]**3),(0,1,2*xcoords[k],3*xcoords[k]**2),(0,1,2*xcoords[k+1],3*xcoords[k+1]**2)])
            N=np.array([(ycoords[k]),(ycoords[k+1]),(dycoords[k]),(dycoords[k+1])])
            O=np.dot(np.linalg.inv(M),N)
            a0,a1,a2,a3=O[0],O[1],O[2],O[3]
            P=[a0,a1,a2,a3]
        fin2=time.time()
        time2+=[fin2-début2]
        plt.clf
    sigma1,sigma2=0,0
    n=len(time1)
    print(time1,time2)
    for k in range(n):
        sigma1+=time1[k]
        sigma2+=time2[k]
    print(sigma1/n,sigma2/n)
# On trouve que la méthode 1 est plus rapide que le seconde en moyenne sur ces quelques fonctions de référence. Changer n n'influe en aucun cas sur ce résultat. Le choix de a et b non plus. La méthode 1 est plus efficace pour toutes les fonctions.
  
## Comparaison des séries entières tronquées, de Lagrange-Tchebychev et des splines pour l'approximation d'intégrales connues:
n=int(input("entrer n: "))
a=float(input("entrer a: "))
b=float(input("entrer b: "))
#SE usuelles tronquées à l'ordre n:
#exponentielle: l'intégrale vaut exp(b)-exp(a)
exp=[1/factorielle(k) for k in range(n)]
I1=integrale_polynome(exp,a,b)
I1_ref=np.exp(b)-np.exp(a)
#ch: l'intégrale vaut sh(b)-sh(a)
ch=[]
for k in range(n):
    if k%2==0:
        ch.append(1/factorielle(k))
    else:
        ch.append(0)
I2=integrale_polynome(ch,a,b)
I2_ref=np.sinh(b)-np.sinh(a)
#sh: l'intégrale vaut ch(b)-ch(a)
sh=[]
for k in range(n):
    if k%2==1:
        sh.append(1/factorielle(k))
    else:
        sh.append(0)
I3=integrale_polynome(sh,a,b)
I3_ref=np.cosh(b)-np.cosh(a)
#cos: l'intégrale vaut sin(b)-sin(a)
cos=[]
for k in range(n):
    if k%2==0:
        cos.append((-1)**(k/2)/factorielle(k))
    else:
        cos.append(0)
I4=integrale_polynome(cos,a,b)
I4_ref=np.sin(b)-np.sin(a)
#sin: l'intégrale vaut cos(a)-cos(b)
sin=[]
for k in range(n):
    if k%2==1:
        sin.append((-1)**((k-1)/2)/factorielle(k))
    else:
        sin.append(0)
I5=integrale_polynome(sin,a,b)
I5_ref=np.cos(a)-np.cos(b)

#méthode de Tchebychev à l'ordre n:
I1T=integrale_interp(np.exp,a,b,n)
I2T=integrale_interp(np.cosh,a,b,n)
I3T=integrale_interp(np.sinh,a,b,n)
I4T=integrale_interp(np.cos,a,b,n)
I5T=integrale_interp(np.sin,a,b,n)

#méthode des splines cubiques à l'ordre n:
I1S=interp_splines(np.exp,a,b,n)
I2S=interp_splines(np.cosh,a,b,n)
I3S=interp_splines(np.sinh,a,b,n)
I4S=interp_splines(np.cos,a,b,n)
I5S=interp_splines(np.sin,a,b,n)

print("erreur SE :")
print(" -exp: ",abs(I1_ref-I1))
print(" -ch: ",abs(I2_ref-I2))
print(" -sh: ",abs(I3_ref-I3))
print(" -cos: ",abs(I4_ref-I4))
print(" -sin: ",abs(I5_ref-I5))
print("erreur Tchebychev :")
print(" -exp: ",abs(I1_ref-I1T))
print(" -ch: ",abs(I2_ref-I2T))
print(" -sh: ",abs(I3_ref-I3T))
print(" -cos: ",abs(I4_ref-I4T))
print(" -sin: ",abs(I5_ref-I5T))
print("erreur splines :")
print(" -exp: ",abs(I1_ref-I1S))
print(" -ch: ",abs(I2_ref-I2S))
print(" -sh: ",abs(I3_ref-I3S))
print(" -cos: ",abs(I4_ref-I4S))
print(" -sin: ",abs(I5_ref-I5S))






