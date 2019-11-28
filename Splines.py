##Splines cubiques
#Lors de la résolution de ce problème, résolution d'un sytème matriciel à 4 inconnues. Utilisation d'une résolution "à la main" et d'une résolution en utilisant les outils de Python:

from scipy import misc
#complexité: linéaire
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
#complexité: linéaire
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