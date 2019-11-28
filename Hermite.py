##Interpolation d'Hermite
#Calcul du polynôme d'interpolation d'Hermite pour une liste de points et leurs dérivées associées données
#complexité: n**5 opérations
def polynome_lhermite(xcoords,ycoords,dycoords):
    P=[]
    n=len(xcoords)
    for i in range(n):
        Li=polynome_Lagrange_i(xcoords,i)
        dLi=derivee_polynome(Li)
        dLi_xi=image_polynome(dLi,xcoords[i])
        hi=produit_polynomes(somme_polynomes([1],polynome_x_constante([2*xcoords[i],-2],dLi_xi)),produit_polynomes(Li,Li))
        ki=produit_polynomes([-xcoords[i],1],produit_polynomes(Li,Li))
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
    coordsx=np.linspace(a-0.5,b+0.5,1000)
    coordsy=[image_polynome(P,x) for x in coordsx]
    fcoords=[f(x) for x in coordsx]
    plt.clf
    plt.plot(coordsx,fcoords,"b-",label="fonction à interpoler")
    plt.plot(coordsx,coordsy,"r-",label="polynôme interpolateur d'Hermite")
    plt.legend(loc="best")

def Lhermite_Tchebychev(f,a,b,n):
    n=n+1
    racines=[(a+b)/2+(a-b)/2*np.cos((2*k+1)*(np.pi)/(2*n)) for k in range(0,n)]
    ycoords=[f(x) for x in racines]
    dycoords=[misc.derivative(f,x) for x in racines]
    L=polynome_lhermite(racines,ycoords,dycoords)
    return L
    
def trace_lhermite_Tchebychev(f,a,b,n):
    n=n+1
    racines=[(a+b)/2+(a-b)/2*np.cos((2*k+1)*(np.pi)/(2*n)) for k in range(0,n)]
    ycoords=[f(x) for x in racines]
    dycoords=[misc.derivative(f,x) for x in racines]
    L=polynome_lhermite(racines,ycoords,dycoords)
    coordsx=np.linspace(a,b,1000)
    coords=[f(x) for x in coordsx]
    coordsy=[image_polynome(L,x) for x in coordsx]
    plt.clf
    plt.plot(coordsx,coords,"r-",linewidth=2,label="fonction à interpoler")
    plt.plot(coordsx,coordsy,"b-",label="polynôme interpolateur de Tchebychev")
    plt.legend(loc="best")
    plt.show()

def random_vs_TnLH(f,a,b,n):
    liste=[]
    LT=Lhermite_Tchebychev(f,a,b,n)
    for k in range(0,n+1):
        liste+=[uniform(a,b)]
    ycoords=[f(x) for x in liste]
    dycoords=[misc.derivative(f,x) for x in liste]
    L=polynome_lhermite(liste,ycoords,dycoords)
    coordsx=np.linspace(a,b,1000)
    coords=[f(x) for x in coordsx]
    coordsy=[image_polynome(L,x) for x in coordsx]
    coordsy2=[image_polynome(LT,x) for x in coordsx]
    plt.clf
    plt.plot(coordsx,coords,"b-",label="fonction à interpoler")
    plt.plot(coordsx,coordsy,"r--",label="polynôme interpolateur d'Hermite")
    plt.plot(coordsx,coordsy2,"g--",label="polynôme interpolateur Tchebychev")
    plt.legend(loc="best")
    plt.show()

def integrale_interp2(f,a,b,n):
    n=n+1
    racines=[(a+b)/2+(a-b)/2*np.cos((2*k+1)*(np.pi)/(2*n)) for k in range(0,n)]
    ycoords=[f(x) for x in racines]
    dycoords=[misc.derivative(f,x) for x in racines]
    L=polynome_lhermite(racines,ycoords,dycoords)
    return integrale_polynome(L,a,b)




