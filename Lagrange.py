## Interpolation de Lagrange, utilisation des abscisses de Tchebychev, formule de Newton-Côtes:
#Calcul de Li(dans la formule de l'analyse-synthèse du problème de Lagrange)
#Complexité: n**3 opérations
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
#Complexité: n**4 opérations
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
    coordsx=np.linspace(a-1,b+1 ,1000)
    coords=[f(x) for x in coordsx]
    coordsy=[image_polynome(L,x) for x in coordsx]
    plt.clf
    plt.plot(coordsx,coords,"r-",linewidth=2,label="fonction à interpoler")
    plt.plot(coordsx,coordsy,"b-",label="polynôme interpolateur de Lagrange")
    plt.legend(loc="best")
    plt.show
#Calcul du polynôme de Lagrange de DEGRE n (n+1 points d'interpolation) par les abscisses de Tchebychev
#complexité: n**4 opérations
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
    plt.plot(coordsx,coords,"b-",linewidth=2)
    plt.plot(coordsx,coordsy,"r-")
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
    trace_Lagrange_Tchebychev(f,a,b,n)
    
#Formule de Newton-Côtes utilisant la base des polynômes de Lagrange pour une approximation intégrale
#complexité:n**4 opérations
def Newton_Cotes(f,a,b,n):
    xcoords=[a+k*(b-a)/n for k in range(0,n)]
    sigma=0
    for i in range (0,n):
        li=polynome_Lagrange_i(xcoords,i)
        sigma+=f(xcoords[i])*integrale_polynome(li,a,b)
    return sigma
    #méthode divergente: reste correct pour n<25 en général. Méthode tout de même assez instable
    #on peut retrouver l'intégrale d'une fonction non primitivable usuellement, commme l'intégrale de Gauss, mais il ne faut pas prendre un intervalle trop grand. En effet, pour la primitive de Gauss, le fonction devient constante sur les bords de l'intervalle et le polynôme devient trop irrégulier pour être constant de la sorte tout en étant de degré 10,20...

#Approximation intégrale en traçant le polynôme de Tchebychev d'ordre n et en calculant son intégrale comme approximation de celle de f:
#complexité: n**4 opérations
def integrale_interp(f,a,b,n):
    n=n+1
    racines=[(a+b)/2+(a-b)/2*np.cos((2*k+1)*(np.pi)/(2*n)) for k in range(0,n)]
    ycoords=[f(x) for x in racines]
    L=polynome_Lagrange(racines,ycoords)
    return integrale_polynome(L,a,b)