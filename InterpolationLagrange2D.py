import imageio
import numpy as np
import matplotlib.pyplot as plt
Im=imageio.imread("/Users/luciendavid/Documents/ANIT.png")
plt.imshow(Im)
plt.show
##
#Interpoler les pixels manquants selon les trois teintes de couleurs
def interp_bilin(mat,m): #on calcule la valeur sur np points dans chaque rectangle de l'échantillon
    n,p=np.shape(mat)
    res=np.zeros((n*(m+1)-m,p*(m+1)-m))
    for i in range(n):
        for j in range(p):
            res[(m+1)*i][(m+1)*j]=mat[i][j]
    k0=0
    while k0<=n-2:
        l0=0
        while l0<=p-2:
            a,b,c,d=mat[k0][l0],mat[k0][l0+1],mat[k0+1][l0],mat[k0+1][l0+1]
            for i in range(k0,k0+m+1):
                for j in range(l0,l0+m+1):
                    if (i,j)!=(k0,l0) and (i,j)!=(k0,l0+m) and (i,j)!=(k0+m,l0) and (i,j)!=(k0+m,l0+m):
                        res[i][j]=f(i,j,a,b,c,d)
            l0+=1
        for i in range(k0,k0+m):
            if i!=k0 and i!=k0+m:
                res[i][p*(m+1)-m-1]=f(i,p*(m+1)-m-1,a,b,c,d)
        k0+=1
    for j in range(l0,l0+m):
        if j!=l0 and j!=l0+m:
            res[n*(m+1)-m-1][j]=f(n*(m+1)-m-1,j,a,b,c,d)
    return res

n=int(input("quelle largeur pour le carré?"))
res=[[[0,0,0] for k in range(n)] for l in range(n)]
for k in range(n):
    res[0][k],res[n-1][k]=[100,0,0],[0,0,100]
    
f0,f1,f2,f3=res[0][0],res[0][n-1],res[n-1][0],res[n-1][n-1]
def f(i,j,a,b,c,d):
    res=((a*(i-(n-1))-c*i)*(j-(n-1))+(b*(i-(n-1))-d*i)*j)/((n-1)**2)
    return res
w0,x0,y0,z0=f0[0],f1[0],f2[0],f3[0]
w1,x1,y1,z1=f0[1],f1[1],f2[1],f3[1]
w2,x2,y2,z2=f0[2],f1[2],f2[2],f3[2]
for i in range(1,n-1):
    for j in range(n):
        if (i,j)!=(0,0) and (i,j)!=(0,n-1) and (i,j)!=(n-1,0) and (i,j)!=(n-1,n-1):
            res[i][j]=[f(i,j,w0,x0,y0,z0),f(i,j,w1,x1,y1,z1),f(i,j,w2,x2,y2,z2)]
plt.imshow(res)
plt.show

