n=int(input("entrer n: "))
a=float(input("entrer a: "))
b=float(input("entrer b: "))
#SE usuelles tronquées à l'ordre n:
#exponentielle: l'intégrale vaut exp(1)-1
exp=[1/factorielle(k) for k in range(n)]
I1=integrale_polynome(exp,a,b)
I1_ref=np.exp(b)-np.exp(a)
#ch: l'intégrale vaut sh(1)
ch=[]
for k in range(n):
    if k%2==0:
        ch.append(1/factorielle(k))
    else:
        ch.append(0)
I2=integrale_polynome(ch,a,b)
I2_ref=np.sinh(b)-np.sinh(a)
#sh: l'intégrale vaut ch(1)-1
sh=[]
for k in range(n):
    if k%2==1:
        sh.append(1/factorielle(k))
    else:
        sh.append(0)
I3=integrale_polynome(sh,a,b)
I3_ref=np.cosh(b)-np.cosh(a)
#cos: l'intégrale vaut sin(1)
cos=[]
for k in range(n):
    if k%2==0:
        cos.append((-1)**(k/2)/factorielle(k))
    else:
        cos.append(0)
I4=integrale_polynome(cos,a,b)
I4_ref=np.sin(b)-np.sin(a)
#sin: l'intégrale vaut 1-cos(1)
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



