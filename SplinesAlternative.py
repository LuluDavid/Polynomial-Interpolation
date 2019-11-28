def splines_cubiques_bis(f,a,b,n):
    xcoords=[a+k*(b-a)/n for k in range(n+1)]
    ycoords=[f(k) for k in xcoords]
    dycoords=[misc.derivative(f,k) for k in xcoords]
    for k in range(n):
        a0=xcoords[k]-xcoords[k+1]
        x,y,fx,fy,dfx,dfy=xcoords[k],xcoords[k+1],ycoords[k],ycoords[k+1],dycoords[k],dycoords[k+1]
        nu0=(dfx+dfy+2*(fx-fy)/a)/(a**2)
        nu1=(3*(x+y)*(fx-fy)-((x**3-3*x*y**2+2*y**3)*dfx+(y**3-3*y*x**2+2*x**3)*dfy)/a)/(a**3)
        nu2=(6*(y**2*x-x**2*y)*(fx-fy)+(y**4+2*y*x**3-3*(x*y)**2)*dfx+(x**4+2*x*y**3-3*(x*y)**2)*dfy)/(a**4)
        nu3=(-(x*y)(y*dfx+x*dfy)+((y**4-4*x*y**3+3*(x*y)**2)*dfx+(x**4-4*y*x**3+3*(x*y)**2)*dfy)/(a**2))/(a**2)
        P=[nu0,nu1,nu2,nu3]
        trace_polynomes(P,xcoords[k],xcoords[k+1])
    coordsx=np.linspace(a,b,1000)
    coordsy=[f(x) for x in coordsx]
    plt.plot(coordsx,coordsy,"b-",label="fonction originale")
    plt.legend(loc="best")