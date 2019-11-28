def Lagrange_Tchebychev_pts(xcoords,ycoords,n):
    n=n+1
    assert n<len(xcoords)
    a,b=min(xcoords),max(xcoords)
    racines=[(a+b)/2+(a-b)/2*np.cos((2*k+1)*(np.pi)/(2*n)) for k in range(0,n)]
    xint=[]
    yint=[]
    for k in racines:
        compteur=0
        m=len(xcoords)
        minx=xcoords[0]
        miny=ycoords[0]
        deltamin=k-minx
        for i in range (1,m):
            if abs(xcoords[i]-k)<deltamin:
                minx=xcoords[i]
                miny=ycoords[i]
                deltamin=k-minx
                compteur=i
        xint.append(minx)
        yint.append(miny)
        del(xcoords[compteur])
        del(ycoords[compteur])
    print(xint)
    trace_Lagrange(xint,yint)
    
    