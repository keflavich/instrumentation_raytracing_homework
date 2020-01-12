"""
HW 1 - Cassegrain Telescope
Q1 - Where put a 30" focal plane?
at the focal point, centered
Q2 - Where put a 5' focal plane?
at the focal point, centered, as far as I can tell
Q3 - How big a secondary is needed? 71.3mm
abs( firstbounce(-0.5,0,300/206265.,0,sqrt(1-(300/206265.)**2),final=True)[0] )
"""
from pylab import *
import matplotlib
import mpl_toolkits.mplot3d.axes3d as axes3d
import pdb 
import numpy
from scipy.optimize.minpack import fsolve,newton
from scipy.optimize import bisect,brenth,fmin

#def sqrt(x):
#    f = lambda y: y**2 - x**2
#    return fsolve(f,numpy.sqrt(x),xtol=1e-18)

diam1 = float64(1.0)
rad1 = diam1/2.0
fr1 = 3
fl1 = fr1*diam1 # focal length
db = float64(-0.1)# distance back
floff = fl1+abs(db)
default_zoff = db + (floff)/float64(2) # 1.45
zoff = default_zoff
fr2 = float64(20.0)
# tangent slope = 2 k x = d/dx (z=k x^2)
# normal slope = -1 / 2kx
# normal line z = (-1 / (2 k x0)) * x + fl1
default_k = float64(1)/(fl1 + sqrt(fl1**2-2*rad1**2)) # 0.169  sqrt(1/k) = 2.432 sqrt(1000/k) = 76.912
k = default_k = float64(1/(4*fl1))
default_b = float64(0)
zparab = k*0.5**2
# solve z = normal line = -2*fr2 * x + db = (1/(2 k x0)) x + fl1
# x = (db - fl1) / (1/(2 k x0) + 2*fr2)
# (tangent vector is negative wrt x)
slope1 = float64(1)/(2*default_k*rad1)
ang = arctan(1.0/slope1) # cotan(1/(2*k*rad1)) = cot(normal slope)
newslope = 1/tan(ang*2)

xhyp = -1*(abs(db)+fl1) / (newslope + 2*fr2) # -0.067524
rad2 = abs(xhyp)
diam2 = rad2*float64(2) # 0.1349
zhyp = xhyp*newslope + fl1 # 2.598096
zhyp = -xhyp*2*fr2+db
#print "check: zhyp == -xhyp*fr2*2 + db",zhyp==-xhyp*2*fr2+db,-xhyp*2*fr2+db,zhyp,-xhyp*2*fr2+db-zhyp
qb = (zhyp-zoff)**2 + xhyp**2 - floff**2/4.0
qc = -xhyp**2*floff**2 / 4.0
a = sqrt((-qb+sqrt(qb**2-4*qc))/2.0)
default_a = a
d = sqrt((floff/2)**2-a**2)
default_d = d
#default_a = 1.0439613255199183
#default_d = 1.1457070964337681
delta = arctan(-2.*fr2)
#print "alpha: %0.5f  delta: %0.5f  " % ( alpha*180/pi , delta*180/pi )


def do_prob1(overlay=False):
    offsets=[0,1,3,5,10,30,60,120,240,600,1200,3000]
    figure(3)
    for o in offsets:
        if not overlay:
            clf()
        focalplane(o)
        xlabel('m')
        ylabel('m')
        if not overlay:
            savefig("spotdiagram_offset%i.png" % o)


def spot(nr=10,nth=60,square=False):
    if square:
        xx,yy = (indices([nr*2+1,nr*2+1]) - nr) / (2*float64(nr))
        ok = (xx**2+yy**2)**0.5 < 0.5
        x = xx[ok]
        y = yy[ok]
    else:
        r = linspace(0.1,0.5,nr) #arange(1,nr,dtype='float64')
        theta = linspace(0,2*pi,nth)
        x = outer(r,cos(theta))
        y = outer(r,sin(theta))


    #x /= 206265.0
    #y /= 206265.0

    #z = sqrt(float64(1.0)-x**2-y**2)
    

    return x.ravel(),y.ravel()


def projectthru(x0,y0,xi=0.0,yi=0.0,zi=1.0):
    x1,y1,z1 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,initial=True)
    x2,y2,z2 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,final=True)
    xv1,yv1,zv1 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,vector=True)
    x3,y3,z3 = secondbounce(x2,y2,z2,xv1,yv1,zv1)
    return x3,y3,z3



def mirror1(x,y,k=default_k,b=default_b,normal=False,tangent=False):
    """
    z = k(x^2+y^2)+b
    """
    x = float64(x)
    y = float64(y)
    if normal:
        fx = 2*k*x
        fy = 2*k*y
        norm = sqrt(fx**2+fy**2+1)
        xv,yv,zv = fx/norm,fy/norm,-1.0/norm
        if isnan(xv): pdb.set_trace()
        return xv,yv,zv
    if tangent:
        fx = 2*k*x
        fy = 2*k*y
        norm = sqrt(fx**2+fy**2+2)
        xv,yv,zv = 1/norm,1/norm,sqrt(fx**2+fy**2)/norm
        if isnan(xv): pdb.set_trace()
        return xv,yv,zv
    else:
        z = k*(x**2+y**2)+b
        return z

def mirror2(x,y,zoff=default_zoff,k=default_k,b=default_b,a=default_a,d=default_d,normal=False,tangent=False):
    """
    (z-zoff)^2/d^2 = (x^2+y^2)/a^2 + 1
    z = +/- d * sqrt( (x^2+y^2)/a^2 + 1 ) + zoff
    z = +/- d * sqrt( r^2/a^2 + 1 ) + zoff
    dz/dx = +/- d * 1/2 * ( (x^2+y^2)/a^2 + 1 )^(-1/2) * 2*x/a^2
    dz/dr = +/- d * 1/2 * ( r^2/a^2 + 1 )^(-1/2) * 2*r/a^2
    normal z = -dx/dz * x + b
    b = z + dx/dz*x for some point on parabola
      = +/- d * sqrt( (x^2+y^2)/a^2 + 1 ) + zoff +/- dx/dz*x
    """
    x=float64(x)
    y=float64(y)
    r=sqrt(x**2+y**2)
    try:
        if r == 0:
            return float64(0),float64(0),float64(-1)
    except:
        pass
    theta=arctan2(y,x)
    z1 = sqrt((r**2)/a**2 + 1) * d + zoff
    z = sqrt((x**2+y**2)/a**2 + 1) * d + zoff
    #print (z==z1).all()
    # first derivative of circular version of mirror eqn....
    fr = d/2.0 * 2*r/a**2 / sqrt((r**2)/a**2+1) 
    nslope = 1.0/fr
    nb = z + abs(nslope)*r
    # normalization of vector: |vn| = sqrt(z^2+r^2) = sqrt( (nslope*r)^2 + r^2 )
    norm = sqrt( (nslope*r)**2 + r**2 )
    vn = array((r*cos(theta),r*sin(theta),nslope*r))/norm
    vn = array((x,y,nslope*r))/norm
    ##print r*cos(theta)/x
    ##print r*sin(theta),y
    xv,yv,zv = vn

    if normal:
        #norm = sqrt(fx**2+fy**2+1)
        #xv,yv,zv = fx/norm,fy/norm,-1.0/norm
        return xv,yv,zv
    elif tangent:
        fx = d/2.0 * 2*x/a**2 / sqrt((x**2+y**2)/a**2+1) + zoff
        fy = d/2.0 * 2*y/a**2 / sqrt((x**2+y**2)/a**2+1) + zoff
        norm = sqrt(fx**2+fy**2+2)
        xv,yv,zv = 1/norm,1/norm,sqrt(fx**2+fy**2)/norm
        return xv,yv,zv
    else:
        return z

def firstbounce(x,y,xi=0,yi=0,zi=1.0,k=default_k,b=default_b,vector=False,a=default_a,d=default_d,
        initial=False,final=False,zoff=default_zoff):
    """
    assumed incoming parallel light

    returns output vector (no matter what vector=)
    """
    # initial position is on mirror 1
    xo = float64(x)
    yo = float64(y)
    zo = mirror1(x,y,k,b)

    #normal vector from mirror 1
    xn,yn,zn = mirror1(x,y,k,b,normal=True)
    vn = array([xn,yn,zn])

    #incident ray is assumed in z-direction for now
    vi = array([xi,yi,zi])

    crossi = cross(vi,vn)
    doti = dot(vi,vn)
    if xn==0 and yn==0:
        xr,yr,zr = 0,0,-1
    elif abs(xn) > abs(yn):
        xr = (doti + zn/xn*crossi[1] - yn/xn*crossi[2]) / (xn + yn**2/xn + zn**2/xn)
        yr = (crossi[2] + xr*yn) / xn
        zr = (-1*crossi[1] + xr*zn) / xn
    else:
        yr = (doti - zn/yn*crossi[0] + xn/yn*crossi[2]) / (yn + xn**2/yn + zn**2/yn)
        xr = (-1*crossi[2] + yr*xn) / yn
        zr = (crossi[0] + yr*zn) / yn
    vr = array([xr,yr,zr])
    vr /= sqrt(sum(vr**2))
    xr,yr,zr = vr

    # X1 * t*V = X2
    QA = zr**2 - d**2/a**2*(xr**2+yr**2)
    QB = 2*zr*(zo-zoff) - 2*d**2/a**2 * (xr*xo + yr*yo)
    QC = (zo-zoff)**2 - d**2/a**2 * (xo**2+yo**2) - d**2
    t = (-QB + sqrt(QB**2 - 4*QA*QC))/(2*QA)
    t1 = (-QB + sqrt(QB**2 - 4*QA*QC))/(2*QA)
    t2 = (-QB - sqrt(QB**2 - 4*QA*QC))/(2*QA)
    xf = xo + t*xr
    yf = yo + t*yr
    zf = zo + t*zr
 
    """
    Q = (1+(yi/xi)**2)
    cost = dot(vi,vn)
    xr = (cost*xi*(xi*xn+yi*yn)-
            xi*zn*sqrt( Q*xi**2*zn**2 + (xi*xn+yi*yn)**2 - cost**2*Q*xi**2 )
            ) / ( (xi*xn+yi*yn)**2 + Q*xi**2*zn**2 )
    yr = xr * yi/xi
    zr = (1-xr**2-yr**2)**0.5
    # #print dot(vi,vr),dot(vn,vr),dot(vi,vn)
    #print "t1: %0.3f  t2: %0.3f" % (t1,t2)
    #print "t1,t2 / real t",t1/2.6133084013297072,t2/2.6133084013297072
    #print "t/t1,t2       ",2.6133084013297072/t1,2.6133084013297072/t2
    #print "QA,QB,QC",QA,QB,QC,"a,d",a,d
    # #print "normalized?",sqrt(sum(vi**2)),sqrt(sum(vn**2)),sqrt(sum(vr**2))
    # #print "normal vector:     %0.3f %0.3f %0.3f" % (xn,yn,zn)
    # ##print "??: %0.3f %0.3f   %0.3f" % (xv,yv,zv)
    # #print "reflected vector:  %0.3f %0.3f %0.3f" % (xr,yr,zr)
    """
    alpha1 = 180/pi * (arctan2(zi,xi)-arctan2(zn,xn)) - 180
    alpha2 = 180/pi * (arctan2(zn,xn)-arctan2(zr,xr)) + 180
    #print "FIRST BOUNCE"
    #print "mirror 1 position: %0.3f %0.3f %0.3f" % (xi,yi,zi)
    #print "mirror 2 position:   %0.3f %0.3f %0.3f" % (xf,yf,zf)
    #print "crossI, dotI:",crossi,doti
    #print "reflected:   %0.3f %0.3f %0.3f" % (xr,yr,zr),vr,sqrt(sum(vr**2))
    #print "normal:   %0.3f %0.3f %0.3f" % (xn,yn,zn)
    #print "first bounce: normal dot reflected, inicident dot normal",dot(vn,vr),dot(vi,vn)
    #print "first bounce: normal alpha: %0.5f  reflected alpha: %0.5f" % (alpha1,alpha2)
    if isnan(alpha1) or isnan(alpha2):
        import pdb; pdb.set_trace()

    if isnan(xr) or isnan(xo) or isnan(xf) or isnan(xr): pdb.set_trace()
    if vector:
        return xr,yr,zr
    elif initial:
        return xo,yo,zo
    elif final:
        return xf,yf,zf
    else:
        return xr,yr,zr

def secondbounce(x,y,z,xi,yi,zi,a=default_a,b=default_b,c=1,k=default_k,vector=False,db=db):
    """
    Computes location and angle of bounce off of secondary
    x2,y2,z2 = t(xv,yv,zv) + x,y,z
    x,y,z should be given from firstbounce
    """
    
    xo,yo,zo=float64(x),float64(y),float64(z)

    #normal vector from mirror 2
    xn,yn,zn = mirror2(x,y,k,b,normal=True)
    zn = -abs(zn) # force bounce downward
    vn = array([xn,yn,zn])

    #incident ray is assumed in z-direction for now
    #if xi is None or yi is None or zi is None:
    #    vi = -1.0*array(firstbounce(x,y))
    #    xi,yi,zi = vi
    #else:
    vi = -1.0*array([xi,yi,zi])

    crossi = cross(vi,vn)
    doti = dot(vi,vn)
    if xn==0 and yn==0:
        yr = 0
        xr = 0
        zr = -1
    elif abs(xn) > abs(yn):
        xr = (doti + zn/xn*crossi[1] - yn/xn*crossi[2]) / (xn + yn**2/xn + zn**2/xn)
        yr = (crossi[2] + xr*yn) / xn
        zr = (-1*crossi[1] + xr*zn) / xn
    else:
        yr = (doti - zn/yn*crossi[0] + xn/yn*crossi[2]) / (yn + xn**2/yn + zn**2/yn)
        xr = (-1*crossi[2] + yr*xn) / yn
        zr = (crossi[0] + yr*zn) / yn
    vr = array([xr,yr,zr])
    vr /= sqrt(sum(vr**2))
    xr,yr,zr = vr

    zf = db
    t = (zf-zo)/zr
    xf = xo+t*xr
    yf = yo+t*yr

    """
    print "SECOND BOUNCE"
    print "incident dot reflected, normal dot reflected, iniciden dot normal",dot(vi,vr),dot(vn,vr),dot(vi,vn)
    print "crossI, dotI:",crossi,doti
    #print "t1: %0.3f  t2: %0.3f" % (t1,t2)
    print "t:",t
    #print "QA,QB,QC",QA,QB,QC,"a,d",a,d
    # print "normalized?",sqrt(sum(vi**2)),sqrt(sum(vn**2)),sqrt(sum(vr**2))
    print "mirror 2 position: %0.5f %0.5f %0.5f" % (xo,yo,zo)
    print "focal pl position:   %0.5f %0.5f %0.5f" % (xf,yf,zf)
    print "incident vector:     %0.5f %0.5f %0.5f slope: %0.5f" % (xi,yi,zi,zi/xi)
    print "normal vector:     %0.5f %0.5f %0.5f" % (xn,yn,zn)
    # #print "??: %0.5f %0.5f   %0.5f" % (xv,yv,zv)
    print "reflected vector:  %0.5f %0.5f %0.5f slope: %0.5f" % (xr,yr,zr,zr/xr)
    alpha1 = 180/pi * -(arctan2(zi,xi)-arctan2(-zn,-xn))
    alpha2 = 180/pi * -(arctan2(-zn,-xn)-arctan2(zr,xr)) + 180
    print "second bounce: normal alpha: %0.5f  reflected alpha: %0.5f" % (alpha1,alpha2)
    """
    #t= (-c**2*(x*xv+y*yv) + a**2*z*zv-\
    #        0.5*sqrt( 4*(c**2*(x*xv+y*yv)-a**2*z*zv)**2 - 
    #            4*(c**2*(a**2*d**2+x**2+y**2)-a**2*z**2)*
    #            (c**2*(xv**2+yv**2)-a**2*zv**2)))/\
    #            (c**2*(xv**2+yv**2)-a**2*zv**2)
    #xf = t*xv + x
    #yf = t*yv + y
    #zf = t*zv + z
    #figure(3); clf();
    #plot([0,xi],[0,zi],label='incident')
    #plot([0,xn],[0,zn],label='normal')
    #plot([0,xr],[0,zr],label='reflected')
    #legend(loc='best')
    if isnan(xr) or isnan(xf): pdb.set_trace()
    if vector:
        return xr,yr,zr
    else:
        return xf,yf,zf

"""
def gratingbounce(x,y,z,xi,yi,zi,vector=False,db=db,mg=1,lg=1200e-10,dg=1.0/3.600e6):
#    Returns the reflected vector from a grating at the focal plane

    xn,yn,zn = gratingnormal(x,y)
    vn = array([xn,yn,zn])
    vi = -array([xi,yi,zi])

    # error checking
    vn /= sqrt( (vn**2).sum() )
    vi /= sqrt( (vi**2).sum() )

    doti = dot(vi,vn)
    if doti > 1:
        doti = 1
    alpha = arccos(doti)
    beta = arcsin( mg*lg/dg - sin(alpha) )

    # force reflection in x-axis
    yo = 0
    xo = cos(beta)
    zo = sin(beta)
    vo = array([xo,yo,zo])
    vo /= sqrt( (vo**2).sum() )
    xo,yo,zo = vo

    print "alpha: %f  beta: %f" % (alpha*180/pi,beta*180/pi)
    if isnan(alpha) or isnan(beta):
        import pdb; pdb.set_trace()

    if isnan(vo): pdb.set_trace()
    return vo
"""


