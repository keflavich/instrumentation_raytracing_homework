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


diam1 = float64(1.0)
rad1 = diam1/2.0
fr1 = 3
fl1 = fr1*diam1 # focal length
db = -0.1 # distance back
floff = fl1+abs(db)
default_zoff = db + (floff)/2 # 1.45
zoff = default_zoff
fr2 = float64(20.0)
# tangent slope = 2 k x = d/dx (z=k x^2)
# normal slope = -1 / 2kx
# normal line z = (-1 / (2 k x0)) * x + fl1
default_k = 1/(fl1 + sqrt(fl1**2-2*rad1**2)) # 0.169  sqrt(1/k) = 2.432 sqrt(1000/k) = 76.912
k = default_k = 1/(4*fl1)
default_b = 0
zparab = k*0.5**2
# solve z = normal line = -2*fr2 * x + db = (1/(2 k x0)) x + fl1
# x = (db - fl1) / (1/(2 k x0) + 2*fr2)
# (tangent vector is negative wrt x)
slope1 = 1/(2*default_k*rad1)
alpha = arctan(1.0/slope1) # cotan(1/(2*k*rad1)) = cot(normal slope)
newslope = 1/tan(alpha*2)

xhyp = -1*(abs(db)+fl1) / (newslope + 2*fr2) # -0.067524
rad2 = abs(xhyp)
diam2 = rad2*2 # 0.1349
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

def spot(nr=10,nth=60):
    r = linspace(0.1,0.5,nr) #arange(1,nr,dtype='float64')
    theta = linspace(0,2*pi,nth)
    x = outer(r,cos(theta))
    y = outer(r,sin(theta))

    #x /= 206265.0
    #y /= 206265.0

    #z = sqrt(float64(1.0)-x**2-y**2)
    

    return x.ravel(),y.ravel()


def plotmirrors(plotx=True):
    x = linspace(-0.5,0.5,1000)
    fig=figure(0,figsize=(8,8))
    #fig.add_axes([0.05,0.05,0.95,0.95])
    plot(x,mirror1(x,x*0))
    x2 = linspace(-rad2,rad2,1000)
    plot(x2,mirror2(x2,x2*0))
    plot(x2,grating(x2,x2*0))
    #fig.axes[0].set_aspect(1)
    xv,yv,zv = mirror1(-0.5,0,normal=True)
    xv,yv,zv = firstbounce(-0.5,0,vector=True)
    z0 = mirror1(-0.5,0)
    #plot([-0.5,xv*3-0.5],[z0,z0+zv*3])
    #plot([0.5,-xv*3+0.5],[z0,z0+zv*3])
    if not plotx:
        plot(cos(linspace(0,2*pi,5000)),sin(linspace(0,2*pi,5000))-1.1)
    fig.axes[0].axis([-2.55,2.55,-2.1,3])
    fig1=figure(1,figsize=(8,8))
    #clf()
    plot(x2,mirror2(x2,x2*0))
    #plot([-0.5,xv*3-0.5],[z0,z0+zv*3])
    #plot([0.5,-xv*3+0.5],[z0,z0+zv*3])
    fig1.axes[0].axis([-0.15,0.15,2.45,2.75])

def plotray(x0,y0,xi=0.0,yi=0.0,zi=1.0,plotx=True,gorder=1,lam=1200e-10,plot3d=False,f0=None,ax=None):

    x1,y1,z1 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,initial=True)
    x2,y2,z2 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,final=True)
    xv1,yv1,zv1 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,vector=True)
    x3,y3,z3 = secondbounce(x2,y2,z2,xv1,yv1,zv1)
    xv2,yv2,zv2 = secondbounce(x2,y2,z2,xv1,yv1,zv1,vector=True)
    t = focus_to_grating(x3,y3,z3,xv2,yv2,zv2)
    x4 = x3 + t * xv2
    y4 = y3 + t * yv2
    z4 = z3 + t * zv2
    xv3,yv3,zv3 = gratingbounce_jim(x4,y4,z4,xv2,yv2,zv2,lg=lam,mg=gorder,dg=groovespacing(x4,y4))
    #xv3,yv3,zv3 = gratingbounce(x3,y3,z3,xv2,yv2,zv2,lg=lam,mg=gorder)
    xn3,yn3,zn3 = gratingnormal(x4,y4)
    x5,y5,z5 = spectral_detector(x4,y4,z4,xv3,yv3,zv3)
    plotmirrors(plotx=plotx)
    if f0 is None:
        f0=figure(0)
    #plot([0,xhyp],[db,db-xhyp*2*fr2],label='slope40')
    xn,yn,zn = mirror2(x2,y2,normal=True)
    if plot3d:
        if ax is None:
            ax=axes3d.Axes3D(f0)
        x3d,y3d = ( indices([11,11]) - 5.0 ) / 50.0
        z3d = grating(x3d,y3d)
        ax.plot_wireframe(x3d,y3d,z3d)
        x3d,y3d = ( indices([11,11]) - 5.0 ) / 10.0
        z3d = mirror1(x3d,y3d)
        ax.plot_wireframe(x3d,y3d,z3d)
        x3d,y3d = ( indices([11,11]) - 5.0 ) / 50.0
        z3d = mirror2(x3d,y3d)
        ax.plot_wireframe(x3d,y3d,z3d)
        ax.plot3D([x1,x2,x3,x4,x5],[y1,y2,y3,y4,y5],[z1,z2,z3,z4,z5])
        return ax
    elif plotx: 
        plot([x1,x2,x3,x4,x5],[z1,z2,z3,z4,z5],label='bounces')
        plot([x1],[z1],'o')
        plot([x2,x2+xn],[z2,z2-zn])
        plot([x4,x4+xn3*0.5],[z4,z4+zn3*0.5])
        gca().axis([-2.55,2.55,-2.1,3])
        legend(loc='best')
    else:
        plot([y1,y2,y3,y4,y5],[z1,z2,z3,z4,z5],label='bounces')
        plot([y1],[z1],'o')
        plot([y2,y2+yn],[z2,z2-zn])
        plot([y4,y4+yv3],[z4,z4+zv3])
        gca().axis([-2.55,2.55,-2.1,3])
        legend(loc='best')
    #gca().axis([-0.15,0.15,2.45,2.75])
    #plot([x0,x0+xv*

def projectthru(x0,y0,xi=0.0,yi=0.0,zi=1.0):
    x1,y1,z1 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,initial=True)
    x2,y2,z2 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,final=True)
    xv1,yv1,zv1 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,vector=True)
    x3,y3,z3 = secondbounce(x2,y2,z2,xv1,yv1,zv1)
    return x3,y3,z3

def focalplane(offset,fignum=3,grating=False,lam=1200e-10,gorder=1):
    figure(fignum)
    xyin = spot()
    xout,yout,zout = zeros(xyin[0].shape),zeros(xyin[0].shape),zeros(xyin[0].shape)
    xinc,yinc,zinc = (offset/206265.,0,sqrt(1-(offset/206265.)**2))
    if grating: 
        for i,xyi in enumerate(array(xyin).T):
            xout[i],yout[i],zout[i] = projectthru_grating(xyi[0],xyi[1],xinc,yinc,zinc,lam=lam,gorder=gorder,final=True)
    else:
        for i,xyi in enumerate(array(xyin).T):
            xout[i],yout[i],zout[i] = projectthru(xyi[0],xyi[1],xinc,yinc,zinc)
    plot(xout,yout,',')


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
        return xv,yv,zv
    if tangent:
        fx = 2*k*x
        fy = 2*k*y
        norm = sqrt(fx**2+fy**2+2)
        xv,yv,zv = 1/norm,1/norm,sqrt(fx**2+fy**2)/norm
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
    theta=arctan2(y,x)
    z1 = sqrt((r**2)/a**2 + 1) * d + zoff
    z = sqrt((x**2+y**2)/a**2 + 1) * d + zoff
    #print (z==z1).all()
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
    if abs(xn) > abs(yn):
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
    if vector:
        return xr,yr,zr
    else:
        return xf,yf,zf

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

    return vo


def gratingbounce_jim(x,y,z,xi,yi,zi,vector=False,mg=1,lg=1200e-10,dg=1.0/3.600e6):
#	m = order
#	d = line spacing (angstroms/line)
#	l = wavelength (angstroms)
#	N is plane normal, pointing up on side that rays hit
#	G is vector in direction of lines on grating
#	I is incident vector pointing in ray prop. direction
#	o is output vector pointing in ray prop. direction
#	x indicates primed coordinate system
#	xx indicates double primed coordinate system

    n1,n2,n3 = gratingnormal(x,y)
    g1,g2,g3 = gratingline(x,y)
    i1,i2,i3 = xi,yi,zi
    m = mg 
    d = dg
    l = lg

    a = sqrt ( 1.0 - g3*g3 )

    xn1 = ( n1*g2 - n2*g1 )/a
    xn2 = g3*( n1*g1 + n2*g2 )/a  - a*n3

#	R is the transform from basic to xx coordinates
#	in xx coordinates, g is in the z direction
#	RT is the return transform

    r11 = ( xn2*g2 - xn1*g3*g1 )/a
    r12 = -( xn2*g1 + xn1*g2*g3 )/a
    r13 = xn1*a
    r21 = (xn1*g2 + xn2*g3*g1)/a
    r22 = ( -xn1*g1 + xn2*g2*g3 )/a
    r23 = -a*xn2
    r31 = g1
    r32 = g2
    r33 = g3

    rt11 = r11
    rt12 = r21
    rt13 = r31
    rt21 = r12
    rt22 = r22
    rt23 = r32
    rt31 = r13
    rt32 = r23
    rt33 = r33

    xxi1 = i1*r11 + i2*r12 + i3*r13

    if (xxi1 < 0) and abs(xxi1)>1e-15:
        print "error xxi1 is negative: ",xxi1

    xxi3 = i1*r31 + i2*r32 + i3*r33

    xxo1 = (m*l/d) + xxi1
    xxo3 = xxi3
    xxo2 = sqrt ( 1.0 - xxo1*xxo1 - xxo3*xxo3 )

    o1 = xxo1*rt11 + xxo2*rt12 + xxo3*rt13
    o2 = xxo1*rt21 + xxo2*rt22 + xxo3*rt23
    o3 = xxo1*rt31 + xxo2*rt32 + xxo3*rt33

    return o1,o2,o3

grating_r0 = 2.0
rowland_r0 = grating_r0/2.0
ga=1.0
gb=1.0
alpha0 = 0.0
rowland_x0 = rowland_r0 * sin(alpha)
rowland_y0 = 0.0
rowland_z0 = grating_z0 + rowland_r0*cos(alpha0)
grating_x0 = grating_r0 * sin(alpha0)
grating_y0 = 0.0
grating_z0 = db - grating_r0*cos(alpha0)

def grating(x,y,grating_z0=grating_z0,grating_r0=grating_r0,ga=ga,gb=gb,grating_x0=grating_x0,grating_y0=grating_y0):
    """
    Defines a spherical mirror.  Default is r=2.0m, positioned 0.1m behind primary
    x^2/a^2 + y^2/b^2 + (z-z0)^2 = r^2
    z = - sqrt(r^2-x^2/a^2-y^2/a^2) + z0
    """

    return -sqrt(grating_r0**2 - (x-grating_x0)**2/ga**2 - (y-grating_y0)**2/gb**2) + grating_z0

def focus_to_grating(x,y,z,xv,yv,zv,ga=ga,gb=gb,grating_r0=grating_r0,grating_z0=grating_z0,grating_x0=grating_x0,grating_y0=grating_y0):
    """
    Returns the time a vector must travel from the focal point to hit the grating
    """

    xi = x - grating_x0
    yi = y - grating_y0

    term1 = (2*gb**2*xi*xv + 2*ga**2*yi*yv + 2*ga**2*gb**2*z*zv - 2*ga**2*gb**2*grating_z0*zv)
    term2 = -1.0* sqrt((-2*gb**2*xi*xv - 2*ga**2*yi*yv - 2*ga**2*gb**2*z*zv + 2*ga**2*gb**2*grating_z0*zv)**2 - 
      4*(ga**2*gb**2*grating_r0**2 - gb**2*xi**2 - ga**2*yi**2 - ga**2*gb**2*z**2 + 
         2*ga**2*gb**2*z*grating_z0 - ga**2*gb**2*grating_z0**2)*(-gb**2*xv**2 - ga**2*yv**2 - 
         ga**2*gb**2*zv**2))
    term3 = 1.0/(2*(-gb**2*xv**2 - ga**2*yv**2 - ga**2*gb**2*zv**2))

    return (term1+term2)*term3
    

def gratingnormal(x,y,grating_r0=grating_r0,ga=ga,gb=gb,grating_x0=grating_x0,grating_y0=grating_y0):
    """
    returns the grating normal.  Test case is xn,yn,zn = 0,0,1
    x^2/a^2 + y^2/b^2 + (z-z0)^2 = r^2
    z = - sqrt(r^2-x^2/a^2-y^2/a^2) + z0
    dzdx = x / a^2 * 1/sqrt(r^2-x^2/a^2-y^2/a^2)
    """

    dzdx = (x-grating_x0)/ga**2 * 1.0/sqrt(grating_r0**2 - (x-grating_x0)**2/ga**2 - (y-grating_y0)**2/gb**2)
    dzdy = (y-grating_y0)/gb**2 * 1.0/sqrt(grating_r0**2 - (x-grating_x0)**2/ga**2 - (y-grating_y0)**2/gb**2)
    xn = -dzdx
    yn = -dzdy
    zn = 1.0 #1.0/dzdx * xn
    vn = array([xn,yn,zn])
    vn /= sqrt( (vn**2).sum() )
    xn,yn,zn=vn

    return vn #xn,yn,zn

def gratingline(x,y,grating_r0=grating_r0,ga=ga,gb=gb):
    """
    returns the grating tangent along a line parallel to the x axis
    """
    dzdx = x/ga**2 * 1.0/sqrt(grating_r0**2 - (x-grating_x0)**2/ga**2 - (y-grating_y0)**2/gb**2)
    xn = 1.0
    yn = 0
    zn = dzdx
    vn = array([xn,yn,zn])
    vn /= sqrt( (vn**2).sum() )
    xn,yn,zn=vn
    return vn #xn,yn,zn

def groovespacing(x,y,d0=1.0/3.6e6,grating_r0=grating_r0):
    """
    computes d(x,y) - the groove spacing projected on the grating
    """
    dz = (grating_r0**2-x**2-(y-d0/2.0)**2)**0.5 - (grating_r0**2-x**2-(y+d0/2.0)**2)**0.5
    newD = sqrt(d0**2 + dz**2)
    return newD

def projectthru_grating(x0,y0,xi=0.0,yi=0.0,zi=1.0,lam=1200e-10,gratingpos=False,final=False,gorder=1):
    x1,y1,z1 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,initial=True)
    x2,y2,z2 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,final=True)
    xv1,yv1,zv1 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,vector=True)
    x3,y3,z3 = secondbounce(x2,y2,z2,xv1,yv1,zv1)
    xv2,yv2,zv2 = secondbounce(x2,y2,z2,xv1,yv1,zv1,vector=True)
    t = focus_to_grating(x3,y3,z3,xv2,yv2,zv2)
    x4 = x3 + t * xv2
    y4 = y3 + t * yv2
    z4 = z3 + t * zv2
    xv3,yv3,zv3 = gratingbounce_jim(x4,y4,z4,xv2,yv2,zv2,lg=lam,dg=groovespacing(x4,y4),mg=gorder)
    if gratingpos: 
        return x4,y4,z4
    elif final:
        return spectral_detector(x4,y4,z4,xv3,yv3,zv3)
    else: 
        return xv3,yv3,zv3

def plotgrating(x,y,ax=None,linlen=0.025):
    if ax is None:
        f=figure(2)
        ax=axes3d.Axes3D(f)
    x1 = linspace(-0.1,0.1,100)
    y1 = x1 * 0
    z1 = grating(x1,y1)
    x2 = linspace(-0.1,0.1,100)
    y2 = x2 * 0 + 0.1
    z2 = grating(x2,y2)
    
    x3d,y3d = ( indices([11,11]) - 5.0 ) / 50.0
    z3d = grating(x3d,y3d)
    ax.plot_wireframe(x3d,y3d,z3d)

    xg,yg,zg = projectthru_grating(x,y,gratingpos=True)
    xv,yv,zv = projectthru_grating(x,y)
    #ax.plot3D(x1,y1,z1)
    #ax.plot3D(x2,y2,z2)
    ax.plot3D([xg,xg+xv*linlen],[yg,yg+yv*linlen],[zg,zg+zv*linlen])
    draw()

    return ax

def spectral_detector(x,y,z,xv,yv,zv,rowland_r0=rowland_r0):
    """
    matches incident vector to point on rowland circle (should be a focal point)
    rowland circle defined by:
    x**2 + y**2 + (z-z0)**2 = rowland_r0**2
    centered at -1.1m (db - rowland_r0)
    x + xv*t = xf
    y + yv*t = yf
    z + zv*t = zf = +/- sqrt(r0**2 - x**2 - y**2) + z0
    (x+xv*t)**2 + (y+yv*t)**2 + (z-z0+zv*t)**2 = r0**2
    (xv**2 + yv**2 + zv**2) * t**2 + (2x+2y+2*(z-z0)) * t + (x**2 + y**2 + (z-z0)**2 - r0)**2 = 0

    Whoops!  want a CYLINDRICAL detector

    x^2+(z-z0)^2 = r0^2
    (x+xv*t)**2 + (z-z0+zv*t)**2 = r0**2
    (xv**2 + zv**2) * t**2 + (2*x*xv + 2*(z-z0)*zv) * t + (x**2 + (z-z0)**2 - r0**2) = 0

    """
    
    z0 = db-rowland_r0
    r0 = rowland_r0

    QA = (yv**2+zv**2)
    QB = (2*y*yv+2*(z-z0)*zv)
    QC = (y**2 + (z-z0)**2 - r0)**2
    t1 = (-QB + sqrt(QB**2 - 4*QA*QC))/(2*QA)
    t2 = (-QB - sqrt(QB**2 - 4*QA*QC))/(2*QA)
    t=t1
    xf = x + t*xv
    yf = y + t*yv
    zf = z + t*zv

    return xf,yf,zf
