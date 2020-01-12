from cassegrain import *
#from numpy import sqrt
#import mathutil

#def sqrt(x):
#    myf = lambda y: y**2 - x
#    return fsolve(myf,numpy.sqrt(x),xtol=1e-15)

grating_r0 = float64(2.0)
rowland_r0 = grating_r0/2.0
#grating_z0 = db#-rowland_r0 #+1.0*grating_r0-db-rowland_r0
#rowland_z0 = db-rowland_r0
#rowland_y0 = float64(0)
rowland_x0 = float64(0)
ga=float64(1.0)
gb=float64(1.0)
gorder = float64(1)
if not globals().has_key('toroid_c'): toroid_c = float64(0)
if not globals().has_key('toroid_pm1'): toroid_pm1 = -1
toroid_pm1 = -1 # it HAS to be, that makes it on the right side of the system....
if not globals().has_key('toroid_pm2'): toroid_pm2 = -1


def setup(a,c=0):
    global alpha
    alpha = float64(a)
    global rowland_y0,rowland_z0,zb,grating_y0,grating_z0,toroid_c
    toroid_c = c
    rowland_y0 = -rowland_r0*sin(alpha*pi/180.)
    rowland_z0 = db - (rowland_r0*cos(alpha*pi/180.))
    zb = db - grating_r0 * cos(alpha*pi/180.) #-2.1
    grating_y0 = -grating_r0*sin(alpha*pi/180.)
    #grating_z0 = db - grating_r0 - toroid_pm1*sqrt(grating_r0**2+toroid_c**2-grating_y0**2+ toroid_pm2*toroid_c*sqrt(grating_r0**2-grating_y0**2))
    grating_z0 = zb + grating_r0*cos(alpha*pi/180.)

# toroid_c = 0.31675 for lambda=1500

#print "alpha: %f   beta(1200): %f   beta(2000): %f" % (alpha,arcsin(1200e-10*3.6e6-sin(alpha*pi/180.))*180/pi,arcsin(2000e-10*3.6e6-sin(alpha*pi/180.))*180/pi)

def do_prob2(alpha=0,toroid_c=0):
    setup(alpha,toroid_c)
    lamarr = arange(1200e-10,2001e-10,100e-10)
    for gorder in [1,2]:
        for lam in [1200e-10,1500e-10,2000e-10]:
            if gorder * lam < 2777e-10:
                focalplane(0,grating=True,lam=lam,nth=200,gorder=gorder)
                focalplane(0,grating=True,lam=lam*(1+1/2e4),nth=200,gorder=gorder)
                savefig("grating_focalplane_restest_alpha%i_gorder%i_c%0.2f_%i.png" % (alpha,gorder,toroid_c,lam*1e10))
                clf()
                focalplane(0,grating=True,lam=lam,nth=200,plotall=True,gorder=gorder)
                focalplane(0,grating=True,lam=lam*(1+1/2e4),nth=200,plotall=True,gorder=gorder)
                savefig("grating_focalplane_restest_alpha%i_plotall_gorder%i_c%0.2f_%i.png" % (alpha,gorder,toroid_c,lam*1e10))
                clf()
    for lam in lamarr:
        focalplane_rays(lam=lam)
        savefig("grating_trace3d_alpha%i_c%0.2f_%i.png" % (alpha,toroid_c,lam*1e10))
        clf()
        plotray(0,0.5,plotx=False,lam=lam)
        plotray(0,-0.5,plotx=False,lam=lam)
        savefig("grating_trace_y_alpha%i_c%0.2f_%i.png" % (alpha,toroid_c,lam*1e10))
        clf()
        focalplane(0,grating=True,lam=lam,nth=200)
        xlabel('mm')
        ylabel('mm')
        savefig("grating_focalplane_alpha%i_c%0.2f_%i.png" % (alpha,toroid_c,lam*1e10))
        clf()


def plotmirrors(plotx=True):
    x = linspace(-0.5,0.5,1000)
    fig=figure(0,figsize=(8,8))
    #fig.add_axes([0.05,0.05,0.95,0.95])
    plot(x,mirror1(x,x*0))
    x2 = linspace(-rad2,rad2,1000)
    plot(x2,mirror2(x2,x2*0))
    x3 = linspace(-rad2*4,rad2*4,1000)
    if plotx: plot(x3,grating(x3,x3*0),label='grating')
    else: plot(x3,grating(x3*0,x3),label='grating')
    #fig.axes[0].set_aspect(1)
    xv,yv,zv = mirror1(-0.5,0,normal=True)
    xv,yv,zv = firstbounce(-0.5,0,vector=True)
    z0 = mirror1(-0.5,0)
    plot([-0.1,0.1],[db,db])
    #plot([-0.5,xv*3-0.5],[z0,z0+zv*3])
    #plot([0.5,-xv*3+0.5],[z0,z0+zv*3])
    if not plotx:
        # plotcircle
        #plot(cos(linspace(0,2*pi,5000))+rowland_y0,sin(linspace(0,2*pi,5000))+rowland_z0,label='rowland')
        plot(linspace(-2,2,5000),rowlandcircle(linspace(-2,2,5000)))
        plot(linspace(-2,2,5000),rowlandcircle(linspace(-2,2,5000),neg=True))
        plot(rowland_y0,rowland_z0,'x')
    fig.axes[0].axis([-2.55,2.55,-2.1,3])
    #fig1=figure(1,figsize=(8,8))
    #clf()
    #plot(x2,mirror2(x2,x2*0))
    #plot([-0.5,xv*3-0.5],[z0,z0+zv*3])
    #plot([0.5,-xv*3+0.5],[z0,z0+zv*3])
    #fig1.axes[0].axis([-0.15,0.15,2.45,2.75])
    return fig

def rowlandcircle(y,rowland_y0=rowland_y0,rowland_z0=rowland_z0,rowland_r0=rowland_r0,neg=False):
    """
    (y-y0)^2 + (z-z0)^2 = rowland_r0^2
    z = sqrt(r0^2 - (y-y0)^2) + z0
    """
    z = sqrt(rowland_r0**2 - (y-rowland_y0)**2) 
    if neg: return -z+rowland_z0
    else:   return z+rowland_z0

def focalplane_rays(ax=None,gorder=1,lam=1200e-10,dospot=False):
    if ax is None:
        f0 = figure(0)
        ax=axes3d.Axes3D(f0)
    xyin = spot()
    #xinc,yinc,zinc = (offset/206265.,0,sqrt(1-(offset/206265.)**2))
    if dospot: 
        for i,xyi in enumerate(array(xyin).T):
            plotray(*xyi,gorder=gorder,plot3d=True,ax=ax,lam=lam)
    else:
        plotray(0.5,0.5,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(-0.5,0.5,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(0.5,-0.5,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(-0.5,-0.5,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(0.5,0,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(-0.5,0,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(0,-0.5,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(0,0.5,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(0.25,0.25,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(-0.25,0.25,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(0.25,-0.25,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(-0.25,-0.25,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(0.25,0,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(-0.25,0,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(0,-0.25,gorder=gorder,plot3d=True,ax=ax,lam=lam)
        plotray(0,0.25,gorder=gorder,plot3d=True,ax=ax,lam=lam)
    return ax

def zoom_detector(ax):
    ax.set_xlim3d(-0.02,0.02)
    ax.set_ylim3d(0.65,0.85)
    ax.set_zlim3d(-0.52,-0.45)
    draw()


def plotray(x0,y0,xi=0.0,yi=0.0,zi=1.0,plotx=True,gorder=1,lam=1200e-10,plot3d=False,f0=None,ax=None,**kwargs):

    x1,y1,z1 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,initial=True)
    x2,y2,z2 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,final=True)
    xv1,yv1,zv1 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,vector=True)
    x3,y3,z3 = secondbounce(x2,y2,z2,xv1,yv1,zv1)
    xv2,yv2,zv2 = secondbounce(x2,y2,z2,xv1,yv1,zv1,vector=True)
    t = focus_to_grating(x3,y3,z3,xv2,yv2,zv2,**kwargs)
    x4 = x3 + t * xv2
    y4 = y3 + t * yv2
    z4 = z3 + t * zv2
    xv3,yv3,zv3 = gratingbounce_jim(x4,y4,z4,xv2,yv2,zv2,lg=lam,mg=gorder,dg=groovespacing(x4,y4,xv2,yv2,zv2,**kwargs),**kwargs)
    #xv3,yv3,zv3 = gratingbounce(x3,y3,z3,xv2,yv2,zv2,lg=lam,mg=gorder)
    xn3,yn3,zn3 = gratingnormal(x4,y4,**kwargs)
    x5,y5,z5 = spectral_detector(x4,y4,z4,xv3,yv3,zv3)
    if f0 is None:
        f0=figure(0,[8,8])
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
        x3d,y3d = ( indices([11,11]) - 5.0 ) / 500.0
        y3d += 0.78
        z3d = detector(x3d,y3d)
        ax.plot_wireframe(x3d,y3d,z3d)
        ax.plot3D([x1,x2,x3,x4,x5],[y1,y2,y3,y4,y5],[z1,z2,z3,z4,z5])
        return ax
    elif plotx: 
        fig=plotmirrors(plotx=plotx)
        plot([x1,x2,x3,x4,x5],[z1,z2,z3,z4,z5],label='bounces')
        plot([x1],[z1],'o')
        plot([x2,x2+xn],[z2,z2-zn])
        plot([x4,x4+xn3*0.5],[z4,z4+zn3*0.5])
        gca().axis([-2.6,2.6,-2.2,3])
        legend(loc='best')
    else:
        fig=plotmirrors(plotx=plotx)
        plot([y1,y2,y3,y4,y5],[z1,z2,z3,z4,z5],label='bounces')
        plot([y1],[z1],'o')
        plot([y2,y2+yn],[z2,z2-zn])
        plot([y4,y4+yn3],[z4,z4+zn3])
        gca().axis([-2.6,2.6,-2.2,3])
        legend(loc='best')
    #gca().axis([-0.15,0.15,2.45,2.75])
    #plot([x0,x0+xv*

def focalplane(offset,fignum=3,grating=False,lam=1200e-10,gorder=1,square=False,nr=10,nth=60,plotall=False,gspace=False,**kwargs):
    figure(fignum)
    xyin = spot(square=square,nr=nr,nth=nth)
    xout,yout,zout = zeros(xyin[0].shape[0]),zeros(xyin[0].shape[0]),zeros(xyin[0].shape[0])
    dg = zeros(xyin[0].shape[0])
    xinc,yinc,zinc = (offset/206265.,0,sqrt(1-(offset/206265.)**2))
    if grating: 
        for i,xyi in enumerate(array(xyin).T):
            if gspace: dg[i] = projectthru_grating(xyi[0],xyi[1],xinc,yinc,zinc,lam=lam,gorder=gorder,final=True,gspace=gspace,**kwargs)
            else: xout[i],yout[i],zout[i] = projectthru_grating(xyi[0],xyi[1],xinc,yinc,zinc,lam=lam,gorder=gorder,final=True,**kwargs)
        #xcen,ycen,zcen = projectthru_grating(0,0,0,0,1,lam=1500e-10,gorder=gorder,final=True)
        if gspace:
            plot(dg,',')
        else:
            if plotall:
                subplot(311)
                plot(yout*1000.0,zout*1000.0,',')
                xlabel('y'); ylabel('z')
                subplot(312)
                plot(xout*1000.0,zout*1000.0,',')
                xlabel('x'); ylabel('z')
                subplot(313)
                plot(xout*1000.0,yout*1000.0,',')
                xlabel('x'); ylabel('y')
            else:
                xo,yo = project_specdec(xout,yout,zout,**kwargs)
                try:
                    plot(xo*1000.0,yo*1000.0,',',label=kwargs['label'])
                except:
                    plot(xo*1000.0,yo*1000.0,',')
                xlabel('mm'); ylabel('mm')
    else:
        for i,xyi in enumerate(array(xyin).T):
            xout[i],yout[i],zout[i] = projectthru(xyi[0],xyi[1],xinc,yinc,zinc,**kwargs)
        plot(xout*1000.0,yout*1000.0,',')

def gratingbounce_jim(x,y,z,xi,yi,zi,vector=False,mg=1,lg=1200e-10,dg=1.0/3.600e6,**kwargs):
#	m = order
#	d = line spacing (angstroms/line)
#	l = wavelength (angstroms)
#	N is plane normal, pointing up on side that rays hit
#	G is vector in direction of lines on grating
#	I is incident vector pointing in ray prop. direction
#	o is output vector pointing in ray prop. direction
#	x indicates primed coordinate system
#	xx indicates double primed coordinate system

    n1,n2,n3 = gratingnormal(x,y,**kwargs)
    g1,g2,g3 = gratingline(x,y,**kwargs)
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
        if negative: 
            pass
            # this is bad!
            #print "error xxi1 is negative no matter what: ",xxi1
            #return 0,0,0
        else:
            return gratingbounce_jim(x,y,z,xi,yi,zi,vector=vector,mg=1,lg=lg,dg=dg,negative=True)

    xxi3 = i1*r31 + i2*r32 + i3*r33

    xxo1 = (m*l/d) + xxi1
    xxo3 = xxi3
    xxo2 = sqrt ( 1.0 - xxo1*xxo1 - xxo3*xxo3 )

    o1 = xxo1*rt11 + xxo2*rt12 + xxo3*rt13
    o2 = xxo1*rt21 + xxo2*rt22 + xxo3*rt23
    o3 = xxo1*rt31 + xxo2*rt32 + xxo3*rt33

    if isnan(o1): pdb.set_trace()
    return o1,o2,o3

def grating(x,y,grating_z0=grating_z0,grating_r0=grating_r0,toroid_c=toroid_c,ga=ga,gb=gb,grating_y0=grating_y0):
    """
    Defines a spherical mirror.  Default is r=2.0m, positioned 0.1m behind primary
    x^2/a^2 + y^2/b^2 + (z-z0)^2 = r^2
    z = - sqrt(r^2-x^2/a^2-y^2/a^2) + z0

    Toroid:
    (c-(x**2+y**2)**0.5)**2 + z**2 = a**2   (symmetric about z-axis)
    we want:
    (c-(x**2+z**2)**0.5)**2 + y**2 = a**2 or x/y switched
    c**2 - 2 * c * sqrt(x**2 + z**2) + x**2 + y**2 + z**2 = a**2

    z = +/- sqrt( a^2 + c^2 - x^2 - y^2 +/- 2 sqrt( a^2 c^2 - c^2 y^2 ))

    need m020 = 1/r + 1/r' - 2 a02 (cos alpha + cos beta) = 0
    r, r' are chords from grating to detector
    a02 is the radius of curvature in one dimension

    m200 = cos^2 alpha / r + cos^2 beta / r' - 2 a20 (cos alpha + cos beta) 
    may also have to be zero
    """
    return toroid_pm1*sqrt( grating_r0**2 + toroid_c**2 - x**2 - (y-grating_y0)**2
            + toroid_pm2*2*toroid_c*sqrt(grating_r0**2-(y-grating_y0)**2) ) + grating_z0 + toroid_pm2*toroid_c

def focus_to_grating(x,y,z,xv,yv,zv,ga=ga,gb=gb,grating_r0=grating_r0,grating_z0=grating_z0,toroid_c=toroid_c,grating_y0=grating_y0,**kwargs):
    """
    Returns the time a vector must travel from the focal point to hit the grating
    Done by solving xfoc+xv*t = xg for all 3 positional arguments and the equation
    of the (spherical or toroidal) grating
    """

    x = float64(x)
    y = float64(y)
    xv = float64(xv)
    yv = float64(yv)
    zv = float64(zv)
    toroid_c = float64(toroid_c)

    def fz(t):
        term1 = (toroid_c + toroid_pm2*sqrt( (x+xv*t)**2 + (z+zv*t-grating_z0+toroid_pm2*toroid_c)**2 ) )**2 
        term2 = (y+yv*t-grating_y0)**2 - grating_r0**2
        return float64(term1)+float64(term2)
    def afz(t):
        return abs(fz(t))
        #return  (z+zv*t) - sqrt( grating_r0**2 + toroid_c**2 - (x+xv*t)**2 - (y+yv*t)**2 - 2*sqrt( grating_r0**2*toroid_c**2 - toroid_c**2 * (y+yv*t)**2 ))
    #t1 = fsolve(fz,1.0,warning=False,xtol=1e-35,maxfev=3000)
    t = bisect(fz,0,3,xtol=float64(1e-35),rtol=float64(1e-35),maxiter=3000)
    #t = findroot(fz,float64(t1),tol=float64(1e-18),solver='bisect',maxsteps=10)
    #t = mathutil.root(fz,(float64(0.5),float64(2.5)),accuracy=1e-35)
    #print "rootfind: %30.25f, %30.25f, t-t1: %30.25f" % (t1, t, t-t1)
    #term1 = (2*gb**2*x*xv + 2*ga**2*y*yv + 2*ga**2*gb**2*z*zv - 2*ga**2*gb**2*grating_z0*zv)
    #term2 = -1.0* sqrt((-2*gb**2*x*xv - 2*ga**2*y*yv - 2*ga**2*gb**2*z*zv + 2*ga**2*gb**2*grating_z0*zv)**2 - 
    #  4*(ga**2*gb**2*grating_r0**2 - gb**2*x**2 - ga**2*y**2 - ga**2*gb**2*z**2 + 
    #     2*ga**2*gb**2*z*grating_z0 - ga**2*gb**2*grating_z0**2)*(-gb**2*xv**2 - ga**2*yv**2 - 
    #     ga**2*gb**2*zv**2))
    #term3 = 1.0/(2*(-gb**2*xv**2 - ga**2*yv**2 - ga**2*gb**2*zv**2))
    #t = (term1+term2)*term3
    #t2 = brenth(fz,0,3,xtol=float64(1e-35),rtol=float64(1e-35),maxiter=3000)
    #t3 = fmin(afz,1.5,xtol=float64(1e-35))
    #print "fsolve: %30.25g  bisect: %30.25g  brenth: %30.25g " % (t,t1,t2)
    #print "f-bisect: %30.25g   f-brent: %30.25g  bisect-brent:%30.25g" % (t-t1,t-t2,t1-t2)

    return t
    

def gratingnormal(x,y,grating_r0=grating_r0,ga=ga,gb=gb,tangent=False,toroid_c=toroid_c,grating_y0=grating_y0,**kwargs):
    """
    returns the grating normal.  Test case is xn,yn,zn = 0,0,1
    x^2/a^2 + y^2/b^2 + (z-z0)^2 = r^2
    z = - sqrt(r^2-x^2/a^2-y^2/a^2) + z0
    dzdx = x / a^2 * 1/sqrt(r^2-x^2/a^2-y^2/a^2)

    z = +/- sqrt( grating_r0^2 + toroid_c^2 - x^2 - y^2 +/- 2 sqrt( grating_r0^2 toroid_c^2 - toroid_c^2 y^2 ))
    """

    denom = sqrt(grating_r0**2 + toroid_c**2 - x**2 - (y-grating_y0)**2 + toroid_pm2*2*sqrt(grating_r0**2 - (y-grating_y0)**2)*toroid_c )
    dzdx = toroid_pm1*(-1)*x/denom
    dzdy = toroid_pm1*(-1)*(y-grating_y0)*(1+toroid_pm2*toroid_c/sqrt(grating_r0**2-(y-grating_y0)**2)) / denom

    # analytic, non-parametric normal
    xn = -dzdx
    yn = -dzdy
    zn = float64(1.0) #1.0/dzdx * xn
    vn = array([xn,yn,zn])
    vn /= sqrt( (vn**2).sum() )
    xn,yn,zn=vn
    if tangent:
        zt = sqrt(yn**2/(yn**2 + zn**2))
        yt = sqrt(1-zt**2)
        xt = 0
        vt = array([xt,yt,zt])
        #vt = cross(vn,gratingline(x,y,**kwargs))
        vt /= sqrt( (vt**2).sum() )
        return -vt

    if isnan(xn): pdb.set_trace()
    return vn #xn,yn,zn

def gratingline(x,y,grating_r0=grating_r0,ga=ga,gb=gb,toroid_c=toroid_c,negative=False,grating_y0=grating_y0,**kwargs):
    """
    returns the grating tangent along a line parallel to the x axis
    note that I define the grating lines to be along the x-axis,
    which means that the other tangent vector can have components in
    all 3 directions
    """
    # parametric definition of the torus
    #sinnu = ((y-grating_y0)/grating_r0)
    #cosnu = sqrt(float64(1) - sinnu**2)
    #cosmu = ( x / (toroid_c + grating_r0*cosnu)  )
    #sinmu = sqrt(float64(1) - cosmu**2)
    # wrong "tangent" xn = dxdnu = -1*grating_r0*sinnu * cosmu
    # wrong "tangent" yn = dydnu = -grating_r0 * cosnu
    # wrong "tangent" zn = dzdnu = -1*grating_r0*sinnu * sinmu
    #xn = -(toroid_c+grating_r0*cosnu)*-1*sinmu
    #yn = 0
    #zn = (toroid_c+grating_r0*cosnu)*cosmu
    #vn = array([xn,yn,zn])
    vn = cross(gratingnormal(x,y,grating_r0=grating_r0,toroid_c=toroid_c,grating_y0=grating_y0),
            gratingnormal(x,y,grating_r0=grating_r0,toroid_c=toroid_c,grating_y0=grating_y0,tangent=True))
    vn /= sqrt( (vn**2).sum() )
    xn,yn,zn=vn
    if isnan(xn): pdb.set_trace()
    if negative: return -vn
    else:        return vn

def groovespacing(x,y,xv,yv,zv,d0=1.0/3.6e6,grating_r0=grating_r0,corr_in=False,corr_flat=True,**kwargs):
    """
    computes d(x,y) - the groove spacing projected on the grating
    """
    # these three lines are old and now kept for debug purposes only
    #dz = (grating_r0**2-x**2-(y-d0/2.0)**2)**0.5 - (grating_r0**2-x**2-(y+d0/2.0)**2)**0.5
    dz = grating(x,y-d0/2.0) - grating(x,y+d0/2.0)
    oldD = sqrt(d0**2 + dz**2)
    v1 = gratingnormal(x,y,tangent=True,**kwargs)
    v2 = gratingnormal(0,0,tangent=True,**kwargs)
    cosbeta = dot(v1,v2)
    newD = d0 / abs( cosbeta )
    if corr_in:
        # apparently this is already done in Jim's code # need to project into incoming-light angle
        gn = gratingnormal(x,y)
        iv = array([xv,yv,zv])
        dgi = abs(dot(gn,iv))  # this is exactly 1 for the alpha=0 case
        if corr_flat: return newD*dgi
        else:  return d0*dgi 
    elif corr_flat: 
        return newD 
    else: return d0

def projectthru_grating(x0,y0,xi=0.0,yi=0.0,zi=1.0,lam=1200e-10,gratingpos=False,final=False,gorder=1,gspace=False,**kwargs):
    x1,y1,z1 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,initial=True)
    x2,y2,z2 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,final=True)
    xv1,yv1,zv1 = firstbounce(x0,y0,xi=xi,yi=yi,zi=zi,vector=True)
    x3,y3,z3 = secondbounce(x2,y2,z2,xv1,yv1,zv1)
    xv2,yv2,zv2 = secondbounce(x2,y2,z2,xv1,yv1,zv1,vector=True)
    t = focus_to_grating(x3,y3,z3,xv2,yv2,zv2,**kwargs)
    x4 = x3 + t * xv2
    y4 = y3 + t * yv2
    z4 = z3 + t * zv2
    #xv3,yv3,zv3 = gratingbounce_jim(x4,y4,z4,xv2,yv2,zv2,lg=lam,dg=1.0/3.6e6,mg=gorder)
    dg = groovespacing(x4,y4,xv2,yv2,zv2,**kwargs)
    if gspace: return dg
    #debug print x4,y4,xv2,yv2,zv2,dg
    xv3,yv3,zv3 = gratingbounce_jim(x4,y4,z4,xv2,yv2,zv2,lg=lam,dg=dg,mg=gorder,**kwargs)
    if gratingpos: 
        if isnan(x4): pdb.set_trace()
        return x4,y4,z4
    elif final:
        if z4 > -0.1: pdb.set_trace()
        return spectral_detector(x4,y4,z4,xv3,yv3,zv3)
    else: 
        if isnan(xv3): pdb.set_trace()
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

def detector(x,y):
    z0 = db-rowland_r0
    z = sqrt(rowland_r0**2 - y**2 ) + z0
    return z

def spectral_detector(x,y,z,xv,yv,zv,rowland_r0=rowland_r0):
    """
    matches incident vector to point on rowland circle (should be a focal point)
    rowland circle defined by:
    y**2 + (z-z0)**2 = rowland_r0**2
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
    
    r0 = rowland_r0

    QA = float64(  (yv**2+zv**2)  )
    QB = float64(  (float64(2)*(y-rowland_y0)*yv+float64(2)*(z-rowland_z0)*zv)  )
    QC = float64(  (y-rowland_y0)**2 + (z-rowland_z0)**2 - r0**2  )
    """
    QA = (xv**2+zv**2)
    QB = (2*x*xv+2*(z-rowland_z0)*zv)
    QC = (x**2 + (z-rowland_z0)**2 - r0)**2
    """
    t1 = (-QB + sqrt(QB**2 - 4*QA*QC))/(2*QA)
    t2 = (-QB - sqrt(QB**2 - 4*QA*QC))/(2*QA)
    t=t1
    #pdb.set_trace();
    #def fz(t):
    #    return (y+yv*t-rowland_y0)**2 + (z+zv*t-rowland_z0)**2 - rowland_r0**2
    #t = mathutil.root(fz,(float64(1.5),float64(2.5)),accuracy=1e-25,max_iterations=1000)
    #t = fsolve(fz,2,xtol=float64(1e-35),warning=False)
    #t3 = bisect(fz,1.5,2.5,xtol=float64(1e-35),rtol=float64(1e-35),maxiter=5000)
    #t3 = bisect(fz,t3-.001,t3+.001,xtol=float64(1e-35),rtol=float64(1e-35),maxiter=5000)
    #t3 = findroot(fz,float64(t1),tol=float64(1e-18),solver='bisect',maxsteps=10)
    #print "t: %30.25f , t3: %30.25f, t-t3: %30.25f,  fz: %30.25f,%30.25f" % (t,t3,t-t3,fz(t),fz(t3))
    #print t1,t2,fz(t1),fz(t2)

    xf = x + t*xv
    yf = y + t*yv
    zf = z + t*zv

    #if isnan(xf): pdb.set_trace()
    return xf,yf,zf


# center of detector in x,y,z space
xcen,ycen,zcen = projectthru_grating(0,0,0,0,1,lam=float64(1500e-10),gorder=gorder,final=True)

def project_specdec(x,y,z,**kwargs):

    # geometry - how far along the cylinder away from its (arbitrary) center?
    dy = (y-rowland_y0)
    dz = (z-rowland_z0)
    dx = (x-rowland_x0)
    th1 = ( arccos(dy/rowland_r0) - arccos((ycen-rowland_y0)/rowland_r0)  )
    th2 = ( arcsin(dz/rowland_r0) - arcsin((zcen-rowland_z0)/rowland_r0)  )  # numerical error
    th3 = ( arctan2(dz,dy) - arctan2((zcen-rowland_z0),(ycen-rowland_y0)) )
    arcdist = rowland_r0 * th1
    #print (th1[0:5],th2[:5],th3[:5],dx[:5],dy[:5],dz[:5],(dx**2+dz**2)[:5])  

    yo = arcdist
    #yo = y-ycen
    xo = x-xcen
    #xo = y-ycen
    #yo = x-xcen
    #yo = y-ycen
    #zo = z-zcen
    #xo = x 
    #yo = sqrt((y-ycen)**2+arcdist**2)
    #yo = z - zcen

    # x+proj
    #xo = x - xcen
    #yo = arcdist

    return xo,yo

