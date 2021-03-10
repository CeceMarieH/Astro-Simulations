from vpython import *
#GlowScript 2.7 VPython
"""
ioorbit.py
simulates orbit of Jupiter's moons
"""

#define constants
G = 6.67E-11 #gravitational constant
AU = 1.5E11 #one astronomical unit (io-jup distance)
YEAR = 365.25*24*3600 #one year, in seconds
DAY = 24*3600
C = 6.1E-13
Q = 100
Y = 6.5E11 #Young's Modulus, dyne/cm^2 (g*cm/(s^2*cm^2))
Y = 6.5E10 #Young's Modulus, N/m^2

#create objects for jup and io
jup = sphere(pos = vec(0,0,0), mass = 1.898E27, radius = 1E8, color = color.orange)
io = sphere(pos = vec(-16102634.7828827,-423475341.974609,-15286319.908124), mass = 8.933E22, radius = 5E7, color = color.yellow)
#io = box(pos = vec(-16102634.7828827,-423475341.974609,-15286319.908124), axis=vec(0,0,1E8), size=vec(1E8,1E8,1E8), mass = 8.933E22, color = color.yellow)
eu = sphere(pos = vec(-620616074.114273,242577301.122559,4278514.38963509), mass = 4.797E22, radius = 5E7, color = color.blue)
gan = sphere(pos = vec(1067920461.08904,49048414.3837891,14978820.9460497), mass = 1.4819E23, radius = 5E7, color = color.white)
cal = sphere(pos = vec(93610714.4340668,-1878370620.62354,-59448290.5347843), mass = 1.076E23, radius = 5E7, color = color.green)

#initial conditions
t = 0
counter = 0
rmin = mag(io.pos - jup.pos)
rmax = mag(io.pos - jup.pos)
rmin2 = mag(eu.pos - jup.pos)
rmax2 = mag(eu.pos - jup.pos)
rmin3 = mag(gan.pos - jup.pos)
rmax3 = mag(gan.pos - jup.pos)

io.vel = vec(17231.2106442922,-583.810327162546,229.56352413948)
io.trail = curve(pos = io.pos, radius = 0.2*io.radius, color = io.color)

eu.vel = vec(-4934.54951488675,-12921.7664690181,-461.78766730855)
eu.trail = curve(pos = eu.pos, radius = 0.2*eu.radius, color = eu.color)

gan.vel = vec(-488.433498943854,10874.3256200888,399.654407591813)
gan.trail = curve(pos = gan.pos, radius = 0.2*gan.radius, color = gan.color)

cal.vel = vec(8195.61933756686,465.386437689664,122.821035464139)
cal.trail = curve(pos = cal.pos, radius = 0.2*cal.radius, color = cal.color)


#define time step
h = 5E2
scene.autoscale = 1

moon = io
a = eu
b = gan
c = cal
#define acceleration as function
def acc(iopos): #iopos is dummy variable
    r = mag(iopos - jup.pos)
    ra = mag(iopos - a.pos)
    rb = mag(iopos - b.pos)
    rc = mag(iopos - c.pos)
    return -G*jup.mass*(iopos - jup.pos)/r**3 + -G*a.mass*(iopos - a.pos)/ra**3 + -G*b.mass*(iopos - b.pos)/rb**3 + -G*c.mass*(iopos - c.pos)/rc**3
    
def dn_dt(n):
    return -C*n*(1 - 14*ecc**2)
    

"""
Runge-Kutta Method
k_1 = h*f(t_n,y_n)
k_2 = h*f(t_n+h/2,y_n+k_1/2)
k_3 = h*f(t_n+h/2,y_n+k_2/2)
k_4 = h*f(t_n+h,y_n+k_3)
y_(n+1) = Y_n + 1/6(k_1 + 2*k_2 + 2*k_3 + k_4) + O[h^5]
"""
        
    
#define RK4 algorithm
def rk4(moon):
    k1v = acc(moon.pos)*h
    k1x = moon.vel*h
    
    k2v = acc(moon.pos + k1x/2.0)*h
    k2x = (moon.vel + k2v/2.0)*h
    
    k3v = acc(moon.pos + k2x/2.0)*h
    k3x = (moon.vel + k2v/2.0)*h
    
    k4v = acc(moon.pos + k3x)*h
    k4x = (moon.vel + k3v)*h
    
    moon.vel += (k1v + 2*k2v + 2*k3v + k4v)/6.0
    moon.pos += (k1x + 2*k2x + 2*k3x + k4x)/6.0

def rkn(n):
    k1n = dn_dt(n)*h
    k2n = dn_dt(n + k1n/2)*h
    k3n = dn_dt(n + k2n/2)*h
    k4n = dn_dt(n+k3n/2)*h
    n += (k1n + 2*k2n + 2*k3n + k4n)/6.0

#set up plotting
plot = gdisplay(x=0,y=400, height=400, width=400, title="Laplace Check", xtitle = "time (d)", ytitle = "L")
data = gdots(color = color.blue)
data2 = gdots(color = color.red)
    
while True: #loop through calcs to animate
    
    #draw path of io
    moon = io
    a = eu
    b = gan
    c = cal
    rk4(io)
    io.trail.append(pos=io.pos,color=io.color)
    
    #draw path of Europa
    moon = eu
    a = io
    b = gan
    c = cal
    rk4(eu)
    eu.trail.append(pos = eu.pos, color = eu.color)
    
    #calculate mean orbital motion n2 of Europa
    r2 = mag(eu.pos - jup.pos)
    if r2>rmax2:
        rmax2 = r2
    if r2<rmin2:
        rmin2 = r2
    #semi-major axis
    sma2 = (rmax2 + rmin2)/2
    
    #period (using Kepler's Third Law)
    per2 = sqrt(4*pi**2*sma2**3/(G*(eu.mass+jup.mass)))/DAY
    
    n2 = 2*pi/per2
    
    #draw path of Ganymede
    moon = gan
    a = io
    b = eu
    c = cal
    rk4(gan)
    gan.trail.append(pos = gan.pos, color = gan.color)
    
    #calculate mean orbital motion n3 of Ganymede
    r3 = mag(gan.pos - jup.pos)
    if r3>rmax3:
        rmax3 = r3
    if r3<rmin3:
        rmin3 = r3
    #semi-major axis
    sma3 = (rmax3 + rmin3)/2
    
    #period (using Kepler's Third Law)
    per3 = sqrt(4*pi**2*sma3**3/(G*(gan.mass+jup.mass)))/DAY
    
    n3 = 2*pi/per3
    
    #draw path of Callisto
    moon = cal
    a = io
    b = eu
    c = gan
    rk4(cal)
    cal.trail.append(pos = cal.pos, color = cal.color)

    #various physical quantities of Io
    
    r = mag(io.pos - jup.pos)
    if r>rmax:
        rmax = r
    if r<rmin:
        rmin = r
    #semi-major axis
    sma = (rmax + rmin)/2
    
    #period (using Kepler's Third Law)
    per = sqrt(4*pi**2*sma**3/(G*(io.mass+jup.mass)))/DAY
    
    #eccentricity
    ecc = (rmax - rmin)/(rmax + rmin)
    
    #mean orbital motion
    if counter < 500:
        n1 = 2*pi/per
    else:
        rkn(n1)

    #tidal force
    Rio = 3.6426E6 #radius of io
    tf = G*jup.mass*io.mass*((r - Rio)**(-2) - (r+Rio)**(-2))
    
    #rotation of io
    omega = mag(io.vel)/sma*h
    io.rotate(angle = omega, axis = vec(0,0,1), origin = io.pos)
    
    #energy buildup in Io
    E = -27*(G*jup.mass*io.mass)**2*r*ecc/(16*Y*sma**6)
    W = 27*(G*jup.mass*io.mass)**2*pi*r*ecc/(8*Q*per*Y*sma**6)   
    
    #test plot:
    data.plot(t/DAY,L)
    
    
    #test values with Laplace relation: n1 - 3*n2 + 2*n3 = 0
    L = n1 - 3*n2 + 2*n3
    if counter == 1000:
        
        print("Galilean Resonance Laplace Relation: n1 - 3*n2 + 2*n3 = 0")
        print("n1 = {0:1.3f}, n2 = {1:1.3f}, n3 = {2:1.3f}".format(n1,n2,n3))
        print("L = {0:1.3f}. L should be 0.\n".format(L))
        print("Tidal force = {0:1.3f}".format(tf))
        print("Energy buildup = {0:1.3f}".format(W))

        counter = 0
    
    
    t+=h
    counter+=1
    rate(1000)



#end script