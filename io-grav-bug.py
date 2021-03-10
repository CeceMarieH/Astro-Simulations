from vpython import *
#GlowScript 2.7 VPython
"""
ioorbit.py
simulates orbit of io or other planet around star
"""

#define constants
G = 6.67E-11 #gravitational constant
AU = 1.5E11 #one astronomical unit (io-jup distance)
YEAR = 365.25*24*3600 #one year, in seconds
DAY = 24*3600

#create objects for jup and io
jup = sphere(pos = vec(0,0,0), mass = 1.898E27, radius = 5E7, color = color.orange)
io = sphere(pos = vec(4.217E8,0,0), mass = 8.9319E22, radius = 2E7, color = color.yellow)
eu = sphere(pos = vec(6.71034E8,0,0), mass = 4.7998E22, radius = 2E7, color = color.blue)
gan = sphere(pos = vec(1.070412E9,0,0), mass = 1.4819E23, radius = 2E7, color = color.white)

#initial conditions
t = 0
L = vec(0,0,0)
rmax = 4.217E8
rmin = 0
counter = 0

r = mag(io.pos - jup.pos)
v = sqrt(G*jup.mass/r*(1 - 0.0043**2))
io.vel = vec(0,v,0)
io.trail = curve(pos = io.pos, radius = 0.2*2E7, color = io.color)

r = mag(eu.pos - jup.pos)
v = sqrt(G*jup.mass/r*(1 - 0.0094**2))
eu.vel = vec(0,v,0)
eu.trail = curve(pos = eu.pos, radius = 0.2*eu.radius, color = eu.color)

r = mag(gan.pos - jup.pos)
v = sqrt(G*jup.mass/r*(1 - 0.0011**2))
gan.vel = vec(0,v,0)
gan.trail = curve(pos = gan.pos, radius = 0.2*gan.radius, color = gan.color)


#define time step
h = 5E3
scene.autoscale = 1


#define acceleration as function
def acc(iopos): #iopos is dummy variable
    r = mag(iopos - jup.pos)
    ra = mag(iopos - a.pos)
    rb = mag(iopos - b.pos)
    return -G*jup.mass*(iopos - jup.pos)/r**3 + -G*a.mass*(iopos - a.pos)/ra**3 + -G*b.mass*(iopos - b.pos)/rb**3

    
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


#set up plotting
plot = gdisplay(x=0,y=400, height=400, width=400, title="Gravitational Potential of Io", xtitle = "Time (s)", ytitle = "Potential Energy (J)")
data = gdots(color = color.blue)
    
    
while True: #loop through calcs to animate
    
    #draw paths of each moon
    
    #draw path of io
    moon = io
    a = eu
    b = gan
    rk4(io)
    io.trail.append(pos=io.pos,color=io.color)
    
    #draw path of Europa
    moon = eu
    a = io
    b = gan
    rk4(eu)
    eu.trail.append(pos = eu.pos, color = eu.color)
    
    #draw path of Ganymede
    moon = gan
    a = io
    b = eu
    rk4(gan)
    gan.trail.append(pos = gan.pos, color = gan.color)
    
    
    #various physical quantities of Io
    
    #find major axis of orbit
    r = mag(io.pos - jup.pos)
    re = mag(io.pos - eu.pos)
    rg = mag(io.pos - gan.pos)
    if r>rmax:
        rmax=r
        tmax = t
    elif r<rmin:
        rmin=r
        tmin = t
    
    L=io.mass*cross(io.pos,io.vel)
    
    E = 0.5*io.mass*dot(io.vel,io.vel) - G*jup.mass*io.mass/r - G*eu.mass*io.mass/re - G*gan.mass*io.mass/rg
    
    ecc = sqrt(sqrt(abs(2*dot(L,L)*E/(G*(jup.mass+io.mass)**2*io.mass**3))))
    
    a = (rmax-rmin)/2.0 #semi-major axis
    
    period = sqrt(2.0*pi*a**3/(G*(jup.mass+io.mass)))/DAY
    
    
    #check Laplace relation: n1 - 3*n2 - 2*n3 = 0
    n1 = 2*pi/period
    
    p2 = sqrt(2.0*pi*mag(eu.pos - jup.pos)**3/(G*(jup.mass+eu.mass)))/DAY
    n2 = 2*pi/p2
    
    p3 = sqrt(2.0*pi*mag(gan.pos - jup.pos)**3/(G*(jup.mass+gan.mass)))/DAY
    n3 = 2*pi/p3
    
    Laplace = n1 - 3*n2 - 2*n3
    
    data.plot(t,E)
    
    #print every 1000 calculations
    if counter>=1000:
        print("Laplace check value: {0:8.3f} - should be 0.".format(Laplace))
        print("Mag. of Ang. Mom = {0:8.3E}, Energy = {1:8.3E}, Period = {2:8.3f} days, Eccentricity = {3:8.3f}".format(mag(L), E, period, ecc))
        counter = 0
    counter += 1

    
    t+=h
    
    rate(200)
    
    





