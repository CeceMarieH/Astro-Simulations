from vpython import *
#GlowScript 2.7 VPython
"""
earthorbit.py
simulates orbit of earth or other planet around star
"""

#keyboard functions for changing Earth's velocity
def keyInput(evt):
    s = evt.key
    if s == "b":
        earth.vel = earth.vel * 1.5
        print(earth.vel)
    elif s == "f":
        earth.vel = earth.vel/2.0
        print(earth.vel)

scene.bind("keydown",keyInput)
def keyIRel(evt):
    s = evt.key
    if s == "keyup":
        earth.vel = earth.vel

scene.bind("keyup",keyIRel)
#define constants
G = 6.67E-11 #gravitational constant
AU = 1.5E11 #one astronomical unit (earth-sun distance)
YEAR = 365.25*24*3600 #one year, in seconds

#create objects for sun and earth
sun = sphere(pos = vec(0,0,0), mass = 2E30, radius = 1E10, color = color.yellow)
earth = sphere(pos = vec(AU,0,0), mass = 6E24, radius = 2E9, color = color.blue)

#initial conditions
earth.vel = vec(0.5*sqrt(G*sun.mass/AU),sqrt(G*sun.mass/AU),0)
earth.vel = vec(0,2*pi*AU/YEAR,0)
earth.acc = vec(-G*sun.mass/AU**2,0,0)
earth.trail = curve(pos = earth.pos, color = earth.color)
counter = 0 #for printing values
L = vec(0,0,0) #angular momentum of earth
rmin = AU
rmax = 0
t = 0

#define time step
h = 1E5
scene.autoscale = 1


#define acceleration as function
def acc(earthpos): #earthpos is dummy variable
    r = mag(earthpos - sun.pos)
    return -G*sun.mass*(earthpos - sun.pos)/r**3
    

"""
Runge-Kutta Method
k_1 = h*f(t_n,y_n)
k_2 = h*f(t_n+h/2,y_n+k_1/2)
k_3 = h*f(t_n+h/2,y_n+k_2/2)
k_4 = h*f(t_n+h,y_n+k_3)
y_(n+1) = Y_n + 1/6(k_1 + 2*k_2 + 2*k_3 + k_4) + O[h^5]
"""
        
    
#define RK4 algorithm
def rk4(earth):
    k1v = acc(earth.pos)*h
    k1x = earth.vel*h
    
    k2v = acc(earth.pos + k1x/2.0)*h
    k2x = (earth.vel + k2v/2.0)*h
    
    k3v = acc(earth.pos + k2x/2.0)*h
    k3x = (earth.vel + k2v/2.0)*h
    
    k4v = acc(earth.pos + k3x)*h
    k4x = (earth.vel + k3v)*h
    
    earth.vel += (k1v + 2*k2v + 2*k3v + k4v)/6.0
    earth.pos += (k1x + 2*k2x + 2*k3x + k4x)/6.0

#set up plotting
plot = gdisplay(x=0,y=400, height=400, width=400, title="Test Plot", xtitle = "t", ytitle = "test value")
data = gdots(color = color.blue)
data2 = gdots(color = color.red)

while True: #loop through calcs to animate

    #find major axis of orbit
    r = mag(earth.pos - sun.pos)
    if r>rmax:
        rmax=r
    elif r<rmin:
        rmin=r
    
    #various physical quantities
    L=earth.mass*cross(earth.pos,earth.vel)
    E = 0.5*earth.mass*dot(earth.vel,earth.vel) - G*sun.mass*earth.mass/r
    ecc = sqrt(1.0 + (2*dot(L,L)*E)/((G*sun.mass)**2*(earth.mass)**3))
    a = (rmax+rmin)/2.0 #semi-major axis
    ecc = (rmax - rmin)/(rmax + rmin)*2
    period = sqrt(2.0*pi*a**3/(G*(sun.mass+earth.mass)))/YEAR
    
    data.plot(t,ecc)
    
    #draw path of Earth
    rk4(earth)
    earth.trail.append(pos=earth.pos,color=earth.color)
    #print every 1000 calculations
    if counter>=1000:
        print("Mag. of Ang. Mom = {0:8.3E}, Energy = {1:8.3E}, Period = {2:8.3f}, Eccentricity = {3:8.3f}".format(mag(L), E, period, ecc))
        counter = 0
    counter += 1
    t+=h
    rate(200)
    
    
    
    
    
    
    
    