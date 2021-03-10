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

#create objects for jup and io
jup = sphere(pos = vec(0,0,0), mass = 1.898E27, radius = 5E7, color = color.orange)
io = sphere(pos = vec(4.217E8,0,0), mass = 8.9319E22, radius = 2E7, color = color.yellow)
eu = sphere(pos = vec(6.71034E8,0,0), mass = 4.7998E22, radius = 2E7, color = color.blue)
gan = sphere(pos = vec(1.070412E9,0,0), mass = 1.4819E23, radius = 2E7, color = color.white)
cal = sphere(pos = vec(1.882709E9,0,0), mass = 1.079E23, radius = 2E7, color = color.green)

#initial conditions
t = 0

io.vel = vec(0,17334.74,0)
io.trail = curve(pos = io.pos, radius = 0.2*io.radius, color = io.color)

eu.vel = vec(0,1.374E4,0)
eu.trail = curve(pos = eu.pos, radius = 0.2*eu.radius, color = eu.color)

gan.vel = vec(0,1.088E4,0)
gan.trail = curve(pos = gan.pos, radius = 0.2*gan.radius, color = gan.color)

cal.vel = vec(0,8.204E3,0)
cal.trail = curve(pos = cal.pos, radius = 0.2*cal.radius, color = cal.color)


#define time step
h = 1E3
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
    
    #draw path of Ganymede
    moon = gan
    a = io
    b = eu
    c = cal
    rk4(gan)
    gan.trail.append(pos = gan.pos, color = gan.color)
    
    #draw path of Callisto
    moon = cal
    a = io
    b = eu
    c = gan
    rk4(cal)
    cal.trail.append(pos = cal.pos, color = cal.color)
    
    rate(200)





