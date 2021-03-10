from vpython import *
#GlowScript 2.7 VPython
'''
6-1-a.py
simulate a binary star system
'''

#define constants
G = 6.67E-11 #gravitational constant
AU = 1.5E11 #one astronomical unit (earth-sun distance)
YEAR = 365.25*24*3600 #one year, in seconds


#create objects for stars a and b
astar = sphere(pos = vec(-5*AU,0,0), mass = 4E30, radius = 1E10, color = color.yellow)
bstar = sphere(pos = vec(5*AU,0,0), mass = 2E30, radius = 1E10, color = color.red)

com = (astar.mass*astar.pos + bstar.mass*bstar.pos)/(astar.mass+bstar.mass)


#initial conditions
astar.vel = vec(0,2*2.3E3,0)
bstar.vel = vec(0,2*(-4.6E3),0)

astar.trail = curve(pos = astar.pos, color = astar.color)
bstar.trail = curve(pos = bstar.pos, color = bstar.color)

counter = 0 #for printing values
a_r = 0
b_r = 0
a_rmin = 0
a_rmax = 0
b_rmin = 0
b_rmax = 0
a_tmax = 0
a_tmin = 0
b_tmax = 0
b_tmin = 0
a_axis = 0
a_period = 0
a_ecc = 0
b_axis = 0
b_period = 0
b_ecc = 0
t = 0

#define time step
h = 1E6
scene.autoscale = 1

while True:
    
    a_r = mag(astar.pos - com)
    b_r = mag(bstar.pos - com)
    if a_r<a_rmin:
        a_rmin = a_r
        a_tmin = t
    if a_r>a_rmax:
        a_rmax = a_r
        a_tmax = t
    if b_r<b_rmin:
        b_rmin = b_r
        b_tmin = t
    if b_r>b_rmax:
        b_rmax = b_r
        b_tmax = t
        
    
    F = -G*astar.mass*bstar.mass*(astar.pos - bstar.pos)/mag(astar.pos - bstar.pos)**3 #force on astar
    #force on bstar = -force astar

    #implement Euler
    astar.vel += F/astar.mass*h    
    bstar.vel -= F/bstar.mass*h
    astar.pos += astar.vel*h
    bstar.pos += bstar.vel*h
    
    astar.trail.append(pos = astar.pos)
    bstar.trail.append(pos = bstar.pos)
    
    #print semi-major axis, period, and eccentricity of stars
    if counter>=500:
        
        a_axis = (a_rmax - a_rmin)/2
        a_period = abs(a_tmax - a_tmin)*2
        a_ecc = (a_rmax - a_axis)/a_axis
        
        b_axis = (b_rmax - b_rmin)/2
        b_period = abs(b_tmax - b_tmin)*2
        b_ecc = (b_rmax - b_axis)/b_axis
        print("Star A (Yellow):\nSemi-major axis = {0:8.3f} AU, Period = {1:8.3f} years, Eccentricity = {2:8.3f}".format((a_axis/AU), (a_period/YEAR), a_ecc))
        print("Star B (Red):\nSemi-major axis = {0:8.3f} AU, Period = {1:8.3f} years, Eccentricity = {2:8.3f}".format((b_axis/AU), (b_period/YEAR), b_ecc))
        
        counter = 0
    
    counter += 1
    
    acheck = 1
    bcheck = 1
    
    t += h
    rate(100)


#end script







