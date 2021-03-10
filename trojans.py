from vpython import *
#GlowScript 2.7 VPython
'''
E3-1: trojan.py
model the orbit of 2010 TK7 around the L4 Lagrange point of Earth's orbit
use asteroid simulation developed in class (Euler)
'''

#define constants
G = 6.67E-11
AU = 1.5E11
year = 365.25*24*3600


#define objects
sun = sphere(pos = vec(0,0,0), mass = 2E30, radius = 1E10, color = color.yellow)

earth = sphere(pos = vec(AU,0,0), mass = 6E24, radius = 5E9, color = color.blue)
earth.vel = vec(0,sqrt(G*sun.mass/AU),0) #starting velocity, based on Earth's position relative to sun
earth.trail = curve(pos = earth.pos, radius = 0.08*earth.radius, color = earth.color)

L4 = sphere(pos = vec(AU*cos(pi/3),AU*sin(pi/3),0), mass = 6E24, radius = 3E9, color = color.white) #L4 Lagrange point
L4.vel = vec(-sqrt(G*sun.mass/(AU))*sin(pi/3.0),sqrt(G*sun.mass/(AU))*cos(pi/3.0), 0) #same magnitude as earth.vel, just tweaked to be perpendicular to distance between L4 and sun

v = mag(L4.vel) #assigning a variable to the magnitude of L4's velocity for easier inclusion below
asteroid = sphere(pos = 1.0004*L4.pos, radius = 3E9, color = color.yellow) #2010 TK7, a Trojan asteroid at L4
asteroid.vel = vec(-v*(1 + 0.00004*random())*cos(pi/6.2), v*(1 + 0.00004*random())*sin(pi/6.2), 0)
    

#set up plotting
plot = gdisplay(x=0,y=400, height=400, width=400, title="2010 TK7: Orbital Path", xtitle = "x", ytitle = "y")
data = gdots(color = color.blue)


#simulation
h = 1.0E4
t = 0


#using simply Euler method to meausure the objects' velocities and positions over time
while True: #means it will run forever, until the program is stopped by the user.
    asteroid.vel += (-G*sun.mass*(asteroid.pos-sun.pos)/mag(asteroid.pos-sun.pos)**3 + -G*earth.mass*(asteroid.pos - earth.pos)/mag(asteroid.pos - earth.pos)**3)*h
    asteroid.pos += asteroid.vel*h
    
    earth.vel += (-G*sun.mass*(earth.pos-sun.pos)/mag(earth.pos-sun.pos)**3)*h
    earth.pos += earth.vel*h
    earth.trail.append(pos = earth.pos, color = earth.color)
    
    L4.vel += (-G*sun.mass*(L4.pos-sun.pos)/mag(L4.pos-sun.pos)**3)*h
    L4.pos += L4.vel*h 
    
    data.plot(asteroid.pos.x/mag(L4.pos),asteroid.pos.y/mag(L4.pos)) #generating a plot depicting 2010 TK7's orbit
    
    t+=h
    rate(1000)

#end script
