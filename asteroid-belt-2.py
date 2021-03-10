from vpython import *
#GlowScript 2.7 VPython
'''
asteroid.py
simulate the motion of an asteroid at resonance with Jupiter
'''

#define constants
G = 6.67E-11
AU = 1.5E11
year = 365.25*24*3600

'''
1:2 Resonance
Resonance Condition: Pa/Pj = 1/n
P^2 = 4pi^2a^3/(G*Mtot)
aa = 0.25^(1/3) * aj
ra = 2^(-2/3) * rj
'''

scene = canvas(background = color.white)

#initial position  CHANGE FOR DIFFERENT RESONANCES
a_A = 5.2*AU*pow(2,-2./3) #initial position of asteroid for 2:1 resonance
a_B = 5.2*AU*pow(3./2,-2./3) #initial position of asteroid for 3:2 resonance
a_C = 5.2*AU*pow(1.8,-2./3)
a_D = 5.2*AU*pow(2.2,-2./3)
a_E = 5.2*AU*pow(1.6,-2./3)
a_F = 5.2*AU*pow(2.5,-2./3)

#define objects

sun = sphere(pos = vec(0,0,0), radius = 7E10, mass = 2E30, color = color.yellow)
jup = sphere(pos = vec(5.2*AU,0,0), radius = 1E10, mass = 1.898E27, color = color.green)#5E10
jup.vel = vec(0,sqrt(G*sun.mass/(5.2*AU)),0)

asteroid_list = [] #define list to hold asteroid objects
for i in range(0,101):
    theta = i*2*pi/100
    asteroid = ellipsoid(pos = vec(a_A*cos(theta),a_A*sin(theta),0), length = 1E9, width = 7E9, height = 5E9, color = color.black)
    asteroid.vel = vec(-sqrt(G*sun.mass/a_A)*sin(theta), sqrt(G*sun.mass/a_A)*cos(theta), 0)
    asteroid_list.append(asteroid)
    
for i in range(101,201):
    theta = i*2*pi/100
    asteroid = ellipsoid(pos = vec(a_B*cos(theta),a_B*sin(theta),0), length = 1E9, width = 7E9, height = 5E9, color = color.black)
    asteroid.vel = vec(-sqrt(G*sun.mass/a_B)*sin(theta), sqrt(G*sun.mass/a_B)*cos(theta), 0)
    asteroid_list.append(asteroid)
    
for i in range(201,301):
    theta = i*2*pi/100
    asteroid = ellipsoid(pos = vec(a_C*cos(theta),a_C*sin(theta),0), length = 1E9, width = 7E9, height = 5E9, color = color.black)
    asteroid.vel = vec(-sqrt(G*sun.mass/a_C)*sin(theta), sqrt(G*sun.mass/a_C)*cos(theta), 0)
    asteroid_list.append(asteroid)
    
for i in range(301,401):
    theta = i*2*pi/100
    asteroid = ellipsoid(pos = vec(a_D*cos(theta),a_D*sin(theta),0), length = 1E9, width = 7E9, height = 5E9, color = color.black)
    asteroid.vel = vec(-sqrt(G*sun.mass/a_D)*sin(theta), sqrt(G*sun.mass/a_D)*cos(theta), 0)
    asteroid_list.append(asteroid)
    
for i in range(401,501):
    theta = i*2*pi/100
    asteroid = ellipsoid(pos = vec(a_E*cos(theta),a_E*sin(theta),0), length = 1E9, width = 7E9, height = 5E9, color = color.black)
    asteroid.vel = vec(-sqrt(G*sun.mass/a_E)*sin(theta), sqrt(G*sun.mass/a_E)*cos(theta), 0)
    asteroid_list.append(asteroid)
    
for i in range(501,601):
    theta = i*2*pi/100
    asteroid = ellipsoid(pos = vec(a_F*cos(theta),a_F*sin(theta),0), length = 1E9, width = 7E9, height = 5E9, color = color.black)
    asteroid.vel = vec(-sqrt(G*sun.mass/a_F)*sin(theta), sqrt(G*sun.mass/a_F)*cos(theta), 0)
    asteroid_list.append(asteroid)

#set up plotting
plot = gdisplay(x=0,y=400, height=400, width=600, title="asteroid: r vs t", xtitle = "t", ytitle = "r")
data = gdots(color = color.blue)
data2 = gdots(color = color.red)
#simulation
h = 1.0E6
t = 0
while True:
    for asteroid in asteroid_list:
        asteroid.vel += (-G*sun.mass*(asteroid.pos-sun.pos)/mag(asteroid.pos-sun.pos)**3 + -G*jup.mass*(asteroid.pos - jup.pos)/mag(asteroid.pos - jup.pos)**3)*h
        asteroid.pos += asteroid.vel*h
        
        if asteroid == asteroid_list[0]:
            data.plot(pos = (t/year,mag(asteroid.pos/a_A) - 1))
        else if asteroid == asteroid_list[101]:
            data2.plot(pos = (t/year,mag(asteroid.pos/a_A) - 1))
    
    jup.vel += (-G*sun.mass*(jup.pos-sun.pos)/mag(jup.pos-sun.pos)**3)*h
    
    jup.pos += jup.vel*h    
    
    t+=h
    rate(1000)
    
    
#for elliptical orbit with vmax at perigee and vmin at apogee: e = (vmax-vmin)/(vmax+vmin)



        
        