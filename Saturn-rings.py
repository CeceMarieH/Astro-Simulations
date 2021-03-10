from vpython import *
#GlowScript 2.7 VPython
'''
6-5: Thats_no_moon.py
simulate the motion of Mimas at 2:1 resonance with Saturn
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

#scene = canvas(background = color.white)

#initial position  CHANGE FOR DIFFERENT RESONANCES
a_A = 5.2*AU*pow(2,-2./3) #initial position of rings for 2:1 resonance

#define objects

saturn = sphere(pos = vec(0,0,0), radius = 7E10, mass = 5.69E26, color = color.yellow)
mimas = sphere(pos = vec(5.2*AU,0,0), radius = 1E10, mass = 3.8E19, color = color.blue)#5E10
mimas.vel = vec(0,sqrt(G*saturn.mass/(5.2*AU)),0)

rings_list = [] #define list to hold rings objects
for i in range(0,101):
    theta = i*2*pi/100
    rings = ellipsoid(pos = vec(a_A*cos(theta),a_A*sin(theta),0), length = 8E9, width = 5E10, height = 4E10, color = color.white)
    rings.vel = vec(-sqrt(G*saturn.mass/a_A)*sin(theta), sqrt(G*saturn.mass/a_A)*cos(theta), 0)
    rings_list.append(rings)

#set up plotting
plot = gdisplay(x=0,y=400, height=400, width=600, title="Rings of Saturn: r vs t", xtitle = "t", ytitle = "r")
data = gdots(color = color.blue)
data2 = gdots(color = color.red)
#simulation
h = 1.0E7
t = 0
while True:
    for rings in rings_list:
        rings.vel += (-G*saturn.mass*(rings.pos-saturn.pos)/mag(rings.pos-saturn.pos)**3 + -G*mimas.mass*(rings.pos - mimas.pos)/mag(rings.pos - mimas.pos)**3)*h
        rings.pos += rings.vel*h
        
        if rings == rings_list[0]:
            data.plot(pos = (t/year,mag(rings.pos/a_A) - 1))
        else if rings == rings_list[101]:
            data2.plot(pos = (t/year,mag(rings.pos/a_A) - 1))
    
    mimas.vel += (-G*saturn.mass*(mimas.pos-saturn.pos)/mag(mimas.pos-saturn.pos)**3)*h
    
    mimas.pos += mimas.vel*h    
    
    t+=h
    rate(1000)
    
    
#for elliptical orbit with vmax at perigee and vmin at apogee: e = (vmax-vmin)/(vmax+vmin)


#end script
        
        
