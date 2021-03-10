from vpython import *
#GlowScript 2.7 VPython
"""
ioalone.py
simulates orbit of Jupiter's moon Io by itself
"""

#define constants
G = 6.67E-11 #gravitational constant
AU = 1.5E11 #one astronomical unit (io-jup distance)
YEAR = 365.25*24*3600 #one year, in seconds
DAY = 24*3600
C = 6.1E-13
D = 4300

#create objects for jup and io
jup = sphere(pos = vec(0,0,0), mass = 1.898E27, radius = 1E8, color = color.orange)
io = sphere(pos = vec(421700000,0,0), mass = 8.933E22, radius = 5E7, color = color.yellow)
#io = sphere(pos = vec(-16102634.7828827,-423475341.974609,-15286319.908124), mass = 8.933E22, radius = 5E7, color = color.yellow)


#initial conditions
t = 0
counter = 0
r = mag(io.pos - jup.pos)
rmin = mag(io.pos - jup.pos)
rmax = mag(io.pos - jup.pos)

io.vel = vec(-1.3,2.0*pi*r/152853.5047-7.3,0)
io.trail = curve(pos = io.pos, radius = 0.2*io.radius, color = io.color)
#io.vel = vec(17231.2106442922,-583.810327162546,229.56352413948)

#define time step
h = 1E2
scene.autoscale = 1

moon = io
#a = eu
#define acceleration as function
def acc(iopos): #iopos is dummy variable
    r = mag(iopos - jup.pos)
    return -G*jup.mass*(iopos - jup.pos)/r**3
    
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
plot = gdisplay(x=0,y=400, height=400, width=400, title="Eccentricity of Io's orbit, no other moons present", xtitle = "time (d)", ytitle = "eccentricity")
data = gdots(color = color.blue)
data2 = gdots(color = color.red)
    
while True: #loop through calcs to animate
    
    #draw path of io
    moon = io
    rk4(io)
    io.trail.append(pos=io.pos,color=io.color)

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
    
    #test plot:
    data.plot(t/DAY,ecc)
    

    if counter == 1000:
        print("ecc = {0:1.4f}".format(ecc))
        counter = 0
    
    
    t+=h
    counter+=1
    rate(1000)



#end script