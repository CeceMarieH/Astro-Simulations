from vpython import *
#GlowScript 2.7 VPython
'''
rndwlk in 3D
'''
#define moving in a lattice

def move():
    for a in arange(no_particles):
        rand = random()
        
        if rand <= 1/6:
            part_list[a].pos.x = part_list[a].pos.x + step_length
        else if rand > 1/6 and rand <= 1/3:
            part_list[a].pos.x = part_list[a].pos.x - step_length
        else if rand > 1/3 and rand <= (1/2-1/100):
            part_list[a].pos.y = part_list[a].pos.y + step_length
        else if rand > (1/2-1/100) and rand <= 2/3:
            part_list[a].pos.y = part_list[a].pos.y - step_length
        else if rand > 2/3 and rand <= 5/6:
            part_list[a].pos.z = part_list[a].pos.z + step_length
        else if rand > 5/6:
            part_list[a].pos.z = part_list[a].pos.z - step_length
            
        if tails:
            part_list[a].trail.append(pos=part_list[a].pos)
            
scene = display(title='Random Walk 2D', x=300, y=0, width = 800, height = 800)
graph1 = gdisplay(x=0,y=200, width =300, height=200, title='Average Distance vs. time')
graph2 = gdisplay(x=0,y=0, width =300, height=200, title='Diffusion Constant vs. time')

no_particles = 50
part_radius = 0.5
step_length = 1
avg_distance = 0
distance_squared = 0
t = 0
tails = True

x_axis = cylinder(pos=vec(-60,0,0), axis = vec(120,0,0), radius=0.3)
y_axis = cylinder(pos=vec(0,-60,0), axis = vec(0,120,0), radius=0.3)
z_axis = cylinder(pos=vec(0,0,-60), axis = vec(0,0,120), radius=0.3)
x_axislabel = label(pos=x_axis.axis/2, text ='X', xoffset=5, yoffset=0, space=x_axis.radius, height=20, box=0, line=0, opacity=0)
y_axislabel = label(pos=y_axis.axis/2, text ='Y', xoffset=5, yoffset=0, space=y_axis.radius, height=20, box=0, line=0, opacity=0)
z_axislabel = label(pos=z_axis.axis/2, text ='Z', xoffset=5, yoffset=0, space=z_axis.radius, height=20, box=0, line=0, opacity=0)

funct1 = gcurve(gdisplay=graph1, color=color.yellow)
funct2 = gcurve(gdisplay = graph1, color=color.green)
funct3 = gcurve(gdisplay=graph2, color=color.blue)

part_list = []

lamp = local_light(pos=vec(0,0,0),color=color.yellow)
for i in arange(no_particles): #initialize
    hue=vec(random(), random(), random())
    part = sphere(color=hue, radius=part_radius)
    part.pos=vec(0,0,0)
    part_list.append(part)
    part_list[i].trail = curve(pos=[part_list[i].pos], color=hue, radius=0.1)

for t in arange(1000):
    funct2.plot(pos=(t,sqrt(t)))
    t = 0

while True:
    move()
    
    for a in arange(no_particles):
        avg_distance += mag(part_list[a].pos)
        distance_squared += mag(part_list[a].pos)**2
        t += 1
        avg_distance /= no_particles
        distance_squared /= no_particles
        diffusion = (distance_squared - avg_distance**2)/t
        funct1.plot(pos=(t, avg_distance))
        funct3.plot(pos=(t,diffusion))
    rate(1000)
