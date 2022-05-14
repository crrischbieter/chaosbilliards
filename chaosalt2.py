import numpy as np
import matplotlib.pyplot as plt

#Initialize variables
N = 250
M = 10

#Define barrier as four corners of the box
b1 = (1,1)
b2 = (1,-1)
b3 = (-1,-1)
b4 = (-1,1)

#Generate file names
filename=[]
for I in range(M):
    filestr = 'chaosx'+str(I+1)+'.dat'
    filename.append(filestr)

#We need the Euclidean norm
def eunorm(v):
    return np.sqrt(v[0]**2+v[1]**2)
#and the infty norm
def inftynorm(v):
    return max(np.abs(min(v)),np.abs(max(v)))
#Function that finds the intersection between any two line segments
def intersection1(v1,v2,v3,v4):
    b21 = (v2[1]-v1[1])/(v2[0]-v1[0])
    b43 = (v4[1]-v3[1])/(v4[0]-v3[0])
    xout = (v3[1]-v1[1]+b21*v1[0]-b43*v3[0])/(b21-b43)
    yout = b21*(xout-v1[0])+v1[1]
    return (xout,yout)
#Intersection1 cannot deal with vertical lines, so here is a function specifically for that
def intersection2(v1,v2,v3,v4):
    b21 = (v2[1]-v1[1])/(v2[0]-v1[0])
    xout = v3[0]
    yout = b21*(xout-v1[0])+v1[1]
    return (xout,yout)

#Loop through different initial values
for K in range(M):
    #Reinitialize variables
    xi = 0
    yi = (K/M)
    alpha0 = 60 * ( np.pi / 180 )
    xf = 10 * np.cos(alpha0) + xi
    yf = 10 * np.sin(alpha0) + yi
    vi = (xi,yi)
    vf = (xf,yf)
    x=[]
    y=[]
    x.append(xi)
    y.append(yi)
    #For the initial ray, we will need to determine the "direction" of the ball.
    #By convention, I choose towards the positive quadrant.
    #Determine two intersections
    z1 = intersection2(vi,vf,b1,b2)
    z1 = (round(z1[0],1),z1[1])
    z2 = intersection1(vi,vf,b4,b1)
    z2 = (z2[0],round(z2[1],1))
    #Determine which of the two intersections happens first.
    if eunorm(z1) < eunorm(z2):
        vi = z1
        xf = -10 * np.cos(alpha0) + vi[0]
        yf = 10 * np.sin(alpha0) + vi[1]
    else:
        vi = z2
        xf = 10 * np.cos(alpha0) + vi[0]
        yf = -10 * np.sin(alpha0) + vi[1]

    x.append(vi[0])
    y.append(vi[1])
    vf = (xf,yf)

    for I in range(N):
        #Determine intersection between ray and barriers.
        #Rounding to keep numbers nice looking.
        z1 = intersection2(vi,vf,b1,b2)
        z1 = (round(z1[0],1),z1[1])
        z2 = intersection1(vi,vf,b2,b3)
        z2 = (z2[0],round(z2[1],1))
        z3 = intersection2(vi,vf,b3,b4)
        z3 = (round(z3[0],1),z3[1])
        z4 = intersection1(vi,vf,b4,b1)
        z4 = (z4[0],round(z4[1],1))
        #One of these intersections should be the initial starting point so it needs to be excluded.
        #If any of the intersections are outside of the box, it should also be excluded.
        for J in [z1,z2,z3,z4]:
            if J == vi:
                del J
            elif inftynorm(J) > 1:
                del J
            else: 
                z = J
        #The resulting z will be the intersection point, which we take to be the next xi and yf.
        #This step will be more complex in the next version of the sim.
        vii = vi
        vi = z
        x.append(vi[0])
        y.append(vi[1])
        if vi[0] == 1: #Bounce off right verticle wall.
            if vi[1] > vii[1]: #From below
                xf = -10 * np.cos(alpha0) + vi[0]
                yf = 10 * np.sin(alpha0) + vi[1]
            if vi[1] < vii[1]: #From above
                xf = 10 * np.cos(alpha0) + vi[0]
                yf = 10 * np.sin(alpha0) + vi[1]
        elif vi[0] == -1: #Bounce off left verticle wall.
            if vi[1] > vii[1]: #From below
                xf = 10 * np.cos(alpha0) + vi[0]
                yf = 10 * np.sin(alpha0) + vi[1]
            if vi[1] < vii[1]: #From above
                xf = -10 * np.cos(alpha0) + vi[0]
                yf = 10 * np.sin(alpha0) + vi[1]
        elif vi[1] == 1: #Bounce off top wall
            if vi[0] > vii[0]: #From left
                xf = 10 * np.cos(alpha0) + vi[0]
                yf = -10 * np.sin(alpha0) + vi[1]
            elif vi[0] < vii[0]: #From right
                xf = (10 * np.cos(alpha0) + vi[0])
                yf = 10 * np.sin(alpha0) + vi[1]
        elif vi[1] == -1: #Bounce off bottom wall
            if vi[0] > vii[0]: #From left
                xf = 10 * np.cos(alpha0) + vi[0]
                yf = 10 * np.sin(alpha0) + vi[1]
            elif vi[0] < vii[0]: #From right
                xf = (10 * np.cos(alpha0) + vi[0])
                yf = (-10 * np.sin(alpha0) + vi[1])
        vf = (xf,yf)
    
    #Write intersection points to .dat files
    f = open(filename[K],'w')
    for I in range(len(x)):
        xstr = str(x[I])
        ystr = str(y[I])
        f.write(xstr)
        f.write(' ')
        f.write(ystr)
        f.write('\n')
    f.close
    plt.plot(x,y)

##Commeting all of this out because I have no idea what any of this is doing.
##I don't know if this is even a valid way to calculate the Lyapunov exponent.
#Calculate Lyapunov exponents from the saved .dat files (must be done outside the loops)
#h = open('lyapunov2.dat','w')
#h.write('Variation of y:\n')
#for I in range(len(filename)-1):
#    f = open(filename[I])
#    g = open(filename[I+1])
#
#    flist = f.readlines()
#    glist = g.readlines()
#
#    f.close
#    g.close
#    #Break the data into x and y
#    x1 = []
#    y1 = []
#    x2 = []
#    y2 = []
#    for K in range(len(flist)):
#        temp = flist[K].split(' ')
#        x1.append(float(temp[0]))
#        y1.append(float(temp[1]))
#    for K in range(len(glist)):
#        temp = glist[K].split(' ')
#        x2.append(float(temp[0]))
#        y2.append(float(temp[1]))
#    #Now to actually calculate the Lyapunov
#    d0 = (eunorm((x1[0]-x2[0],y1[0]-y2[0])))
#    lyap = []
#    for K in range(len(x1)-1):
#        lyap.append(np.log(np.abs((eunorm((x1[I+1]-x2[I+1],y1[I+1]-y2[I+1])))/d0)))
#    #Write calculated exponents to data file
#    strg = 'Average Lyapunov for  ' + str(I+1) + ' and ' + str(I+2) + '\n'
#    h.write(strg)
#    avelya = sum(lyap) / len(lyap)
#    h.write(str(avelya) + '\n')
#h.close



boxx = [-1,1,1,-1,-1]
boxy = [1,1,-1,-1,1]
plt.plot(boxx,boxy,color="#339900")
plt.show()