import numpy as np
import matplotlib.pyplot as plt

"""
#   -----1------
#   .          .
#   .          .
#   2          2
#   .          .
#   .          .
#   -----7------
#        .
#   -----4------
#   .    .     .
#   .    .     .
#   5    3     5
#   .          .
#   .    .     .
#   -----8------
#        .
#        .
#        6   

# Shaft Demensions |-------------[]-------[]
#                  |-----------L-----------|
#                                |---l---|
"""
g = 386.4 # in/sec2 -- gravity

# Linkage Demensions
L7 = 6 # in -- top shaft
l7 = 3 # in --top shaft
c7 = 1.18/2 # in -- radius of shaft 7 -- 6206-2RS SKF
A7 = np.pi*c7**2 # in2 -- Area of shaft 3
I7 = (np.pi/4)*c7**4 # in7 -- Second Moment Area 

L8 = 6 # in -- lower shaft
l8 = 3 # in -- lower shaft
c8 = 1.18/2 # in -- radius of shaft 8 -- 6206-2RS SKF
A8 = np.pi*c8**2 # in2 -- Area of shaft 3
I8 = (np.pi/4)*c8**4 # in8 -- Second Moment Area 

L3 = 5 # in -- shaft to motor 3
l3 = 3 # in -- bearing seperation on shaft 3
c3 = 1.18/2 # in -- radius of shaft 3 -- 6206-2RS SKF
A3 = np.pi*c3**2 # in2 -- Area of shaft 3
I3 = (np.pi/4)*c3**4 # in4 -- Second Moment Area 

L6 = 5 # in -- shaft to motor 6
l6 = 3 # in -- bearing seperation on shat 6
c6 = 1.18/2 # in -- radius of shat 6 -- 6206-2RS SKF Bearing Bore Diameter
A6 = np.pi*c6**2 # in -- Area of shaft 3
I6 = (np.pi/4)*c6**4 # in4 -- Second Moment Area 

L2 = 12 # in -- shaft to motor 6
c2 = 1.18 # in -- radius of shat 6 -- 6206-2RS SKF
A2 = np.pi*c2**2 # in -- Area of shaft 3
I2 = (np.pi/4)*c2**4 # in4 -- Second Moment Area 

L5 = 12 # in -- shaft to motor 6
c5 = 1.18 # in -- radius of shat 6 -- 6206-2RS SKF
A5 = np.pi*c5**2 # in -- Area of shaft 3
I5 = (np.pi/4)*c5**4 # in4 -- Second Moment Area

L4 = 6 # in -- shaft to motor 6
c4 = 1.18 # in -- radius of shat 6 -- 6206-2RS SKF
A4 = np.pi*c4**2 # in -- Area of shaft 3
I4 = (np.pi/4)*c4**4 # in4 -- Second Moment Area



# Material Properties
Ee = 29e6 # psi -- Youngs Modulus -- Steel
rhos = .285 # lb/ft3 -- density -- Steel
#Ee = 10e6 # psi -- Aluminum, YM
rhoa = .098 # lb/in3 -- density -- Aluminum


# Steppper Motor Characteristics
T = 31 # lb-in -- Holding Torque

# Weight of Appliances
def weight(rho, A, L):
    V = A*L
    return rho*V

S = 3.527 # lbs -- Motor weight 
mb = .1 # lbs -- Bearing Weight
mcl = 0 # lbs -- lower cassing Weight
mcu = 0 # lbs -- upper casing Weight

m2 = weight(rhoa, A2, L2) # lbs -- Weight of link 2
m3 = weight(rhoa, A3, L3) # lbs -- Weight of link 3
m4 = weight(rhoa, A4, L4) # lbs -- Weight of link 4
m5 = weight(rhoa, A5, L5) # lbs -- Weight of link 5
m6 = weight(rhoa, A6, L6) # lbs -- Weight of link 6
m7 = weight(rhoa, A7, L7) # lbs -- Weight of link 7
m8 = weight(rhoa, A8, L8) # lbs -- Weight of link 8


# Combined masses per section -- Extra Mass
em3 = ((m3 + 2*m2 + m7) + 2*S + 4*mb + mcu) # lbf
em6 = ((m4 + 2*m5 + m8) + 3*S + 6*mb + mcl) # lbf



# Motion Analysis

# known
r1 = L5 # in -- bottom link
r2 = L2 # in -- top link
th3 = 0 # rad
th4 = np.pi/2 # rad

# input IC
th1 = 0 # rad
s = 0 # rad -- step off of th1
th2 = th1 + s # rad

#unknown IG
r3 = 3 # in
r4 = 2 # in 

# Dynamics
mass1 = em6/g # lbf.s2/in -- equivalent mass of bottom linkage
mass2 = em3/g # lbf.s2/in -- equivalent mass of top linkage

Ii1 = (1/12)*(m5/g)*r2 # lbf.in.s2 -- Inertia of bottom linkage
Ii2 = (1/12)*(m2/g)*r2 # lbf.in.s2 -- Inertia of top linkage

P2x = 0 # lbs -- applied force x component
P2y = 0 # lbs -- applied force y component 


def E(r3,r4,th1):
    Ex = -(r1*np.cos(th1) + r2*np.cos(th1+s) - r3*np.cos(th3) - r4*np.cos(th4))
    Ey = -(r1*np.sin(th1) + r2*np.sin(th1+s) - r3*np.sin(th3) - r4*np.sin(th4))
    return np.array([Ex, Ey])

def J(r3,r4,th1):
    return np.array([(-np.cos(th3), -np.cos(th4)), (-np.sin(th3), -np.sin(th4))])

def KI(th1):
    return np.array([r1*np.sin(th1)+r2*np.sin(th1+s), -r1*np.cos(th1)-r2*np.cos(th1+s)])

def KII(th1):
    return np.array([r1*np.cos(th1)+r2*np.cos(th1+s), r1*np.sin(th1)+r2*np.sin(th1+s)])

def changeofBasis(Fx, Fy, th):
    return np.array([Fx*np.cos(th) + Fy*np.sin(th), -Fx*np.sin(th) + Fy*np.cos(th)]) # FX, FY

def linkForce_inPLane(Fx, Fy, L, l, th):
    G = np.array([(1,1,0,0),(0,0,1,1),(0,0,1,-1),(1,0,0,0)])
    cob = changeofBasis(Fx,Fy,th)
    return np.linalg.inv(G)@np.array([-cob[1], -cob[0], 0, -cob[1]*L/l]) # Ay By Ax Bx

def linkForce_outOfPlane(Fy, Fz, L, l):
    H = np.array([(1,1,0,0),(0,0,1,1),(0,0,1,0),(1,0,0,0)])
    return np.linalg.inv(H)@np.array([-Fz, -Fy, -Fy*L/l, -Fz*L/l]) # Az Bz Ay By

def maxMoment(A,B,l,L,Ee,I,st=.1):
    
    C1 = -B*l**2/6
    C3 = -B*l**2/6 + A*l**2/2
    C4 = -B*l**3/6 + A*l**3/3 - C3*l

    defl1 = np.zeros(len(np.linspace(0,l,20)))
    defl2 = np.zeros(len(np.linspace(0,l,20)))
    m1 = np.zeros(len(np.linspace(0,l,20)))
    m2 = np.zeros(len(np.linspace(0,l,20)))
    z = -1
    for x in np.linspace(0,l,20):
        z +=1
        defl1[z] = 1/(Ee*I)*(B*x**3/6 + C1*x)
        m1[z] = B*x

    z = -1
    for x in np.linspace(l,L,20):
        z +=1
        defl2[z] = 1/(Ee*I)*(B*x**3/6 + A*(x**3/6 - l*x**2/2) + C3*x + C4) 
        m2[z] = B*x + A*(x-l)
    
    x1 = np.linspace(0,l,20)
    x2 = np.linspace(l,L,20)
    x = np.hstack((x1,x2))
    cat = np.hstack((defl1,defl2))
    mcat = np.hstack((m1,m2))
    wm = np.where(mcat == 0)
    w = np.argmax(abs(cat))

    maxi = max(mcat, key = abs)
    p = np.where(mcat == maxi)
    ip = p[0]
    idx = ip[0]

    xatMaxY = x[idx]
    if xatMaxY <= l:
        M = B*xatMaxY
    else:
        M = B*xatMaxY + A*(xatMaxY - l)

    return M

def maxBendingStress(M,c,I):
    return M*c/I

def maxTorsionStress(T,d):
    return 16*T/(np.pi*d**3)

def maxAxialStress(F,A):
    return F/A

# initiating Variable
a = -1
step = .1
size  = len(np.arange(0,np.pi,step))


# ----------------------- Garbage
r33 = np.zeros(size)
r44 = np.zeros(size)

f3 = np.zeros(size)
ff3 = np.zeros(size)
f4 = np.zeros(size)
ff4 = np.zeros(size)

f1x = np.zeros(size)
f1y = np.zeros(size)
f2x = np.zeros(size)
f2y = np.zeros(size)
ff1x = np.zeros(size)
ff1y = np.zeros(size)
ff2x = np.zeros(size)
ff2y = np.zeros(size)

F01x = np.zeros(size)
F01y = np.zeros(size)
F12x = np.zeros(size)
F12y = np.zeros(size)
T1 = np.zeros(size)
T2 = np.zeros(size)

A7z = np.zeros(size)
B7z = np.zeros(size)
A7y = np.zeros(size)
B7y = np.zeros(size)

A8z = np.zeros(size)
B8z = np.zeros(size)
A8y = np.zeros(size)
B8y = np.zeros(size)

A3y = np.zeros(size)
B3y = np.zeros(size)
A3x = np.zeros(size)
B3x = np.zeros(size)

A6y = np.zeros(size)
B6y = np.zeros(size)
A6x = np.zeros(size)
B6x = np.zeros(size)

Mmax3y = np.zeros(size)
Mmax6y = np.zeros(size)
Mmax7y = np.zeros(size)
Mmax8y = np.zeros(size)
Mmax7z = np.zeros(size)
Mmax8z = np.zeros(size)

sig3 = np.zeros(size)
sig6 = np.zeros(size)
sig7 = np.zeros(size)
sig8 = np.zeros(size)

tou3 = np.zeros(size)
tou6 = np.zeros(size)
tou7 = np.zeros(size)
tou8 = np.zeros(size)


F3x = np.zeros(size)
F6x = np.zeros(size)

ns3 = np.zeros(size)
ns6 = np.zeros(size)

ag1x = np.zeros(size)
ag1y = np.zeros(size)
ag2x = np.zeros(size)
ag2y = np.zeros(size)
# ------------------------ Garbage


# CG Accelerations
al1 = 0 # rad/s2
w1 = 5 # rad/s
al2 = al1 # rads/s2


for i in np.arange(0,np.pi,step):
    a += 1
    while True:
        if abs(E(r3,r4,i).sum()) <= .001:
            r33[a] = r3
            r44[a] = r4
            
            f = np.linalg.inv(J(r33[a],r44[a],i))@KI(i)
            f3[a] = f[0]
            f4[a] = f[1]
            ff = np.linalg.inv(J(r33[a],r44[a],i))@KII(i)
            ff3[a] = ff[0]
            ff4[a] = ff[1]

            f1x[a] = -r1/2*np.sin(i) # in
            f1y[a] = r1/2*np.cos(i) # in
            ff1x[a] = -r1/2*np.cos(i) # in
            ff1y[a] = -r1/2*np.sin(i) # in

            f2x[a] = -r1*np.sin(i) - r2/2*np.sin(i+s) # in
            f2y[a] = r1*np.cos(i) + r2/2*np.cos(i+s) # in
            ff2x[a] = -r1*np.cos(i) - r2/2*np.cos(i+s) # in
            ff2y[a] = -r1*np.sin(i) - r2/2*np.sin(i+s) # in

            ag1x[a] = f1x[a]*al1 + ff2x[a]*w1**2 # in/s2
            ag1y[a] = f1y[a]*al1 + ff1y[a]*w1**2 # in/s2
            ag2x[a] = f2x[a]*al1 + ff2x[a]*w1**2 # in/s2
            ag2y[a] = f2y[a]*al1 + ff2y[a]*w1**2 # in/s2


            # P 'applied forces' must be accounted for the positive convention... must specify direction when you express them
            D = np.array([(1,0,1,0,0,0),(0,1,0,1,0,0),(0,0,-r1*np.sin(i),r1*np.cos(i),1,0),(0,0,-1,0,0,0),(0,0,0,-1,0,0),(0,0,0,0,0,1)])
            y = np.array([mass1*ag1x[a], 
                         mass1*ag1y[a] + mass1*g, 
                         Ii1*al1 + mass1/2*r1*(np.cos(i)*ag1y[a] - np.sin(i)*ag1x[a]) + r1/2*np.cos(i)*mass1*g,
                         mass2*ag2x[a] - P2x, 
                         mass2*ag2y[a] - P2y + mass2*g,
                         Ii2*al2 + r2/2*mass2*(np.cos(i+s)*ag2y[a] - np.sin(i+s)*ag2x[a]) - r2*np.cos(i+s)*P2y + r2*np.sin(i+s)*P2x + r2/2*np.cos(i+s)*mass2*g ])


            F = np.linalg.inv(D)@y
            F01x[a] = F[0]
            F01y[a] = F[1]
            F12x[a] = F[2]
            F12y[a] = F[3]
            T1[a] = F[4]
            T2[a] = F[5]

            F01cob = changeofBasis(F01x[a], F01y[a], i)
            F6x[a] = F01cob[0]

            F12cob = changeofBasis(F12x[a], F12y[a], i+s)
            F3x[a] = F12cob[0]

            # All 7 and 8 components have been broken into half due to there beng two shafts
            linkF7 = linkForce_outOfPlane(F12y[a]/2, F12x[a]/2, L7, l7)
            A7z[a] = linkF7[0]
            B7z[a] = linkF7[1]
            A7y[a] = linkF7[2]
            B7y[a] = linkF7[3]

            linkF8 = linkForce_outOfPlane(F01y[a]/2, F01x[a]/2, L8, l8)
            A8z[a] = linkF8[0]
            B8z[a] = linkF8[1]
            A8y[a] = linkF8[2]
            B8y[a] = linkF8[3]


            linkF3 = linkForce_inPLane(F12x[a], F12y[a], L6, l6, i+s) # Ay, By, Ax, Bx
            A3y[a] = linkF3[0]
            B3y[a] = linkF3[1]
            A3x[a] = linkF3[2]
            B3x[a] = linkF3[3]

            linkF6 = linkForce_inPLane(F01x[a], F01y[a], L3, l3, i)
            A6y[a] = linkF6[0]
            B6y[a] = linkF6[1]
            A6x[a] = linkF6[2]
            B6x[a] = linkF6[3]

            Mmax3y[a] = maxMoment(A3y[a], B3y[a], l3, L3, Ee, I3)
            Mmax6y[a] = maxMoment(A6y[a], B6y[a], l6, L6, Ee, I6)
            Mmax7y[a] = maxMoment(A7y[a], B7y[a], l7, L7, Ee, I7)
            Mmax7z[a] = maxMoment(A7z[a], B7z[a], l7, L7, Ee, I7)
            Mmax8y[a] = maxMoment(A8y[a], B8y[a], l8, L8, Ee, I8)
            Mmax8z[a] = maxMoment(A8z[a], B8z[a], l8, L8, Ee, I8)

            sig3[a] = maxBendingStress(Mmax3y[a],c3,I3)
            sig6[a] = maxBendingStress(Mmax6y[a],c6,I6)
            sig7y = maxBendingStress(Mmax7y[a],c7,I7)
            sig7z = maxBendingStress(Mmax7z[a],c7,I7)
            sig7[a] = abs(sig7y) + abs(sig7z)
            sig8y = maxBendingStress(Mmax8y[a],c8,I8)
            sig8z = maxBendingStress(Mmax8z[a],c8,I8)
            sig8[a] = abs(sig8y) + abs(sig8z)

            tou3[a] = maxTorsionStress(T,2*c3) # Will assume max torque since we are not accounting for it in this analysis
            tou6[a] = maxTorsionStress(T,2*c6) # same here
            tou7[a] = maxTorsionStress(T2[a],2*c7)
            tou8[a] = maxTorsionStress(T1[a],2*c8)

            ns3[a] = maxAxialStress(F3x[a],A3)
            ns6[a] = maxAxialStress(F6x[a],A6)

            break
        else:
            x = np.linalg.inv(J(r3,r4,i))@E(r3,r4,i)
            r3 = r3 + x[0]
            r4 = r4 + x[1]


# Kinematic Coeficents were checked by comparing their slops with one another considering each one is the derivative of the other 
plt.plot(np.arange(0,np.pi,.1),r33,'.k', label = 'r33')
plt.plot(np.arange(0,np.pi,.1),r44,'.r', label = 'r44')
plt.plot(np.arange(0,np.pi,.1),f3,'--', label = 'f3')
plt.plot(np.arange(0,np.pi,.1),ff3,'--', label = 'ff3')
#plt.plot(np.arange(0,np.pi,.1),f4,'.', label = 'f4')
#plt.plot(np.arange(0,np.pi,.1),ff4,'.', label = 'ff4')
plt.legend()
plt.title('Kinematics of Linkages')
plt.xlabel('Input Arm ANgle [radians]')
plt.ylabel('Kinematic Coefficents')
plt.show()

plt.plot(np.arange(0,np.pi,.1), F01x, label = 'F01x')
plt.plot(np.arange(0,np.pi,.1), F01y, label = 'F01y')
plt.plot(np.arange(0,np.pi,.1), F12x, label = 'F12x')
plt.plot(np.arange(0,np.pi,.1), F12y, label = 'F12y')
plt.plot(np.arange(0,np.pi,.1), T1, label = 'T1')
plt.plot(np.arange(0,np.pi,.1), T2, label = 'T2')
plt.legend()
plt.xlabel('Input Arm Angle [radians]')
plt.ylabel('Dynamc Loads Applied [lbf] or [lbf.in]')
plt.show()
"""
plt.plot(np.arange(0,np.pi,.1),Mmax3y)
plt.plot(np.arange(0,np.pi,.1),Mmax6y)
plt.plot(np.arange(0,np.pi,.1),Mmax7y)
plt.plot(np.arange(0,np.pi,.1),Mmax7z)
plt.plot(np.arange(0,np.pi,.1),Mmax8y)
plt.plot(np.arange(0,np.pi,.1),Mmax8z)
plt.legend(['Mmax3y','Mmax6y','Mmax7y','Mmax7z','Mmax8y','Mmax8z'])
plt.xlabel('Input Arm Angle [radians]')
plt.ylabel('Maximum Moments on Connecting Links [lbf.in]')
plt.show()
"""



    












