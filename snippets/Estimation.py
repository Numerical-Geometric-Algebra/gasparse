#!/usr/bin/env python3
import gasparse
import numpy as np
from gasparse import multivector
from sympy import poly, nroots
from sympy.abc import x

ga = gasparse.GA(4,1) # Use 3D conformal geometric algebra
blades = ga.blades()
#locals().update(blades) # Update all of the blades variables into the environment

bivectors = ['e12','e13','e14','e15','e23','e24','e25','e34','e35','e45']
vanilla_vecs = ['e1','e2','e3']

e1 = blades['e1']
e2 = blades['e2']
e3 = blades['e3']


einf = (1.0*blades['e5'] - 1.0*blades['e4'])*(1/np.sqrt(2))
eo = (1.0*blades['e5'] + 1.0*blades['e4'])*(1/np.sqrt(2))

I = blades['e1']*blades['e2']*blades['e3']
i = I*(eo^einf)


def get_float(X):
    return X.list()[0]

def inv(X):
    scalar = 1/((X*~X).list()[0])
    return X*scalar

def normalize(X):
    scalar = 1/np.sqrt((X*~X).list()[0])
    return X*scalar

# generate random bivector
def rdn_biv():
    return ga.multivector(list(np.random.rand(10)),blades = bivectors)

def rdn_vanilla_vec():
    return ga.multivector(list(np.random.rand(3)),blades = vanilla_vecs)

def rdn_rotor():
    a = normalize(rdn_vanilla_vec())
    b = normalize(rdn_vanilla_vec())
    return a*b

def rdn_translator():
    t = 10*rdn_vanilla_vec()
    return 1 + (1/2)*einf*t

def P_I(X):
    return (X|I)*~I

def get_coeffs(X):
    a = -P_I(einf|X)(1)
    b = -P_I(eo|X)(1)
    mu = (X|(eo^einf))(0)
    D = P_I(X)

    #return a*eo + b*einf + mu*(eo^einf) + D

    return [a,b,mu,D]

def polynomial(P,Q):
    a,b,mu,D = get_coeffs(P)
    e,f,gamma,H = get_coeffs(Q)

    h = I|H
    d = I|D

    c1 = get_float(e|e)
    c2 = -get_float(h|inv(e))
    c3 = get_float(a|a)
    c4 = get_float(d|d)
    c5 = get_float(a|d)
    
    return poly(c1*((x+c2)**2)*(c3 + c4*x**2 + 2*x*c5) - (x*c4 + c5)**2)

def sqrt_rotor(e,u):
    e = normalize(e)
    u = normalize(u)
    scalar = 1/np.sqrt(get_float(e|u) + 1)
    return (1/np.sqrt(2))*scalar*(e*u + 1)


def mag(X):
    return ((X*~X).list()[0])

def dist(X,Y):
    A = -P_I(einf|X)
    C = (X|(eo^einf))
    D = P_I(X)

    E = -P_I(einf|Y)
    G = (Y|(eo^einf))
    H = P_I(Y)

    return np.sqrt(mag(A-E) + mag(C-G) + mag(D-H))

def proj(B,X):
    return (X|B)*inv(B)

def proj_perp(B,X):
    return X - proj(B,X)


def rotor_sqrt(R):
    theta = np.arccos(get_float(R(0)))
    B = normalize(R(2))
    return np.cos(theta/2) + B*np.sin(theta/2)
'''
def random_point_in_circle(radius,center,normal):
    theta = np.pi*np.random.rand()
    B = normalize(I*(normal))
    u = normalize(proj(B,rdn_vanilla_vec())) + center
    u = eo + u + (1/2)*mag(u)*einf
    T = 1 + (1/2)*einf*center
    U = T*(np.cos(theta) + B*np.sin(theta))*~T
    x = U*u*~U
    #x = eo + normalize(P_I(x)) + (1/2)*einf
    return ~T*x*T
    #x = T*radius*U*u*~U*~T
    #return P_I(x) + center
'''
'''
def test_rot_freedom(a,b):
    Plane = (b-a)*I
    u = rdn_vanilla_vec()
    c = Proj(Plane,u)
    a_para = 
'''

def random_point_in_circle(center,normal):
    # computes a ramdom point in a circle that is the intersection of a sphere 
    # with a plane
    circle = (eo + center)^normal

    s = circle|einf
    s = 1/get_float(s|s)
    radius_sq = get_float(-circle*circle*s)
    d = circle*inv(-einf|circle)

    B = normalize(I*(normal))
    u = np.sqrt(radius_sq)*normalize(proj(B,rdn_vanilla_vec()))
    return u + P_I(d)

def random_common_rotation(a,b):
    # Determines a random rotation that rotates a to b
    d = random_point_in_circle((1/2)*b,a-b)
    a_perp = (a|d)*inv(d)
    a_para = (a^d)*inv(d) 

    R_sq = (b-a_perp)*inv(a_para)
    R = rotor_sqrt(R_sq)
    return R    

'''
radius = 5
a = radius*normalize(3*e2 + 0.7*e3 + 10*e1)
b = radius*normalize(7*e2 + 0.5*e3 + 9*e2)

# generate random rotor that rotates a point a to a point b
# normal = a-b


# x = random_point_in_circle((1/2)*a,normal)
d = random_point_in_circle((1/2)*b,a-b)
a_perp = (a|d)*inv(d)
a_para = (a^d)*inv(d)

R_sq = (b-a_perp)*inv(a_para)
R = rotor_sqrt(R_sq)

print(R*a*~R - b)
'''




'''
# The common plane of rotation experiment
# Generate rdn rotor
R = rdn_rotor()
a1 = rdn_vanilla_vec()
a2 = rdn_vanilla_vec()
b1 = R*a1*~R
b2 = R*a2*~R
# common plane of rotation
n = I*((b1-a1)^(b2 - a2))
R_sq = (b1 - proj(n,a1))*inv(proj_perp(n,a1))
R_est = rotor_sqrt(R_sq)

print(R_est*a1*~R_est - b1)
print(R_est*a2*~R_est - b2)
'''

# plane of rotation
#D = normalize(I*d)
#a_para = proj(D,a)
#a_perp = a - a_para
#b_para = proj(D,b)

#print(mag(b_para), mag(a_para))



#R_sq = b_para*inv(b_perp)
#print(R*b_para)


#a_para = -(1/2)*(b-a) + x
#a_perp = a - a_para

# Determine the plane of rotation
#B = normalize(a_perp*I)
#b_para = proj(B,b)
#b_perp = b - b_para

#R_sq = (b-a_perp)*inv(a_para)
#R = rotor_sqrt(R_sq)

#print(R*a*~R - b)

#u = rdn_vanilla_vec()
#u = eo + u + (1/2)*mag(u)*einf

#v = (u|circle)*circle
#s = 1/get_float(-v|einf)
#w = P_I(v)*s


#x = random_point_in_circle(radius,center,normal)




#print('radius',mag(x-center))
#print('plane',x|normal)

#print(rs)

'''
u = a

# determine the suboptimal rotor

slambda = mag(D) + np.sqrt(mag(u)*mag(e))
d = I*D
v = np.sqrt(mag(e))*normalize(u)
z = inv(d)*(slambda - a*v)

print(z)

R_est = sqrt_rotor(e,u)
t_est = (inv(e)*(gamma - mu)) + (inv(e)|(R_est*D*~R_est - H))
T_est = 1 + (1/2)*einf*t_est

# determine the optimal rotor

P_est = ~R_est*~T_est*Q*T_est*R_est

#print('alpha:',alpha)
#print(T_est - T)
#print(P-P_est)
'''



# The optimal rotation and translation
R = rdn_rotor()
t = 10*rdn_vanilla_vec()
T = 1 + (1/2)*einf*t

Q = rdn_biv()
P = ~R*~T*Q*T*R

#print(Q)
#print(P)

a,b,mu,D = get_coeffs(P)
e,f,gamma,H = get_coeffs(Q)

'''
p = polynomial(P,Q)
rs = np.array(nroots(p))

for r in rs:
    alpha = float(r)
    if abs(alpha) > 1E-10:
        u = a + alpha*(I|D)
        slambda = mag(D) + np.sqrt(mag(u)*mag(e))
        d = I*D
        v = np.sqrt(mag(e))*normalize(u)
        z = inv(d)*(slambda - a*v)
        n = I*((z-d-e)^(normalize(u) - normalize(e)))
        #B = (z-d-e)^(normalize(u) - normalize(e))
        #print(mag(z-d) - mag(alpha*e))
        R_squared = (alpha*e - proj_perp(B,z-d))*inv(proj(B,z-d))
        print(R_squared*~R_squared)
        
        R_est = rotor_sqrt(R_squared)
        #R_est = sqrt_rotor(e,u)
        t_est = (inv(e)*(gamma - mu)) + (inv(e)|(R_est*D*~R_est - H))
        T_est = 1 + (1/2)*einf*t_est
        P_est = ~R_est*~T_est*Q*T_est*R_est

        print('alpha:',alpha)
        print('dist:',dist(P,P_est))
        print('t:',t-t_est)

    #print(T_est - T)
    #print(P-P_est)
    

    #print(R_est - R)
    #print(alpha)
'''

u = a
u_hat = normalize(u)
e_hat = normalize(e)
R_est = random_common_rotation(u_hat,e_hat)
t_est = (inv(e)*(gamma - mu)) + (inv(e)|(R_est*D*~R_est - H))
T_est = 1 + (1/2)*einf*t_est
P_est = ~R_est*~T_est*Q*T_est*R_est

#print(R_est*u_hat*~R_est - e_hat)
#print(dist(P_est,P))
#print(get_coeffs(P_est))
#print(get_coeffs(P))

w = -e*I*H*inv(e) - I*H
y = -a*I*D*inv(a) - I*D
print((w - y)^(a - e))
