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

e12 = blades['e12']
e13 = blades['e13']
e23 = blades['e23']


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

def rdn_multivector():
    return ga.multivector(list(np.random.rand(32)),blades=list(blades.keys()))

def rdn_rotor():
    a = normalize(rdn_vanilla_vec())
    b = normalize(rdn_vanilla_vec())
    return a*b

def rdn_translator():
    t = 10*rdn_vanilla_vec()
    return 1 + (1/2)*einf*t

def P_I(X):
    return X(0) + ((X(1) + X(2) + X(3))|I)*~I

def get_coeffs(X):
    A = -P_I(einf|X)
    B = -P_I(eo|X)
    C = P_I(X|(eo^einf))
    D = P_I(X)
    
    #print(eo*A + einf*B + (eo^einf)*C + D -X)
    
    return [A,B,C,D]


def sqrt_rotor(e,u):
    e = normalize(e)
    u = normalize(u)
    scalar = 1/np.sqrt(get_float(e|u) + 1)
    return (1/np.sqrt(2))*scalar*(e*u + 1)


def mag(X):
    return ((X*~X).list()[0])

def dist(X,Y):
    A,B,C,D = get_coeffs(X-Y)
    return np.sqrt(mag(A) + mag(C) + mag(D))


def pos_cga_dist(X,Y):
    A,B,C,D = get_coeffs(X-Y)
    return np.sqrt(mag(A) +  mag(B) + mag(C) + mag(D))


def comp_dist_matrix(X_list,Y_list):
    matrix = np.zeros([len(X_list),len(Y_list)])
    for i in range(len(X_list)):
        for j in range(len(Y_list)):
            #matrix[i][j] = abs(get_float(X_list[i]|Y_list[j]))
            #matrix[i][j] = dist(X_list[i],Y_list[j])
            matrix[i][j] = pos_cga_dist(X_list[i],Y_list[j])

    return matrix

def trans_list(X_list,U):
    new_list = []
    for i in range(len(X_list)):
        new_list += [U*X_list[i]*~U] 
    return new_list

def proj(B,X):
    return (X|B)*inv(B)

def proj_perp(B,X):
    return X - proj(B,X)


def rotor_sqrt(R):
    theta = np.arccos(get_float(R(0)))
    B = normalize(R(2))
    return np.cos(theta/2) + B*np.sin(theta/2)


# Estimate  correspondences
def compute_magnitudes(X_list):
    mag_vector = np.zeros([len(X_list),1])
    for i in range(len(X_list)):
        mag_vector[i] = mag(X_list[i])
    return mag_vector

def estimate_corrs_index(X_maglist,Y_maglist):
    Matrix = np.abs(X_maglist - Y_maglist.T)
    return np.argmin(Matrix,axis=1)

def get_corr(X_list,Y_list):
    X_magarray = compute_magnitudes(X_list)
    Y_magarray = compute_magnitudes(Y_list)
    corr_index = estimate_corrs_index(X_magarray,Y_magarray)
    return corr_index

def estimate_rot(p_lst,q_lst):
    beta_matrix = np.zeros([4,4])
    basis_rotor = [blades['e'],e12,e13,e23]

    for k in range(len(p_lst)):
        p = p_lst[k]
        q = q_lst[k]

        def Func(Rotor):
            return q*Rotor*p

        for i in range(4):
            for j in range(4):
                beta_matrix[i][j] += get_float(Func(basis_rotor[i])*(~basis_rotor[j]))


    eigenvalues, eigenvectors = np.linalg.eig(beta_matrix)
    u = eigenvectors[:,np.argmax(eigenvalues)]# eigenvalues are column vectors
    R_est = 0
    for i in range(4):
        R_est += u[i]*basis_rotor[i]

    R_est = normalize(R_est)
    
    return R_est
    

def estimate_rbm(P_lst,Q_lst):
    # The optimal rotation and translation
    beta_matrix = np.zeros([4,4])

    basis_rotor = [blades['e'],e12,e13,e23]

    for k in range(len(P_lst)):
        A,B,C,D = get_coeffs(P_lst[k])
        E,F,G,H = get_coeffs(Q_lst[k])

        def F(Rotor):
            return E(1)*Rotor*A(1) - E(2)*Rotor*A(2)
            
        for i in range(4):
            for j in range(4):
                beta_matrix[i][j] += get_float(F(basis_rotor[i])*(~basis_rotor[j]))
        
    eigenvalues, eigenvectors = np.linalg.eig(beta_matrix)
    R_est = 0

    u = eigenvectors[:,np.argmax(eigenvalues)]
    for i in range(4):
        R_est += u[i]*basis_rotor[i]

    R_est = normalize(R_est)
    def Rot(X):
        return R_est*X*~R_est

    scalar = 0
    v = 0
    for i in range(len(P_lst)):
        A,B,C,D = get_coeffs(P_lst[k])
        E,F,G,H = get_coeffs(Q_lst[k])
        
        v += ((Rot(A) + E)*~(G+H-Rot(C+D)))(1)
        scalar += mag(A) + mag(E)

    t_est = v*(1/scalar)
    T_est = 1 + (1/2)*einf*t_est

    return (T_est,R_est)