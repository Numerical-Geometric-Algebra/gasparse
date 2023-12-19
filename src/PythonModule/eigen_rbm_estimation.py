from geo_algebra import *


m_points = 100
mu = 0
sigma = 0

x_lst = generate_rdn_PC(m_points)
theta = 100*np.pi/180
u = normalize(rdn_vanilla_vec())
R = np.cos(theta/2) + I*u*np.sin(theta/2)
t = 10*rdn_vanilla_vec()
T = 1 + (1/2)*einf*t
y_lst = apply_vec_RBM(x_lst,R,t,mu,sigma)

p_lst = vanilla_to_cga_vecs(x_lst)
q_lst = vanilla_to_cga_vecs(y_lst)

cga_basis = list(ga.blades(grades=[0,2]).values())

def get_func(X_lst):
    def F(Y):
        out = 0
        for i in range(len(X_lst)):
            for j in range(i,len(X_lst)):
                out += X_lst[i]*Y*X_lst[j]
        return out
    return F


P_lst,lambda_P = eigen_decomp(get_func(p_lst),cga_basis)
Q_lst,lambda_Q = eigen_decomp(get_func(q_lst),cga_basis)


# Check if Y_ordered are eigenmultivectors of Func with eigenvalues eigenvalues_ordered
Func = get_func(q_lst)
for i in range(len(Q_lst)):
    #print(eigenvalues[i]*eigenvectors[:,i] - beta@eigenvectors[:,i])
    print(np.max(np.abs(np.array((lambda_Q[i]*Q_lst[i] - Func(Q_lst[i])).list()))))

P_lst = normalize_null_mvs(P_lst)
Q_lst = normalize_null_mvs(Q_lst)

Q_est_lst = trans_list(P_lst,T*R)

#T_est,R_est = estimate_rbm(P_lst[0:2],Q_lst[0:2])

#T_est,R_est = estimate_rbm([-P_lst[4],-P_lst[5]],Q_lst[4:6])

# How to determine the right sign for the eigenmultivectors???
# For vectors I can use the eo vector

#Q_est_lst = trans_list(P_lst,~R_est*~T_est)

#for i in range(len(q_lst)):
#    print(q_est_lst[i]|q_lst[i])

'''
matrix = np.zeros([30,30])
for i in range(30):
    for j in range(30):
        matrix[i][j] = abs(abs(lambda_P[i]) - abs(lambda_Q[j]))

index = np.argmin(matrix,axis=1)
P_ordered = [P_lst[i] for i in index]
Q_est_lst = trans_list(P_ordered,T*R)
'''

#Q_est_lst = trans_list(P_lst,T*R)
#Q_est_lst = normalize_null_mvs(Q_est_lst)
#Q_lst = normalize_null_mvs(Q_lst)

#Q_est_lst = trans_list(P_lst,~R*~T)
#P_est_lst = trans_list(Q_lst,~R*~T)
'''
matrix0 = np.zeros([len(Q_lst),len(Q_lst)])
for i in range(len(Q_lst)):
    for j in range(len(Q_lst)):
        matrix0[i][j] = mag(normalize_null(Q_est_lst[i]) - normalize_null(Q_lst[j]))
'''
'''
error_dist = 0
error_mag = 0 
for i in range(len(P_lst)):
    error_dist += dist(Q_est_lst[i],Q_lst[i])
    error_mag += abs(mag(Q_est_lst[i]- Q_lst[i]))

print(error_dist)
print(error_mag)
'''
'''
t = 0
y_lst = apply_vec_RBM(x_lst,R,t,mu,sigma)
X_lst,lambda_X = eigen_decomp_pc(x_lst,grades=[1])
Y_lst,lambda_Y = eigen_decomp_pc(y_lst,grades=[1])

X_est_lst = trans_list(Y_lst,~R)


matrix1 = np.zeros([len(X_lst),len(X_lst)])
for i in range(len(X_lst)):
    for j in range(len(X_lst)):
        matrix1[i][j] = mag(X_est_lst[i] - X_lst[j])

'''
# Reciprocal blades sanity check
basis = list(ga.blades(grades=[0,2]).values())
#basis =  [e1,e2,e3,eo,einf]
#basis += [e1*e2,e1*e3,e1*eo,e1*einf,e2*e3,e2*eo,e2*einf,e3*eo,e3*einf,eo^einf]
rec_basis = reciprocal_blades(basis)
matrix = np.zeros([len(basis),len(basis)])

for i in range(len(basis)):
    for j in range(len(basis)):
        matrix[i][j] = get_float(rec_basis[i]*basis[j])

print(np.sum(abs(matrix - np.eye(len(basis)))))