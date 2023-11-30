# Suboptimal Rotor Estimation
from geo_algebra import *
m_multivectors = 1

Q_lst = []
P_lst = []

R = rdn_rotor()
t = 10*rdn_vanilla_vec()
T = 1 + (1/2)*einf*t

noise = 0

for i in range(m_multivectors):
    Q = rdn_multivector() # + noise*rdn_multivector() # Don't forget to allways add a bunch of noise
    P = ~R*~T*Q*T*R + noise*rdn_multivector()
    #Q_lst = [Q(1)] + [Q(2)] + [Q(3)] + [Q(4)] + Q_lst
    #P_lst = [P(1)] + [P(2)] + [P(3)] + [P(4)] + P_lst
    Q_lst = [Q(1) + Q(2) + Q(3) + Q(4)] + Q_lst
    P_lst = [P(1) + P(2) + P(3) + P(4)] + P_lst


T_est,R_est = estimate_rbm(P_lst,Q_lst)
#print(np.c_[compute_magnitudes(P_lst),compute_magnitudes(Q_lst)])
#print(get_corr(Q_lst,P_lst))

#P_est_lst = trans_list(Q_lst,~R_est*~T_est)
'''
print(np.argmin(comp_dist_matrix(Q_lst,P_lst),axis=1))
print(comp_dist_matrix(Q_lst,P_lst))
print("P estimated:")
print(np.argmin(comp_dist_matrix(P_est_lst,P_lst),axis=1))
print(comp_dist_matrix(P_est_lst,P_lst))

print('Correspondences from negative magnitudes:')
print(get_corr(Q_lst,P_lst))

'''
#print(pos_cga_dist(P_lst[0],Q_lst[0]))
m_dists = 10
if m_multivectors < m_dists:
    m_dists = m_multivectors

error = 0
for i in range(len(P_lst)):
    P = P_lst[i]
    Q = Q_lst[i]
    P_est = ~R_est*~T_est*Q*T_est*R_est
    error += dist(P_est,P)
print(error/m_multivectors)

'''
# Rotation and translation sanity check
for i in range(m_dists):
    P = P_lst[i]
    Q = Q_lst[i]
    P_est = ~R_est*~T_est*Q*T_est*R_est
    print(pos_cga_dist(P_est,P))
    print(dist(P_est,P))
'''

'''
# Rotation Sanity Check
for i in range(10):

    P = P_lst[i]
    Q = Q_lst[i]

    A,B,C,D = get_coeffs(P)
    E,F,G,H = get_coeffs(Q)

    print(mag(R_est*A*~R_est - E),E)
    #print()
'''


'''
# Sanity check!!
for i in range(10):

    P = P_lst[i]
    Q = Q_lst[i]

    A,B,C,D = get_coeffs(P)
    E,F,G,H = get_coeffs(Q)

    print(mag(R*(I*A)*~R - (I*E)),E)
    print()
'''