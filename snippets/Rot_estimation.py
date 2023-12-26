from geo_algebra import *
m_vectors = 5

q_lst = []
p_lst = []

R = rdn_rotor()
noise = 0.1

for i in range(m_vectors):
    q = rdn_vanilla_vec()
    p = ~R*q*R + noise*rdn_vanilla_vec() # Don't forget to allways add a bunch of noise
    q_lst = [q] + q_lst
    p_lst = [p] + p_lst


R_est = estimate_rot(p_lst,q_lst)

p_est_lst = trans_list(q_lst,R_est)

print(np.argmin(comp_dist_matrix(q_lst,p_lst),axis=1))
print(comp_dist_matrix(q_lst,p_lst))
print("P estimated:")
print(np.argmin(comp_dist_matrix(p_est_lst,p_lst),axis=1))
print(comp_dist_matrix(p_est_lst,p_lst))

print('Correspondences from magnitudes:')
print(get_corr(q_lst,p_lst))

for i in range(5):
    p = p_lst[i]
    q = q_lst[i]
    print(mag(R_est*p*~R_est - q),q)
