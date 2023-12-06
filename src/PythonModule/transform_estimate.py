from geo_algebra import *
freqs = [1,2,3,4,5,6]
n_points = 10
mu = 0
sigma = 0

x_lst = generate_rdn_PC(n_points)
R = rdn_rotor()
t = 10*rdn_vanilla_vec()
y_lst = apply_vec_RBM(x_lst,R,t,mu,sigma)
p_lst = vanilla_to_cga_vecs(x_lst)
q_lst = vanilla_to_cga_vecs(y_lst)

#R,t,p_lst,q_lst = generate_rdn_PCs(100)



P_lst = transform(freqs,p_lst)
Q_lst = transform(freqs,q_lst)

T_est,R_est = estimate_rbm(P_lst,Q_lst)

# P_est = ~R_est*~T_est*Q_lst[0]*T_est*R_est
t_est = -2*eo|T_est


error = compute_PC_error(t_est,R,x_lst,y_lst)
#error = compute_error(T_est,R_est,p_lst,q_lst)

print(error)
