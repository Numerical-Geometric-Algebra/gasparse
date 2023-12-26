from geo_algebra import *
from random import shuffle
import matplotlib.pyplot as plt
plt.style.use('dark_background')



#freqs = [0.0001,0.0005,0.001,0.005,0.01,0.05,0.1]
freqs = list(np.arange(10,21,0.1))
#nfreqs = 30
#freqs = list(np.arange(0.001,1,0.001))
#freqs = [100,200,300]
#freqs = np.logspace(-5,5, num=nfreqs)
#freqs = [0.0001,0.001]
#freqs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
#freqs = [1,2,3,4,5,6,7,8,9,10]

nfreqs = len(freqs)
n_points = 100
mu = 0
sigma = 0

s = 5
gamma = 1
f = 0

#np.random.seed(657483)
np.random.seed(5847963)
x_lst = generate_rdn_PC(n_points)
#R = rdn_rotor()
# the rotation angle
theta = 100*np.pi/180
u = normalize(rdn_vanilla_vec())
R = np.cos(theta/2) + I*u*np.sin(theta/2)

t = 0
y_lst = apply_vec_RBM(x_lst,R,t,mu,sigma)

for i in range(len(y_lst)):
    y_lst[i] = normalize(y_lst[i])
    x_lst[i] = normalize(x_lst[i])

y_shuffled = y_lst.copy()
shuffle(y_shuffled)
#corr_array = exp_corrs(x_lst,y_lst,s,gamma,f)
#print(corr_array)

X_lst = x_lst
Y_lst = y_shuffled

#X_lst = convoution(x_lst)
#Y_lst = convoution(y_shuffled)

P_lst = simple_transform(freqs,X_lst)
Q_lst = simple_transform(freqs,Y_lst)

R_est_lst = [0]*len(P_lst)
#angle_lst = np.zeros([len(P_lst)])
angle_lst = []
freqs_lst = []
mag_lst = []
for i in range(len(P_lst)):
    R_est_lst[i] = estimate_rot([P_lst[i]],[Q_lst[i]])
    angle = 2*180*np.arccos(get_float(R*~R_est_lst[i]))/np.pi
    if angle > 180:
        angle -= 180
        #if angle < 90:
    freqs_lst += [freqs[i]]
    angle_lst += [angle]
    #mag_lst += [abs(mag(P_lst[i]) - mag(Q_lst[i]))]
    #print(R_est_lst[i])
    #print(freqs[i],' ',2*180*np.arccos(get_float(R*~R_est_lst[i]))/np.pi)
    #print(' ', normalize(I|R_est_lst[i]))

#print(R_est_lst)
#freq_array = np.array(freqs_lst)
#angle_array = np.arra
plt.plot(freqs_lst,angle_lst)
#plt.plot(freqs_lst,mag_lst)
plt.show()

#fig, ax = plt.subplots()




R_est1 = estimate_rot(P_lst,Q_lst)
P_lst1 = apply_rotation(P_lst,R_est1)
R_est2 = estimate_rot(P_lst1,Q_lst)
R_est = R_est2*R_est1

R_SVD = estimate_rot_SVD(P_lst,Q_lst)

R_est_trad = estimate_rot(x_lst,y_lst)
error = compute_PC_error(0,R_est,x_lst,y_lst)

ground_error = compute_error_euclidean(R,P_lst,Q_lst)
mvs_error = compute_error_euclidean(R_est,P_lst,Q_lst)
ground_PC_error = compute_PC_error(0,R,x_lst,y_lst)



#error = compute_error(T_est,R_est,p_lst,q_lst)

print('Sigma:',sigma)
print('N freqs:',nfreqs)
print('N points:',n_points)
print('PCs Estimation Error:',error)
print('PCs Ground Truth Error:',ground_PC_error)
print()
print('MVs Estimation Error:',mvs_error)
print('MVs Ground Truth Error:',ground_error)
print()
print('Ground Truth Rotor:',R)
print('Estimated Rotor:',R_est)
print('Rotor Difference:',mag(R-R_est))
print('Angle Difference:',2*180*np.arccos(get_float(R*~R_est))/np.pi)
