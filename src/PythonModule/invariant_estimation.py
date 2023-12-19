import numpy as np
import matplotlib.pyplot as plt
from geo_algebra import*
import matplotlib.cm as cm
from scipy.stats import ncx2, norm

plt.style.use('dark_background')


def rdn_rot_matrix(theta):
    basis_vecs = [e1,e2,e3]
    theta = theta*np.pi/180
    u = normalize(rdn_vanilla_vec())
    R = np.cos(theta/2) + I*u*np.sin(theta/2)

    R_matrix = np.zeros([3,3])

    for i in range(3):
        for j in range(3):
            R_matrix[i][j] = get_float(basis_vecs[i]|(R*basis_vecs[j]*~R))
    
    return R_matrix


def estimate_rot(x,y,w):
    matrix = np.einsum('kj,lj,mk->lm',w,y,x)
    U, S, V = np.linalg.svd(matrix, full_matrices=True)
    M = np.array([[1,0,0],
                  [0,1,0],
                  [0,0,np.linalg.det(U)*np.linalg.det(V)]])
    R_est = U@M@V
    return R_est

def estimate_odd(x_array,y_array,sigma_xi,df,epsilon):
    # returns prob[i][j] = P(||x_j|^2 - eps| < |Y_i|^2 < |x_j|^2 + eps)

    x_normsq = np.sum(x_array*x_array,axis=0)
    y_normsq = np.sum(y_array*y_array,axis=0)

    Prob_matrix = np.zeros([len(y_normsq),len(x_normsq)])
    for i in range(len(y_normsq)):
        Prob_matrix[i] = ncx2.cdf((x_normsq + epsilon)/(sigma_xi**2), df, y_normsq[i]/(sigma_xi**2))
        mask = epsilon < x_normsq
        Prob_matrix[i][mask] -= ncx2.cdf((x_normsq[mask] - epsilon)/(sigma_xi**2), df, y_normsq[i]/(sigma_xi**2))
    
    return Prob_matrix

def get_weights(x_array,y_array,sigma_n,sigma_xi,df):
    # Determine the magnitude square of the vectors
    x_normsq = np.sum(x_array*x_array,axis=0)
    y_normsq = np.sum(y_array*y_array,axis=0)

    Prob_matrix = np.zeros([len(y_normsq),len(x_normsq)])
    for i in range(len(y_normsq)):
        Prob_matrix[i] = ncx2.pdf(x_normsq/(sigma_xi**2), df, y_normsq[i]/(sigma_xi**2))/(sigma_xi**2)

    L_vector = np.prod(Prob_matrix,axis=0)
    odd = Prob_matrix - L_vector
    #Q_odd = Prob_matrix/(1-Prob_matrix)
    #odd = np.einsum('j,ij->ij',L_odd,Q_odd)
    

    norm_diff_array = (y_normsq - np.expand_dims(x_normsq,axis=0).T).T/2

    Prob_matrix = np.zeros([len(y_normsq),len(x_normsq)])
    for j in range(len(x_normsq)):
        Prob_matrix[j] = norm.pdf(norm_diff_array[j],loc=0,scale=sigma_n*np.sqrt(x_normsq[j]))

    return (Prob_matrix,odd)


def get_angle(rot_matrix):
    eigenvalues, eigenvectors = np.linalg.eig(rot_matrix)

    # The axis of rotation
    u = np.real(eigenvectors[:,np.argmax(eigenvalues)])
    Kn = np.array([[0,-u[2],u[1]],
                   [u[2],0,-u[0]],
                   [-u[1],u[0],0]])

    cos_theta = (np.trace(rot_matrix) - 1)/2
    sin_theta = -np.trace(Kn@rot_matrix)/2

    return 180*np.angle(cos_theta + 1j*sin_theta)/np.pi


def plot_statistics(mu=0):

    # Sanity check to understand the non central chi-squared distribution 
    # when X is non standardized
    # TODO:
    # - Adjust the linspace better...
    sigma_lst = [1,0.5,0.4,0.3,0.2,0.1,0.01]
    mu_vec = np.array([[0.1,0.2,0.3]]).T

    for sigma in sigma_lst:
        print(sigma)
        x = mu_vec + np.random.normal(mu, sigma, [3,10000])
        x_normsq = np.sum(x*x,axis=0)

        nc = np.sum(mu_vec*mu_vec)/(sigma**2)
        x_sq = np.linspace(0,20,10000)
        #x = np.linspace(ncx2.ppf(0.01, df, nc),
        #        ncx2.ppf(0.99, df, nc), 100)
        plt.plot(x_sq, ncx2.pdf(x_sq/(sigma**2), df, nc)/(sigma**2))

        q25, q75 = np.percentile(x_normsq, [25, 75])
        bin_width = 2 * (q75 - q25) * len(x_normsq) ** (-1/3)
        bins = round((x_normsq.max() - x_normsq.min()) / bin_width)
        plt.hist(x_normsq,bins = bins, density = True)
        plt.title('Sigma=' + str(sigma))

        plt.show()

def noisy_iterations(x_array,y_array,mu,sigma_n,sigma_xi,df,n_iter=10,sigma_iter=0.001):
    for i in range(n_iter):
        Prob_matrix, odd = get_weights(x_array,y_array,sigma_n,sigma_xi,df)
        R_est_p = estimate_rot(x_array,y_array,Prob_matrix.T)

        noise = np.random.normal(mu,sigma_iter,x_array.shape)
        x_array = R_est_p@x_array + noise
    return R_est_p



# Given two point clouds it determines the rotation in a single iteration
def simple_algorithm(x_array,y_array,sigma_xi,df,epsilon=0.01):
    P = estimate_odd(x_array,y_array,sigma_xi,df,epsilon)
    # Estimate rotation
    return estimate_rot(x_array,y_array,P.T)


def dist_estimation_alg(x_array,y_array,sigma_xi,df,epsilon=0.01,n_iter=1):
    rho = 0.001
    P = estimate_odd(x_array,y_array,sigma_xi,df,epsilon)
    R0 = estimate_rot(x_array,y_array,P.T)

    index = np.argmax(P,axis=1)
    a_array = x_array.T[index].T
    z_array = R0@a_array

    zy_norm = np.sum((z_array - y_array)*(z_array - y_array),axis=0)
    mask = zy_norm < rho
    a_array = x_array.T[index[mask]].T

    P = estimate_odd(a_array,y_array,sigma_xi,df,epsilon)
    R1 = estimate_rot(a_array,y_array,P.T)

    return R0,R1


def mean_estimation_alg(x_array,y_array,sigma_xi,df,epsilon=0.01,n_iter=1):
    rho = 0.5
    rho_iter = 1.1

    R_est = np.eye(3)
    
    mean_mag = np.sum(y_array*y_array)/(y_array.shape[0]*y_array.shape[1])
    mean_mag += np.sum(x_array*x_array)/(x_array.shape[0]*x_array.shape[1])
    mean_mag /= 2
    #print(mean_mag)

    mu_array = np.copy(y_array)
    for i in range(n_iter):
        P = estimate_odd(x_array,mu_array,sigma_xi,df,epsilon)
        
        # Estimate rotation
        R_est_i  = estimate_rot(x_array,mu_array,P.T)

        # Estimate correspondences
        index = np.argmax(P,axis=1)
        value = np.max(P,axis=1)
        
        a_array = x_array.T[index].T
        z_array = R_est_i@a_array
        R_est = R_est_i

        # Update the correspondences
        zy_norm = np.sum((z_array - y_array)*(z_array - y_array),axis=0)*value
        zy_mask = zy_norm/(np.sum(zy_norm)*mean_mag) < rho
        #print(zy_norm/(np.sum(zy_norm)*mean_mag))
        rho /= rho_iter
        mu_array.T[zy_mask] = z_array.T[zy_mask]
        mu_array.T[~zy_mask] = np.copy(y_array.T[~zy_mask])
        
    
    return R_est


#np.random.seed(8743)
#np.random.seed(78347)
#np.random.seed(7744)

#np.random.seed(74582)

sigma_list = [0.0001,0.001,0.002,0.003,0.005,0.01,0.05,0.1]

n_points = 100
n_outliers = 0
mu = 0
sigma = 0.001
df = 3
epsilon = 0.1 # Probability Interval

sigma_n = 0.1 # variance between the source and target PC
sigma_xi = 0.1 # Uncertainty between magnitudes of the source PC
R_matrix = rdn_rot_matrix(100)

x_array = np.random.rand(3,n_points) - 0.5

#x_array[:,0] = 10*np.random.rand(3) - 5
#x_array[:,1] = 20*np.random.rand(3) - 10
#x_array[:,1] = 30*np.random.rand(3) - 15

y_array = R_matrix@x_array + np.random.normal(mu, sigma, [3,n_points])
y_array = np.concatenate((y_array,np.random.rand(3,n_outliers)),axis=1)


# Estimate with known correspondences (Sanity Check)
weights = np.zeros([n_points+n_outliers,n_points])
weights[:n_points] = np.eye(n_points)
R_est_sanity = estimate_rot(x_array,y_array,weights.T)


x_shuffled = np.copy(x_array).T
np.random.shuffle(x_shuffled)
x_shuffled = x_shuffled.T


'''

for i in range(10):
    R_est = mean_estimation_alg(x_shuffled,y_array,sigma_xi,df,epsilon=epsilon,n_iter=i)
    print(get_angle(R_est.T@R_matrix))
'''

'''
ux_array = np.copy(x_array)
sigma_iter = 0.03
n_iter = 20
#R_est = noisy_iterations(x_array,y_array,mu,sigma_n,sigma_xi,df)
R_est = np.eye(3)

for i in range(n_iter):
    Prob_matrix, odd = get_weights(x_array,y_array,sigma_n,sigma_xi,df)
    R_est_i = estimate_rot(x_array,y_array,Prob_matrix.T)
    R_est = R_est_i@R_est
    noise = np.random.normal(mu,sigma_iter,x_array.shape)
    x_array = R_est_i@x_array + noise
    print(get_angle(R_est.T@R_matrix))
'''

'''
When adding noise the uncertainty of the correspondences will become greater
'''

ux_array = np.copy(x_array)
for sigma in sigma_list:
    y_array = R_matrix@x_array + np.random.normal(mu, sigma, [3,n_points])
    y_array = np.concatenate((y_array,np.random.rand(3,n_outliers)),axis=1)

    Prob_matrix, odd = get_weights(x_array,y_array,sigma_n,sigma_xi,df)
    odd = estimate_odd(x_array,y_array,sigma_xi,df,epsilon)

    P_total = np.sum(odd,axis=0)/len(odd) # Total chance of matching with [x_1,x_2,...,x_m]
    # The most distinct magnitudes should have smaller chances

    w = np.einsum('ij,j->ij',odd,1/(P_total + 1E-15))


    #comp_array = np.concatenate(([np.max(odd.T,axis=0)],[np.argmax(odd.T,axis=0)],[np.arange(0,n_points)]),axis=0).T


    R_est = estimate_rot(x_array,y_array,(Prob_matrix*odd).T)
    R_est_p = estimate_rot(x_array,y_array,Prob_matrix.T)
    R_est_o = estimate_rot(x_array,y_array,odd.T)
    R_est_w = estimate_rot(x_array,y_array,w.T)
    R_est_z = estimate_rot(x_array,y_array,(Prob_matrix*w).T)


    corr_o = np.sum((np.arange(0,n_points,1) - np.argmax(odd,axis = 0)) == 0)
    corr_p = np.sum((np.arange(0,n_points,1) - np.argmax(Prob_matrix,axis = 0)) == 0)
    corr_w = np.sum((np.arange(0,n_points,1) - np.argmax(w,axis = 0)) == 0)
    print()
    print("Sigma:",sigma)
    print('Best Correspondences:')
    print('From Odd:',corr_o)
    print('From Prob:',corr_p)
    print('From Weight:', corr_w)


    print('With Odd:', get_angle(R_est.T@R_matrix))
    print('Without Odd:', get_angle(R_est_p.T@R_matrix))
    print('Only With Odd:', get_angle(R_est_o.T@R_matrix))
    print('Only With Weight:', get_angle(R_est_w.T@R_matrix))
    print('With Weight and Prob:', get_angle(R_est_z.T@R_matrix))

    






# TODO:
# - [x] Calculate angle difference between the estimated and ground truth rotations
# - Print best correspondences estimates when outlier ratio is greater then inlier
# - Quantitative analysis of noise in a point-cloud of a real data-set