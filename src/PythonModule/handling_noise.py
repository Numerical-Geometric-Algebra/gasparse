import numpy as np
import matplotlib.pyplot as plt
from geo_algebra import*
import matplotlib.cm as cm


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

def generate_mesh(lower_x,upper_x):
    #lower_x = np.max(x_array)
    #upper_x = np.min(x_array)

    nx, ny, nz = (10, 10, 10)
    x = np.linspace(lower_x, upper_x, nx)
    y = np.linspace(lower_x, upper_x, ny)
    z = np.linspace(lower_x, upper_x, nz)
    mesh = np.meshgrid(x, y, z)
    pos = np.vstack(list(map(np.ravel, mesh)))
    return pos


def get_weights(x_mesh,x_array,sigmax_sq):
    n_points = len(y_array.T)
    len_xmesh = len(x_mesh.T)
    weights = np.zeros([len_xmesh])

    for i in range(len(x_mesh.T)):
        norm = np.linalg.norm(x_mesh[:,i] - x_array.T,axis=1)
        weights[i] = np.sum(np.exp(-(1/(2*sigmax_sq))*norm))
    weights /= np.linalg.norm(weights)

    return weights

def estimate_rot(x_mesh,y_array,sigmay_sq,sigmax_sq):
    n_points = len(y_array.T)
    len_xmesh = len(x_mesh.T)
    weights = np.zeros([len_xmesh,n_points])
    for i in range(len(x_mesh.T)):
        norm = np.linalg.norm(x_mesh[:,i] - x_array.T,axis=1)
        weight = np.sum(np.exp(-(1/(2*sigmax_sq))*norm))
        norm = np.linalg.norm(x_mesh[:,i] - y_array.T,axis=1)
        weights[i] = np.exp(-(1/(2*sigmay_sq))*norm)
    weights /= np.linalg.norm(weights)

    matrix = np.einsum('kj,lj,mk->lm',weights,y_array,x_mesh)
    U, S, V = np.linalg.svd(matrix, full_matrices=True)
    M = np.array([[1,0,0],
                  [0,1,0],
                  [0,0,np.linalg.det(U)*np.linalg.det(V)]])
    R_est = U@M@V
    return R_est



def increase_prob(x_mesh,x_array,sigmax_sq):
    Id = np.eye(3)
    for i in range(len(x_mesh.T)):
        vector_diff = x_mesh[:,i]-x_array.T
        magnitude_sq = np.sum(vector_diff*vector_diff,axis=1)
        exp_value = np.exp(-(1/(2*sigmax_sq))*magnitude_sq)
        #matrix = np.einsum('k,kl,km->lm',exp_value,vector_diff,vector_diff)
        #alpha = -sigmax_sq*np.sum(exp_value)
        #if abs(alpha)  > 1E-33:
        f = -sigmax_sq*np.einsum('k,kl->l',exp_value,vector_diff)
        x_mesh[:,i] = 20*f + x_mesh[:,i]
        #u = -np.linalg.inv(alpha*Id + matrix)@f
        #x_mesh[:,i] = u + x_mesh[:,i]

    return x_mesh

n_points = 100
mu = 0
sigma = 1

sigmay_sq = 0.002
sigmax_sq = 0.002

#R = rand_rotation_matrix(0.0001)

R_matrix = rdn_rot_matrix(1)

x_array = np.random.rand(3,n_points) - 0.5
y_array = R_matrix@x_array + np.random.normal(mu, sigma, [3,n_points])
#y_array = x_array + np.random.normal(mu, sigma, [3,n_points])
# Sanity check
# print(np.linalg.norm(y_array,axis=0) - np.linalg.norm(x_array,axis=0))


x_mesh = generate_mesh(np.max(x_array),np.min(x_array))
#x_mesh = np.random.rand(3,5) - 0.5
rot_xarray = x_array

for i in range(100):
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(x_array[0],x_array[1],x_array[2])
    #ax.scatter(x_mesh[0],x_mesh[1],x_mesh[2])
    #plt.show()
    x_mesh = increase_prob(x_mesh,x_array,sigmax_sq)





'''
data = get_weights(x_mesh,x_array,sigmax_sq)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(x_mesh[0],x_mesh[1],x_mesh[2])
ax.scatter(x_mesh[0],x_mesh[1],x_mesh[2], s=300*data, c=data, cmap=cm.Oranges)
ax.scatter(x_array[0],x_array[1],x_array[2],s=500)


#ax.scatter(y_array[0],y_array[1],y_array[2])
#ax.set_zlim(lower_x, upper_x)
plt.show()
'''



data = get_weights(x_mesh,x_array,sigmax_sq)


for i in range(30):
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(rot_xarray[0],rot_xarray[1],rot_xarray[2])
    ax.scatter(x_mesh[0],x_mesh[1],x_mesh[2],s=300*data, c=data, cmap=cm.Oranges)
    ax.scatter(y_array[0],y_array[1],y_array[2])
    #ax.set_zlim(lower_x, upper_x)
    plt.show()
    
    R_est = estimate_rot(x_mesh,y_array,sigmay_sq,sigmax_sq)
    #R_est = estimate_rot_simple(x_mesh,y_array,sigmay_sq)
    x_mesh = R_est@x_mesh
    #print(R_est)
    print(np.linalg.norm(R_matrix-R_est))
    rot_xarray = R_est@x_array

#print(R_matrix)
#print(R_est)
    


'''
nx, ny, nz = (10, 10, 10)
x = np.linspace(lower_x, upper_x, nx)
y = np.linspace(lower_x, upper_x, ny)
z = np.linspace(lower_x, upper_x, nz)
mesh = np.meshgrid(x, y, z)
pos = np.vstack(list(map(np.ravel, mesh)))

weights = np.zeros([nx*ny*nz,n_points])

for i in range(len(pos.T)):
    norm = np.linalg.norm(pos[:,i] - x_array.T,axis=1)
    weight = np.sum(np.exp(-(1/(2*sigmax_sq))*norm))
    norm = np.linalg.norm(pos[:,i] - y_array.T,axis=1)
    weights[i] = np.exp(-(1/(2*sigmay_sq))*norm)
weights /= np.linalg.norm(weights)


# B matrix with known correspondences
#Bmatrix = np.einsum('lj,mj->lm',y_array,x_array)

matrix = np.einsum('kj,lj,mk->lm',weights,y_array,pos)
U, S, V = np.linalg.svd(matrix, full_matrices=True)
M = np.array([[1,0,0],
              [0,1,0],
              [0,0,np.linalg.det(U)*np.linalg.det(V)]])
R_est = U@M@V


#matrix = np.einsum('lj,mk->lm',y_array,pos)
#matrix = np.einsum('i,j->ij',y_array.T[0],x_array.T[0])
'''
'''
m = np.outer(y_array[:,0],x_array[:,0])
#for i in range(1,n_points):
#    m += np.outer(y_array[:,i],x_array[:,i])
'''
'''
xx, yy, zz = mesh
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xx, yy, zz)
ax.set_zlim(lower_x, upper_x)
plt.show()

xx, yy, zz = mesh
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_array[0],x_array[1],x_array[2])
ax.scatter(y_array[0],y_array[1],y_array[2])
ax.set_zlim(lower_x, upper_x)
plt.show()
'''


