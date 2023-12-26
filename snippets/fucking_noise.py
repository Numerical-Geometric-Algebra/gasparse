import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')



def transform(s_array,x_array):
    F = np.zeros([len(s_array)])

    for i in range(len(s_array)):
        exp_value = np.cos(-s_array[i]*abs(x_array))
        F[i] = np.sum(x_array*exp_value)
    
    return F

n_points = 200
mu = 0
sigma = 10

s_array = np.arange(0.01,10.01,0.01)
#s_array = np.arange(1,2,1)
x_array = np.random.rand(n_points)

F = transform(s_array,x_array)

n_noises = 1000

F_noise_array = np.zeros([n_noises,len(s_array)])

for i in range(n_noises):
    noise = np.random.normal(mu, sigma, n_points)
    F_noise = transform(s_array,x_array + noise)
    F_noise_array[i] = F_noise
    #plt.plot(s_array,abs((F_noise-F)/F))
    #plt.plot(s_array,F_noise)


#plt.plot(s_array,F_noise_array.T)
#plt.plot(F_noise_array)
plt.plot(s_array,F,'o-')
plt.plot(s_array,np.mean(F_noise_array,axis=0),'ro-')

#plt.plot(noise)
plt.show()