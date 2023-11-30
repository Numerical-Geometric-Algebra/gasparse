# Estimate  correspondences
import numpy as np

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



