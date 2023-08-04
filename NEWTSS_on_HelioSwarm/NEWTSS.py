# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 11:09:27 2022

@author: teddy
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la
from itertools import combinations

# %% define functions that will be used in the script
def poly_O2(x, a, b, c):
    # a helper function for Prob_alias()
    return a*np.power(x,2) + b*x + c

def Prob_alias(k, chi, B):
    # takes the coefficients B and reports the probability of aliasing for 
    # points with (normalized) magnitude k and shape chi
    b = poly_O2(chi, B[0], B[1], B[2])
    return 0.5*(1 + np.tanh(4*np.log10(k) - b) )

def std_eqn(X, c0, c1, c2):
    # takes the values of the c coefficients and returns the values for the standard deviation
    # of log(error) for points with (normalized) magnitude k and shape chi
    k = X[0]
    chi = X[1]
    return c0 + c1*np.log10(k) + c2*np.power(chi,2)

def Error_fun(k, chi, A):
    # takes the values of the c coefficients and returns the values for the mean value
    # of error for points with (normalized) magnitude k and shape chi
    a0, a1, a2, a3, a4, a5 = A
    
    # fit equation to mean error value
    t1 = np.power(a0,chi - a1) + a2
    exp1 = np.power(a3,chi - a4) + a5
    out = t1*np.power(k,exp1)
    
    # report the value of Background_Error if the (k,chi) combination falls outside of
    # the range of scales that the equations were fitted on
    valid_scales = np.logical_and(k < 17.652, k > 0.016)
    out = np.minimum(Background_Error*valid_scales, out)
    out[out==0] = Background_Error
    return out

def find_polys(loc,n):
    # find_poly creates a list of all n length combinations of numbers of 0 through N
    N = np.size(loc,0)
    temp = list(range(N))
    return list(combinations(temp,n))

def get_positions():
    # get_positions reads the text files containing HelioSwarm node/Hub position and time
    # data and returns them as numpy arrays
    length = np.zeros(9)
    for name in range(9):
        data = np.loadtxt("input_data/HS_config/n%s_clean.txt" % name,dtype = str)
        length[name] = np.size(data[:,4])
    L = np.min(length).astype(int)
    times = data[:L,0:4]
    positions = np.zeros([9,L,3])
    for name in range(9):
        data = np.loadtxt("input_data/HS_config/n%s_clean.txt" % name,dtype = str)    
        positions[name,:,:] = data[:L,4:7].astype(float)
    return [positions, times]

def calc_RLEP(r): 
    # function calc_RLEP takes in array, each row of which represents a point in 3d space
    # it outputs a list containing the volumetric tensor R, size L, elongation E, and planarity P
    
    # find the number of satellites as the number of rows of r
    N = np.size(r[:,0])
    # find the mesocentre rb
    rb = np.mean(r, axis=0)
    # calculate the Volumetric Tensor R
    R = np.zeros([3,3])
    for i in range(N):
        R += np.outer(r[i,:]-rb, r[i,:]-rb)/N
    # find the eigenvalues of R as value in lambdas
    temp = la.eig(R)
    lambdas = temp[0]
    # find semiaxes of quasi-ellipsoid a,b,c
    # check if eigenvalues are real
    if any(np.imag(lambdas) != np.zeros(3)):
        raise ValueError('Eigenvalue has imaginary component')
    lambdas_real = np.real(lambdas)
    #print(lambdas_real)
    [c,b,a] = np.sqrt( np.sort(lambdas_real) )
    # calculate L,E,P
    L = 2*a
    E = 1 - b/a
    P = 1 - c/b
    return [R,L,E,P]

# %% Code Description
'''
This code scans through all subsets of spacecraft of a given HS configuration. For each value in the 
defined range of wavevector magnitudes (k_mag) it selects the subset of spacecraft such that the expression
    mu_error + std_coef*std_dev_error
is minimized. It then identifies the k-magnitudes for which the above expression is below the threshold 
    Er_thres. 

This code is a complement to the publication found here: doi.org/10.3847/1538-4365/acc6c7
'''
# %% define code parameters
hour = 94                       # Configuration hour of HS Phase B DRM mission (0 to 8521 for Phase B mission file)
std_coef = 2                    # see above
Er_thres = 20                   # see above

pts = 500                       # resolution in k-space to scan
Background_Error = 1*10**3      # error value assigned to points outside of valid k-domain
k_mags = np.logspace(-5,0,pts)  # wavevector magnitude values (in km) to scan
N_save = [4,5,6,7,8,9]          # n-gons to break overall configuration into
chi_thres = 1                   # minimum shape value to be included in analysis (must be less/eq to 1)

# %% loop through spacecraft configurations
b_terms = np.load('input_data/sigmoid_b_terms.npy')
print('Using NEWTSS')
print('Estimating errors for HS hour %i' %hour)
print('--------------------------------------------------')
Er_min_save = np.zeros([pts])
Std_min_save = np.zeros([pts])
L_save_N = []
chi_save_N = []
ind_save_N = np.zeros([np.size(N_save),pts]).astype(int)
Er_min_N = np.zeros([np.size(N_save),pts])
Std_min_N = np.zeros([np.size(N_save),pts])
for n in N_save:
    saved = np.load('input_data/saved_coef/coefs_n%i.npy' %n)
    a0,a1,a2,a3,a4,a5,c0,c1,c2 = np.mean(saved,axis=1)
    A = [a0,a1,a2,a3,a4,a5]
    
    [b0_vals, b1_vals, b2_vals] = np.load('input_data/saved_coef/sigma_coefs_n%i.npy' %n)
    b = np.array([b0_vals.mean(), b1_vals.mean(), b2_vals.mean()])
        
    [positions, times] = get_positions()
    r = positions[:,hour,:]
    r = r - np.mean(r,axis=0)
    
    poly_indices = find_polys(r,n)
    N_tetra = len(poly_indices)
    r_tetra = np.zeros([N_tetra,n,3])
    L_save = np.zeros([N_tetra])
    E_save = np.zeros([N_tetra])
    P_save = np.zeros([N_tetra])
    for j in range(N_tetra):
        # j is row of indices array to use points of 
        # pick out all polyhedra of order j (combinations with j points)
        index = list(poly_indices[j])
        r_tetra[j,:,:] = r[index,:]
        [R,L_temp,E,P] = calc_RLEP(r[index,:])
        L_save[j] = L_temp
        E_save[j] = E
        P_save[j] = P
    
    
    # plot E and P of subsets of N-hedra    
    plt.figure()
    plt.scatter(E_save, P_save, c=L_save)  
    plt.plot(np.linspace(0,chi_thres,100), np.sqrt(chi_thres**2 - np.linspace(0,chi_thres,100)**2), '--k', label='$\chi = %.2f$' %chi_thres)
    plt.colorbar(label='L (km)')
    plt.legend(loc='upper right')
    plt.xlabel('Elongation',fontsize=16)
    plt.ylabel('Planarity',fontsize=16)
    plt.title('N=%i Subsets: Hour %i' %(n,hour), fontsize=20)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.savefig('figures/Ngon_properties/HS_hour%i_%igon_shapes.png' %(hour,n), format='png',dpi = 300) 
    
    # %% Use the mean error and prob aliasing equations to compute the effective 
    # error for each subset where chi < 1
    L_save_N.append(L_save)
    Errors = np.zeros([N_tetra,pts]) + Background_Error
    Ef_Error = np.zeros([N_tetra,pts]) + Background_Error
    P_alias = np.zeros([N_tetra,pts])
    std = np.zeros([N_tetra,pts])
    chi_save = np.sqrt(np.power(E_save,2) + np.power(P_save,2))
    chi_save_N.append(chi_save)
    ind_save = []
    for j in range(N_tetra):
        chi = chi_save[j]
        if chi < chi_thres:
            ind_save.append(j)
            # k_s is the 'normalized' wavevector magnitude. This must be done because 
            # we did all of our learning on simulations where L=1
            k_s = k_mags*L_save[j]                  
            Errors[j,:] = Error_fun(k_s, chi, A)
            std[j,:] = std_eqn([k_s,chi], c0, c1, c2)
            #b = b_terms[int(n-4),:]
            P_alias[j,:] = Prob_alias(k_s, chi, b)
            
            # Compute expected value of error
            Ef_Error[j,:] = Errors[j,:]*(1-P_alias[j,:]) + np.maximum(100,400*np.power(k_s,-1))*P_alias[j,:]
    
    # because we only fit the functional form to chi < 1, we have to throw away some of the subsets
    # of spacecarft that do not fit this criteria
    print('   N=%i: %i/%i tetra have chi < %.2f' %(n,np.sum(chi_save<chi_thres),N_tetra,chi_thres))

    # %% make plot of the mean errors for all subsets with N=n
    
    plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    for j in ind_save:
        plt.plot(k_mags, Ef_Error[j,:], alpha=0.5, label=r'$\chi=%.2f$, $L=%.1f$ km' %(chi_save[j], L_save[j]))
        #plt.fill_between(k_mags, np.power(10,np.log10(Ef_Error[j,:]) - std[j,:]), np.power(10,np.log10(Ef_Error[j,:]) + std[j,:]), alpha=0.2)
    
    plt.xscale('log')
    plt.xlabel('$|k|$ ($km^{-1}$)',fontsize=18)
    plt.ylabel('$Error$ ($\%$)',fontsize=18)
    plt.title('$N=%i$' %(n), fontsize=20)
    #plt.legend(fontsize=10)
    plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.2)
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.1)
    plt.ylim(10**(-1),1*10**(2))
    plt.tight_layout()
    plt.savefig('figures/Ngon_errors/HS_hour%i_n%i.png' %(hour,n), format='png',dpi = 300) 
    
    # %% seach through all subsets with n=N and find the one which minimizes the 
    # expression [mu + std_coef*sigma]
    x = k_mags
    y = np.power(10,np.log10(Ef_Error) + std_coef*std)
    ind_min = np.argmin(y,axis=0)
    ind_save_N[int(n-4),:] = ind_min
    
    Er_min = np.zeros(pts)
    Std_min = np.zeros(pts)
    for k in range(pts):
        Er_min[k] = Ef_Error[ind_min[k],k]
        Std_min[k] = std[ind_min[k],k]
    
    Er_min_N[int(n-4),:] = Er_min
    Std_min_N[int(n-4),:] = Std_min
    # %% plot the minimized error curve for subsets with N=n
    
    plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(k_mags, Er_min, alpha=0.8, c='k')
    plt.fill_between(k_mags, np.power(10,np.log10(Er_min) - Std_min), np.power(10,np.log10(Er_min) + Std_min), alpha=0.1, color='k')
    plt.xscale('log')
    plt.xlabel('$|k|$ ($km^{-1}$)',fontsize=18)
    plt.ylabel('$Error$ ($\%$)',fontsize=18)
    plt.title('Min Error: $N=%i$' %(n), fontsize=20)
    #plt.legend(fontsize=10)
    plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.2)
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.1)
    plt.ylim(10**(-1),1*10**(2))
    plt.tight_layout()
    plt.savefig('figures/Ngon_min_errors/HS_hour%i_n%i_min.png' %(hour,n), format='png',dpi = 300) 
    
# %% search the saved minimums for each N and pick out best for each k_mag
N_ind = np.argmin(Er_min_N + std_coef*Std_min_N,axis=0)
N_min = N_ind + np.min(N_save)

Er_min = np.zeros(pts)
Std_min = np.zeros(pts)
Chi_min = np.zeros(pts)
L_min = np.zeros(pts)
for k in range(pts):
    sc_ind = ind_save_N[N_ind[k],k]
    Er_min[k] = Er_min_N[N_ind[k],k]
    Std_min[k] = Std_min_N[N_ind[k],k]
    Chi_min[k] = chi_save_N[N_ind[k]][sc_ind]
    L_min[k] = L_save_N[N_ind[k]][sc_ind]

# %% find region where error is less than Er_thres %
good_ind = np.where(np.power(10,np.log10(Er_min) + std_coef*Std_min) < Er_thres)[0]
if np.size(good_ind) != 0:
    k_min, k_max = k_mags[good_ind[0]], k_mags[good_ind[-1]]

# print region of accuracy and how many orders of k magnitude it is
print('--------------------------------------------------')
print('mu_error + %.1f sigma < %.1f %% for k in (%.2e, %.2e) 1/km' %(std_coef,Er_thres,np.mean(k_min),np.mean(k_max)))
if k_min != 0:
    print('   %.3f Orders of k Magnitude' %np.log10( np.nanmean(k_max/k_min)) )
else:
    print('   0 Orders of k Magnitude' ) 

# %% make figure of minimized Error showing selection 

fig, ax = plt.subplots(2, 2, sharex=True, figsize=(10, 6))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
ax[0,0].set_title('$(k_{min},k_{max})$ = (%.2e, %.2e) km' %(k_min,k_max), fontsize=16)
ax[0,0].plot(k_mags, Er_min, alpha=0.8, c='k', label='$\mu_{Error}$')
ax[0,0].fill_between(k_mags, np.power(10,np.log10(Er_min) - std_coef*Std_min), np.power(10,np.log10(Er_min) + std_coef*Std_min), alpha=0.1, color='k', label='$\mu \pm %.1f\sigma$' %std_coef)
ax[0,0].plot(k_mags, Er_thres*np.ones(np.shape(k_mags)), '--r', linewidth=1, label='%i\%% Error' %Er_thres)
ax[0,0].set_xscale('log')
#ax[0,0].set_yscale('log')
ax[0,0].set_ylabel('$Error$ ($\%$)',fontsize=18)
ax[0,0].grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.2)
ax[0,0].minorticks_on()
ax[0,0].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.1)
ax[0,0].set_ylim(0,1*10**(2))
ax[0,0].legend(fontsize=10,loc='upper left')

ax[1,0].scatter(k_mags, N_min, c='b', s=8)
ax[1,0].set_ylabel('$N$',fontsize=18)
ax[1,0].set_xlabel('$|k|$ ($km^{-1}$)',fontsize=18)
ax[1,0].grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.2)
#ax[1,0].minorticks_on()
ax[1,0].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.1)
ax[1,0].set_ylim(3.8,9.2)

ax[0,1].scatter(k_mags, L_min, c='b', s=8)
#ax[0,1].set_yscale('log')
ax[0,1].set_ylabel('$L$ (km)',fontsize=18)
ax[0,1].grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.2)
ax[0,1].minorticks_on()
ax[0,1].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.1)
ax[0,1].set_ylim(0,np.maximum(2200, np.max(L_min)))

ax[1,1].scatter(k_mags, Chi_min, c='b', s=8)
ax[1,1].set_ylabel('$\chi=\sqrt{E^2 + P^2}$',fontsize=15)
ax[1,1].set_ylim(0,1)
ax[1,1].set_xlabel('$|k|$ ($km^{-1}$)',fontsize=18)
ax[1,1].grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.2)
ax[1,1].minorticks_on()
ax[1,1].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.1)

fig.suptitle('HS Config: Hour %i' %(hour), fontsize=24)
plt.tight_layout()
plt.savefig('figures/Optimal_Ngon/HS_hour%i.png' %(hour), format='png',dpi = 300) 
