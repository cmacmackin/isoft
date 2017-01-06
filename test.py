#!/usr/bin/env python

import scipy.optimize
import scipy.sparse.linalg
import numpy as np

points = 50
n_tot = 3*points
lower_bound = -1.0
upper_bound = 0.0
u_g = 2.0
q_g = 8.0
tol = 1e-8
perturbation = 0.1

xvals, dx = np.linspace(lower_bound,upper_bound,points,True,True)

expected = np.empty(n_tot)
expected[0:points]          = xvals - lower_bound
expected[points:2*points]   = q_g**(1./3.)
expected[2*points:3*points] = q_g**(2./3.)

def grad(arr, delta):
    # Computes upwinded gradient
    g = np.empty_like(arr)
    g[1:] = (arr[1:] - arr[:-1])/delta
    g[0] = 0.0
    return g

def residual(guess):
    assert guess.size%3 == 0
    n = guess.size/3
    resid = np.empty_like(guess)
    resid[0:n]     = np.gradient(guess[0:n]*guess[n:2*n],dx) - guess[n:2*n]
    resid[n:2*n]   = np.gradient(guess[0:n]*guess[n:2*n]**2,dx) - guess[2*n:3*n]
    resid[2*n:3*n] = np.gradient(guess[n:2*n]*guess[2*n:3*n], dx)
    resid[0]   = guess[0]*guess[2*n]
    resid[n]   = guess[n] - u_g
    resid[2*n] = guess[2*n] - q_g/u_g
    #print resid
    return resid

def residual_upwind(guess):
    assert guess.size%3 == 0
    n = guess.size/3
    resid = np.empty_like(guess)
    resid[0:n]     = grad(guess[0:n]*guess[n:2*n],dx) - guess[n:2*n]
    resid[n:2*n]   = grad(guess[0:n]*guess[n:2*n]**2,dx) - guess[2*n:3*n]
    resid[2*n:3*n] = grad(guess[n:2*n]*guess[2*n:3*n], dx)
    resid[0]   = guess[0]*guess[2*n]
    resid[n]   = guess[n] - u_g
    resid[2*n] = guess[2*n] - q_g/u_g
    #print resid
    return resid

assert np.all(np.abs(residual(expected)) < tol)

def precond(vec):
    assert vec.size%3 == 0
    assert vec.size == n_tot
    n = vec.size/3
    result = np.empty_like(vec)
    assert result.size == n_tot
    for i in range(1,n):
        result[i]     = np.trapz(vec[1:i],dx=dx)
        result[n+i]   = np.trapz(vec[n+1:n+i],dx=dx)
        result[2*n+i] = np.trapz(vec[2*n+1:3*n],dx=dx)
    return result

pre = scipy.sparse.linalg.LinearOperator((n_tot,n_tot), precond)

#initial = expected + perturbation*np.concatenate((np.cos(np.pi*2*xvals), 
#                                                  np.sin(np.pi*2*xvals),
#                                                  np.cos(np.pi*2*xvals)))
initial = expected + perturbation#np.random.normal(scale=perturbation,size=n_tot)
result = scipy.optimize.newton_krylov(residual, initial, verbose=True, 
                                      f_tol=tol, method='gmres')#, inner_M=pre)

try:
    assert np.all(np.abs(result - expected) < tol)
except AssertionError:
    print('Convergence to wrong solution')
    print(result-expected)
