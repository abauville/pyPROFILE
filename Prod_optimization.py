#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 11:23:03 2021

@author: abauville


Projet avec Dewi: optimization de la productivite en O2 des organismes en fonction de la profondeur
d'apres le papier de Berg et al 1998
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import scipy.optimize as opt
import pandas as pd

# User defined values
# =============================================
# Physical parameters
# ===================
phi_val = 0.8
alpha_val = 0.0
D = 11.7e-06
Ds_val = phi_val**2 * D
# Db_val = 0.0

# n = 10 # number of grid points
# d_top = -0.02 # surface depth
# d_bot = 0.27 # bottom depth
# x = np.linspace(d_top,d_bot,n)
# Dx = x[1] - x[0]


## Fake data
# =============================================
# df = pd.read_csv("./GT_test.csv")
df = pd.read_csv("./profil_D4_rep.csv", sep=';')
df = df[df['X']>=0.028]
x = df['X'].values
Dx = x[1] - x[0]
phi = df['PHI'].values
alpha = df['ALPHA'].values
Db = df['DB'].values
C_gt = df['C'].values

n = len(C_gt)

# C_gt = (x-0.15)**2
plt.figure(2)
plt.clf()
plt.plot(C_gt)


n_replica = df.shape[1]-4
C_GT = np.zeros((n_replica,df.shape[0]))
C_GT[0,:] = df['C'].values
C_GT[1,:] = df['C2'].values
C_GT[2,:] = df['C3'].values
# C_GT = 25.0*np.andom.randn(4,n)
# C_GT = C_GT+C_gt



# phi = phi_val * np.ones(n)
# alpha = alpha_val * np.ones(n)
Ds = Ds_val * np.ones(n)
# Db = Db_val * np.ones(n)






# Boundary conditions
# ===================
# Type of boundary conditions (1:t=C b=C, 2:t=C t=F, 3:b=C b=F 4:t=C b=F 5:t=F b=C)
BC_type = 1
BC_val1 = C_gt[0]
BC_val2 = C_gt[-1]







# Value to optimize

# R_val = 0.0
def compute_C(R_val):

    R = np.ones(n)
    Db_val = R_val[0]
    Db_val = 0
    Db = Db_val*np.ones(n)
    R_val = R_val[1:]
    Iprevious = 0
    if np.isscalar(R_val):
        R*=R_val
    else:
        for i in range(len(R_val)):
            I = int(np.ceil(n/(len(R_val))*(i+1)))
            R[Iprevious:I] = R_val[i]
            Iprevious = I
    # Initialize parameters array
    AA = np.ones(n)
    BB = np.ones(n)
    CC = np.ones(n)
    DD = np.zeros(n)


    P = phi*(Ds+Db)

    # Set BC
    if BC_type == 1:
        top_boundary_type = 0 # 0: given concentration, 1:given flux
        bot_boundary_type = 0 # 0: given concentration, 1:given flux
        top_boundary_val = BC_val1
        bot_boundary_val = BC_val2

    elif BC_type == 2:
        top_boundary_type = 0 # 0: given concentration, 1:given flux
        bot_boundary_type = 1 # 0: given concentration, 1:given flux
        top_boundary_val = BC_val1
        bot_boundary_val = np.sum(R*Dx) - BC_val2

    elif BC_type == 3:
        top_boundary_type = 1 # 0: given concentration, 1:given flux
        bot_boundary_type = 0 # 0: given concentration, 1:given flux

        top_boundary_val = np.sum(R*Dx) - BC_val2
        bot_boundary_val = BC_val1

    elif BC_type == 4:
        top_boundary_type = 0 # 0: given concentration, 1:given flux
        bot_boundary_type = 1 # 0: given concentration, 1:given flux

        top_boundary_val = BC_val1
        bot_boundary_val = BC_val2

    elif BC_type == 5:
        top_boundary_type = 1 # 0: given concentration, 1:given flux
        bot_boundary_type = 0 # 0: given concentration, 1:given flux

        top_boundary_val = BC_val1
        bot_boundary_val = BC_val2


    C0 = top_boundary_val # /!\ hard coded

    # System of equations


    for j in range(1,n-1):
        # Fill AA, BB, CC, DD

        AA[j] = 1.0/Dx * (2.0*P[j-1] * P[j  ]) / (Dx*P[j  ] + Dx*P[j-1])
        CC[j] = 1.0/Dx * (2.0*P[j  ] * P[j+1]) / (Dx*P[j+1] + Dx*P[j  ])
        BB[j] = -AA[j] - CC[j] - phi[j]*alpha[j]

        DD[j] = -phi[j] *alpha[j] *C0 - R[j]


    # Fill the system of equations
    A_mat = diags(AA[1:], -1, shape=(n,n)) + diags(BB, 0, shape=(n,n)) + diags(CC[:-1],1, shape=(n,n))
    RHS = DD.copy()

    # Apply BC

    if top_boundary_type==0:
        # Dirichlet top
        A_mat[0,0] = 1.0
        A_mat[0,1] = 0.0
        RHS[0] = top_boundary_val
    elif top_boundary_type==1:
        # Neumann top
        BB0 = (phi[1] * (Ds[1]+Db[1]) ) / (Dx/2)
        CC0 = - (phi[1] * (Ds[1]+Db[1])  ) / (Dx/2)
        RHS[0] = top_boundary_val
        A_mat[0,0] = BB0
        A_mat[0,1] = CC0
    else:
        raise ValueError('Unknwon top_boundary_type')

    if bot_boundary_type==0:
        # Dirichlet bottom
        A_mat[-1,-1] = 1.0
        A_mat[-1,-2] = 0.0
        RHS[-1] = bot_boundary_val
    elif bot_boundary_type==1:
        # Neumann top
        AA0 = -(phi[-2] * (Ds[-2]+Db[-2]) ) / (Dx/2)
        BB0 =  (phi[-2] * (Ds[-2]+Db[-2])  ) / (Dx/2)
        RHS[-1] = bot_boundary_val
        A_mat[-1,-1] = BB0
        A_mat[-1,-2] = AA0

    C = spsolve(A_mat, RHS)
    return C


def SSE(R_vals):
    SSE = 0
    for i in range(C_GT.shape[0]):
        gt = C_GT[i,:]
        SSE += np.sum((compute_C(R_vals)-gt)**2)
    return SSE

R_val = -1e-4 # initial guess



# R_vals_0 = np.array([0,0,0,0,0,0,0,0])
n_r = 5
R_vals_0 = np.zeros(n_r)
# opt_result = opt.minimize(SSE,R_vals_0, bounds=opt.Bounds(-100,0))
opt_result = opt.minimize(SSE,R_vals_0)
R_opt = opt_result.x
plt.figure(1)
plt.clf()
plt.subplot(211)

plt.plot(C_GT.T,'x-r',lw=0.3)

plt.plot(compute_C(R_opt),'ok',markerfacecolor='none',markersize=10)
plt.subplot(212)
plt.plot(R_opt)


# plt.title([f'{r:.2f}' for r in R_opt])
