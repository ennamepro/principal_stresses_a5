#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on Wed Mar 20 16:40:02 2024

@author: mjb
"""

# import math library that provides access to the mathematical functions defined by the C standard
# import numpy library that provides checkin calculation

import math
import numpy as np

# definition of stress tensor in 3d space
# fill values manually

sigma_x = 50; tau_xy = 20; tau_xz = -30
tau_yx = tau_xy; sigma_y = 40; tau_yz = 20
tau_zx = tau_xz; tau_zy = tau_yz; sigma_z = 20

# print stress tensor in 3d space

print("")
print("stress tensor:")
print(f"sigma_x={sigma_x:7.2f} \ttau_xy={tau_xy:7.2f} \ttau_xz={tau_xz:7.2f} \
          \ntau_yx={tau_yx:7.2f} \tsigma_y={sigma_y:7.2f} \ttau_yz={tau_yz:7.2f} \
          \ntau_zx={tau_zx:7.2f} \ttau_zy={tau_zy:7.2f} \tsigma_z={sigma_z:7.2f}\n")

# calculations according to Appendix B of Advanced Mechanics of Materials and Applied Elasticity, 5th edition
# invariants of the roots of stress cubic equation

I_1 = sigma_x + sigma_y + sigma_z
I_2 = sigma_x * sigma_y + sigma_x * sigma_z + sigma_y * sigma_z - tau_xy**2 - tau_yz**2 - tau_xz**2
I_3 = sigma_x * sigma_y * sigma_z + 2 * tau_xy * tau_yz * tau_xz - sigma_x * tau_yz**2 - sigma_y * tau_xz**2 - sigma_z * tau_xy**2

print("calculations using pure Python with math library\n")
print("stress tensor invariants: ")
print(f"I_1={I_1:7.2f} I_2={I_2:7.2f} I_3={I_3:7.2f}\n")

# calculation of constants

r_1 = I_1**2 / 3. - I_2
s_1 = math.sqrt(r_1 / 3.)
t_1 = math.sqrt(r_1**3 / 27.)
q_1 = I_1 * I_2 / 3. - I_3 - (2 / 27.) * I_1**3
alpha = math.acos(-q_1 / (2. * t_1))

# the principal stresses

sigma_principal = [0, 0, 0]
sigma_principal[0] = 2 * s_1 * math.cos(alpha / 3.) + I_1 / 3.
sigma_principal[1] = 2 * s_1 * math.cos((alpha / 3.) + 120 * math.pi / 180.) + I_1 / 3.
sigma_principal[2] = 2 * s_1 * math.cos((alpha / 3.) + 240 * math.pi / 180.) + I_1 / 3.

# sort list of the principal stresses
# print of the pricipal stress in 3d space

sigma_principal.sort(reverse=True)

print("principle stress values:")
print(f"sigma_1={sigma_principal[0]:8.3f} \tsigma_2={sigma_principal[1]:8.3f} \tsigma_3={sigma_principal[2]:8.3f}\n")

# the determinant cofactors of calculation matrices

a_1 = [0, 0, 0]
b_1 = [0, 0, 0]
c_1 = [0, 0, 0]
k_1 = [0, 0, 0]

# direction cosines matrices

l_1 = [0, 0, 0]
m_1 = [0, 0, 0]
n_1 = [0, 0, 0]

# calculation of direction cosines

for i in range(3):
    a_1[i] = ((sigma_y - sigma_principal[i]) * (sigma_z - sigma_principal[i]) - tau_yz**2)
    b_1[i] = -(tau_xy * (sigma_z - sigma_principal[i]) - (tau_xz * tau_yz))
    c_1[i] = ((tau_xy * tau_yz) - tau_xz * (sigma_y - sigma_principal[i]))
    k_1[i] = 1. / (math.sqrt((a_1[i]**2) + b_1[i]**2 + c_1[i]**2))
    l_1[i] = a_1[i] * k_1[i]
    m_1[i] = b_1[i] * k_1[i]
    n_1[i] = c_1[i] * k_1[i]
    
# print of direction cosines in 3d space

print("direction cosines:")
for i in range(3):
    print(f"l_1[{i+1}]={l_1[i]:7.4f} \tm_1[{i+1}]={m_1[i]:7.4f} \tn_1[{i+1}]={n_1[i]:7.4f}")
    
# calculations with NumPy library
# stress tensor definition

stress_tensor_1 = np.array([[sigma_x, tau_xy, tau_xz],[tau_yx, sigma_y, tau_yz],[tau_zx, tau_zy, sigma_z]])

print()
print("calculations using Python with NumPy library\n")
print("stress_tensor_1:\n", stress_tensor_1,"\n")

# stress tensor determinant

print("stress_tensor determinant:\n", np.linalg.det(stress_tensor_1),"\n")

# eigenvalues/principal and eigenvectors/direction cosines 

eigen_values, eigen_vectors = np.linalg.eig(stress_tensor_1)
eigen_values = np.sort(eigen_values)
print("eigenvalues:\n", np.flip(eigen_values),"\n")
print("eigenvectors:\n", eigen_vectors)