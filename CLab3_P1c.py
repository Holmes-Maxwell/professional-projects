# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 10:29:53 2025

@author: hmaxwell
"""
# Problem 1c
import numpy as np
import matplotlib.pyplot as plt

V = np.linspace(-40, 140, 2000)

import numpy as np

# Function to Calculate Rates
def rates(Vm):
    # Compute alpha_n and beta_n
    if Vm == -10:  # usually the special case for n gating variable is at Vm = -10 mV
        a_n = 0.1
    else:
        a_n = (10 - Vm) / (100 * (np.exp((10 - Vm) / 10) - 1))
    b_n = 0.125 * np.exp(-Vm / 80)

    # Compute alpha_m and beta_m
    if Vm == -25:  # special case for m gating variable
        a_m = 1
    else:
        a_m = (25 - Vm) / (10 * (np.exp((25 - Vm) / 10) - 1))
    b_m = 4 * np.exp(-Vm / 18)

    # Compute alpha_h and beta_h
    a_h = 0.07 * np.exp(-Vm / 20)
    b_h = 1 / (np.exp((30 - Vm) / 10) + 1)

    return [a_n, b_n, a_m, b_m, a_h, b_h]


ninf = []
minf = []
hinf = []
tau_n = []
tau_m = []
tau_h = []

   
   # Function to Calulate Steady State Gating Values   
def SS (a_n,b_n,a_m,b_m,a_h,b_h):
    minf = a_m/(b_m+a_m)
    hinf = a_h/(b_h+a_h)
    ninf = a_n/(b_n+a_n)
    return [ninf,minf,hinf]
# Function to Caluclate Taus
def tau (a_n,b_n,a_m,b_m,a_h,b_h):
    tau_n = 1/(a_n+b_n)
    tau_m = 1/(a_m+b_m)
    tau_h = 1/(a_h+b_h)
    return [tau_n,tau_m,tau_h]
# Appending Arrays for Plotting
for v in V:
    soln = rates(v)
    a_n = soln[0]
    b_n = soln[1]
    a_m = soln[2]
    b_m = soln[3]
    a_h = soln[4]
    b_h = soln[5]
    
    soln1 = SS(a_n,b_n,a_m,b_m,a_h,b_h)
    soln2 = tau(a_n,b_n,a_m,b_m,a_h,b_h)
    ninf.append(soln1[0])
    minf.append(soln1[1])
    hinf.append(soln1[2])
    tau_n.append(soln2[0])
    tau_m.append(soln2[1])
    tau_h.append(soln2[2])
# Create Arrays and Calulcate m^3 and tau^3
minf = np.array(minf)
tau_m = np.array(tau_m)
minf_3 = (minf**3.0) 
taum_3 = (tau_m**3.0)

# Plot minf,minf^3, and h
plt.plot(V,minf, label = 'minf')
plt.plot(V,minf_3, label = 'minf^3')
plt.plot(V,hinf, label = 'hinf')
plt.title('Steady State Values for m, m^3 and h',fontsize = 16)
plt.xlabel('Voltage (mV)',fontsize = 14)
plt.ylabel('Steady State Values',fontsize = 14)
plt.grid(True)
plt.xlim(-40,140)
plt.legend()
plt.show()
# Plot taum,taum^3,and tauh
plt.plot(V,tau_m, label = 'tau_m')
plt.plot(V,taum_3, label = 'tau_m^3')
plt.plot(V,tau_h, label = 'tau_h')
plt.title('Values for tau_m, tau_m^3 and tau_h',fontsize = 16)
plt.xlabel('Voltage (mV)',fontsize = 14)
plt.ylabel('Steady State Values',fontsize = 14)
plt.grid(True)
plt.xlim(-40,140)
plt.legend()
plt.show()
