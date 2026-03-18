# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 10:06:52 2025

@author: hmaxwell
"""

# Problem 1a
import numpy as np
import matplotlib.pyplot as plt

V = np.linspace(-40, 140, 2000)

import numpy as np

# FUnction to Calulate the Rates
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


a_n = []
b_n = []
a_m = []
b_m = []
a_h = []
b_h = []

# For Loop to Creat Rates arrays
for v in V:
    soln = rates(v)
    a_n.append(soln[0])
    b_n.append(soln[1])
    a_m.append(soln[2])
    b_m.append(soln[3])
    a_h.append(soln[4])
    b_h.append(soln[5])

# Plotting the Rates as functions of Voltage
plt.plot(V, a_n, label = 'alpha_n')
plt.plot(V, b_n, label = 'beta_n')
plt.title('n Gate rate Constants', fontsize = 16)
plt.xlabel('Voltage (mV)', fontsize = 14)
plt.ylabel('Rate Constant (1/msec)', fontsize = 14)
plt.xlim(-40,140)
plt.grid(True)
plt.legend()
plt.show()
plt.plot(V, a_m, label = 'alpha_m')
plt.plot(V, b_m, label = 'beta_m')
plt.title('m Gate rate Constants',fontsize = 16)
plt.xlabel('Voltage (mV)',fontsize = 14)
plt.ylabel('Rate Constant (1/msec)',fontsize = 14)
plt.grid(True)
plt.xlim(-40,140)
plt.legend()
plt.show()
plt.plot(V, a_h, label = 'alpha_h')
plt.plot(V, b_h, label = 'beta_h')
plt.title('h Gate rate Constants',fontsize = 16)
plt.xlabel('Voltage (mV)',fontsize = 14)
plt.ylabel('Rate Constant (1/msec)',fontsize = 14)
plt.grid(True)
plt.xlim(-40,140)
plt.legend()

