# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 10:46:23 2025

@author: hmaxwell
"""

# Problem 2a
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

V = np.linspace(0, 120, 13)
t = np.linspace(0,11,2000)
def odefun (n,t,a_n,b_n):
    dndt = a_n*(1-n)-b_n*n
    return dndt
import numpy as np

def rates(Vm):
    # Compute alpha_n and beta_n
    if Vm == 10:  # usually the special case for n gating variable is at Vm = -10 mV
        a_n = 0.1
    else:
        a_n = (10 - Vm) / (100 * (np.exp((10 - Vm) / 10) - 1))
    b_n = 0.125 * np.exp(-Vm / 80)

    # Compute alpha_m and beta_m
    if Vm == 25:  # special case for m gating variable
        a_m = 1
    else:
        a_m = (25 - Vm) / (10 * (np.exp((25 - Vm) / 10) - 1))
    b_m = 4 * np.exp(-Vm / 18)

    # Compute alpha_h and beta_h
    a_h = 0.07 * np.exp(-Vm / 20)
    b_h = 1 / (np.exp((30 - Vm) / 10) + 1)

    return [a_n, b_n, a_m, b_m, a_h, b_h]

def SS (a_n,b_n,a_m,b_m,a_h,b_h):
    minf = a_m/(b_m+a_m)
    hinf = a_h/(b_h+a_h)
    ninf = a_n/(b_n+a_n)
    return [ninf,minf,hinf]
def tau (a_n,b_n,a_m,b_m,a_h,b_h):
    tau_n = 1/(a_n+b_n)
    tau_m = 1/(a_m+b_m)
    tau_h = 1/(a_h+b_h)
    return [tau_n,tau_m,tau_h]
def n (t, a_n, b_n):
    a_n_0,b_n_0 = rates(0)[0:2]
    n0 = a_n_0/(a_n_0+b_n_0)
    soln = odeint(odefun,n0,t, args = (a_n,b_n))
    return soln[:,0]
def gK (n):
    n = np.array(n)
    gk_bar = 24.01
    gK1 = gk_bar*(n**4.0)
    return gK1

plt.figure(figsize=(8,6))
for i, v in enumerate(V):               
    r = rates(v)
    n1 = n(t, r[0], r[1])
    gK1 = gK(n1)
    
    linestyle = '-' if i % 2 == 0 else '--'   # alternate by index
    plt.plot(t, gK1, linestyle=linestyle, label=f'V = {v:.0f} mV')

plt.xlabel('Time (ms)')
plt.ylabel('gK (mmol/cm²)')
plt.title('Potassium Conductance vs Time at Different Voltages')
plt.legend(title='Voltage', ncol=2, fontsize='small')
plt.tight_layout()
plt.show()

# Problem 2b
gk_final = []
for v in V:
    r = rates(v)
    n1 = n(t, r[0], r[1])
    gK1 = gK(n1)
    gk_final.append(gK1[-1])
    
gK_final = np.array(gk_final)
gK_norm = gk_final/np.max(gk_final)

plt.figure(figsize=(6,5))
plt.plot(V, gK_norm, 'o-', color='purple', lw=2)
plt.xlabel('Voltage (mV)')
plt.ylabel('Normalized gK (Activation)')
plt.title('K⁺ Channel Activation Curve (g–V)')
plt.grid(True)
plt.tight_layout()
plt.show()

# Problem 2c
tNa = np.linspace(0,100,2000)

def odefun2 (m,tNa,a_m,b_m):
    dmdt = a_m*(1-m)-b_m*m
    return dmdt
def odefun3(h,tNa,a_h,b_h):
    dhdt = a_h*(1-h)-b_h*h
    return dhdt
def m (tNa, a_m, b_m):
    a_m_0,b_m_0 = rates(0)[2:4]
    m0 = a_m_0/(a_m_0+b_m_0)
    soln = odeint(odefun2,m0,tNa, args = (a_m,b_m))
    return soln[:,0]
def h (tNa,a_h,b_h):
    a_h_0,b_h_0 = rates(0)[4:6]
    h0 = a_h_0/(a_h_0+b_h_0)
    soln = odeint(odefun3,h0,tNa, args = (a_h,b_h))
    return soln[:,0]
def gNa (m,h):
    m = np.array(m)
    h = np.array(h)
    gNa_bar = 70.7 
    gNa1 = gNa_bar*(m**3)*h
    return gNa1

plt.figure(figsize=(8,6))
for v in V:
    r = rates(v)
    m1 = m(tNa, r[2], r[3])
    h1 = h(tNa, r[4], r[5])
    gNa1 = gNa(m1, h1)
    plt.plot(tNa, gNa1, label=f'{v:.0f} mV')

plt.xlabel('Time (ms)')
plt.ylabel('gNa (mS/cm²)')
plt.title('Sodium Conductance vs Time for Voltage Steps')
plt.legend(ncol=2, fontsize='small')
plt.tight_layout()
plt.show()

# Na Activation Curve
gNa_final = []
for v in V:
    r = rates(v)
    m1 = m(tNa, r[2], r[3])
    h1 = h(tNa, r[4], r[5])
    gNa1 = gNa(m1, h1)
    gNa_final.append(gNa1[-1])

    
gNa_final = np.array(gNa_final)
gNa_norm = gNa_final/np.max(gNa_final)

plt.figure(figsize=(6,5))
plt.plot(V, gNa_norm, 'o-', color='purple', lw=2)
plt.xlabel('Voltage (mV)')
plt.ylabel('Normalized gNa (Activation)')
plt.title('Na⁺ Channel Activation Curve (g–V)')
plt.grid(True)
plt.tight_layout()
plt.show()






