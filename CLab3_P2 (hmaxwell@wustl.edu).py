# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 10:46:23 2025

@author: hmaxwell
"""

# Problem 2a
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import curve_fit

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
plt.xlim(0, 10)
plt.grid(True)
plt.tight_layout()
plt.show()

# Na Activation Curve
gNa_final = []
for v in V:
    r = rates(v)
    m1 = m(tNa, r[2], r[3])
    gNa1 = m1**3
    gNa_final.append(gNa1[-1])

    
gNa_final = np.array(gNa_final)
gNa_norm = gNa_final/np.max(gNa_final)

plt.figure()
plt.plot(V, gNa_norm, 'o-', color='purple', lw=2)
plt.xlabel('Voltage (mV)')
plt.ylabel('Normalized gNa (Activation)')
plt.title('Na⁺ Channel Activation Curve (g–V)')
plt.grid(True)
plt.tight_layout()
plt.show()

# Problem 2D
V_inac = np.linspace(-40,80,13)
t_pulse = np.linspace(0,200,2000)
t_test = np.linspace(0,20,2000)
t_total = np.concatenate([t_pulse, t_pulse[-1] + t_test])  # total 220 ms

E_Na = 110
V_test = 60
INa_Peak = []
gNa_bar = 70.7

plt.figure(figsize=(8,6))
for v in V_inac:
    # Prepulse
    r = rates(v)
    a_h = r[4]
    b_h = r[5]
    a_h_0,b_h_0 = rates(0)[4:6]
    h0 = a_h_0/(a_h_0+b_h_0)
    h_pulse = odeint(odefun3,h0,t_pulse,args = (a_h,b_h))[:,0]

    # Test Pulse
    r_test = rates(V_test)
    a_m_test = r_test[2]
    b_m_test = r_test[3]
    a_h_test = r_test[4]
    b_h_test = r_test[5]
    m0_test = 0
    h0_test = h_pulse[-1]
    m_test = odeint(odefun2,m0_test,t_test,args = (a_m_test,b_m_test))[:,0]
    h_test = odeint(odefun3, h0_test, t_test, args=(a_h_test,b_h_test))[:,0]

    # Na Current
    INa = gNa_bar*(m_test**3)*h_test*(V_test-E_Na)
    INa_Peak.append(np.max(np.abs(INa)))
    plt.plot(t_test, INa, label=f'{v} mV')  # only test pulse

plt.xlabel('Time (ms)')
plt.ylabel('I_Na (µA/cm²)')
plt.title('Na⁺ Current during Test Pulse for Different Prepulses')
plt.legend(ncol=2, fontsize='x-small')
plt.grid(True)
plt.tight_layout()
plt.show()

# Na Innacticvation
INa_Peak = np.array(INa_Peak)
INa_norm = INa_Peak / np.max(INa_Peak)

plt.figure(figsize=(5,4))
plt.plot(V_inac, INa_norm, 'o-', linewidth=2)
plt.xlabel("Prepulse Voltage (mV)")
plt.ylabel("Normalized Peak I_Na")
plt.title("Na⁺ Inactivation Curve")
plt.grid(True)
plt.tight_layout()
plt.show()

# Problem 2e

V1 = 60
t1 = 200
V2 = 0
t2_values = np.linspace(0,1000,2000)
V3 = 60
t3 = 20
def run_pulse(Vm, tspan, m0, h0):
    a_m, b_m, a_h, b_h = rates(Vm)[2:6]
    m_sol = odeint(odefun2, m0, tspan, args=(a_m, b_m))[:, 0]
    h_sol = odeint(odefun3, h0, tspan, args=(a_h, b_h))[:, 0]
    return m_sol, h_sol

INa_peaks = []

# Loop over durations and solve for needed values of n,m and h
for t2 in t2_values:
    # 1️⃣ First pulse (to inactivate)
    tspan1 = np.linspace(0, t1, 1000)
    a_m1, b_m1, a_h1, b_h1 = rates(V1)[2:6]
    m_inf1 = a_m1 / (a_m1 + b_m1)
    h_inf1 = a_h1 / (a_h1 + b_h1)
    m0 = m_inf1  # assume steady state
    h0 = h_inf1

    m1, h1 = run_pulse(V1, tspan1, m0, h0)

    # Interim at 0 mV
    tspan2 = np.linspace(0, t2, max(2, int(t2 * 10 + 2)))  # variable resolution
    m2, h2 = run_pulse(V2, tspan2, m1[-1], h1[-1])

    # Second pulse at 60 mV
    tspan3 = np.linspace(0, t3, 1000)
    m3, h3 = run_pulse(V3, tspan3, m2[-1], h2[-1])

    # Sodium current during second pulse
    INa = gNa_bar * (m3**3) * h3 * (V3 - E_Na)
    INa_peaks.append(np.max(np.abs(INa)))

# Normalize recovery
INa_peaks = np.array(INa_peaks)
INa_norm = INa_peaks / np.max(INa_peaks)

# Exponential fit for recovery 
def exp_recovery(t, tau, A):
    return 1 - A * np.exp(-t / tau)

popt, _ = curve_fit(exp_recovery, t2_values, INa_norm, p0=[10, 1])
tau_fit = popt[0]

#  Plot single trace for t2 = 100 ms 
t2_example = 100
a_m1, b_m1, a_h1, b_h1 = rates(V1)[2:6]
m_inf1 = a_m1 / (a_m1 + b_m1)
h_inf1 = a_h1 / (a_h1 + b_h1)
m0 = m_inf1
h0 = h_inf1

tspan1 = np.linspace(0, t1, 1000)
m1, h1 = run_pulse(V1, tspan1, m0, h0)

tspan2 = np.linspace(0, t2_example, int(t2_example * 10 + 2))
m2, h2 = run_pulse(V2, tspan2, m1[-1], h1[-1])

tspan3 = np.linspace(0, t3, 1000)
m3, h3 = run_pulse(V3, tspan3, m2[-1], h2[-1])
INa = gNa_bar * (m3**3) * h3 * (V3 - E_Na)

plt.figure(figsize=(7,5))
plt.plot(tspan3, INa, color='navy')
plt.xlabel('Time (ms)')
plt.ylabel('I_Na (µA/cm²)')
plt.title('Na⁺ Current Trace (Second Pulse, 100 ms Interim)')
plt.grid(True)
plt.tight_layout()
plt.show()
# Plot recovery curve 
plt.figure(figsize=(6,5))
plt.plot(t2_values, INa_norm, 'o', label='Data', markersize = 4)
plt.plot(t2_values, exp_recovery(t2_values, *popt), '-', label=f'Fit: τ = {tau_fit:.2f} ms')
plt.xlabel('Interim Duration (ms)')
plt.ylabel('Normalized Peak I_Na')
plt.title('Na⁺ Channel Recovery from Inactivation')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()



























