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
# Function to return dn/dt
def odefun (n,t,a_n,b_n):
    dndt = a_n*(1-n)-b_n*n
    return dndt

# Function to Calculate Rates
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
# Function to Calculate SS for n,m, and h
def SS (a_n,b_n,a_m,b_m,a_h,b_h):
    minf = a_m/(b_m+a_m)
    hinf = a_h/(b_h+a_h)
    ninf = a_n/(b_n+a_n)
    return [ninf,minf,hinf]
# Calculate tau for n,m and h
def tau (a_n,b_n,a_m,b_m,a_h,b_h):
    tau_n = 1/(a_n+b_n)
    tau_m = 1/(a_m+b_m)
    tau_h = 1/(a_h+b_h)
    return [tau_n,tau_m,tau_h]
# Solve for n using odeint
def n (t, a_n, b_n):
    a_n_0,b_n_0 = rates(0)[0:2]
    n0 = a_n_0/(a_n_0+b_n_0)
    soln = odeint(odefun,n0,t, args = (a_n,b_n))
    return soln[:,0]
# Solve for gK 
def gK (n):
    n = np.array(n)
    gk_bar = 24.01
    gK1 = gk_bar*(n**4.0)
    return gK1
# Plot gK for every voltage step; Every other line will be dotted
plt.figure(figsize=(8,6))
for i, v in enumerate(V):               
    r = rates(v)
    n1 = n(t, r[0], r[1])
    gK1 = gK(n1)
    
    linestyle = '-' if i % 2 == 0 else '--'   # alternate by index
    plt.plot(t, gK1, linestyle=linestyle, label=f'V = {v:.0f} mV')

plt.xlabel('Time (ms)',fontsize = 14)
plt.ylabel('gK (mmol/cm²)',fontsize = 14)
plt.title('Potassium Conductance vs Time at Different Voltages',fontsize = 16)
plt.legend(title='Voltage', ncol=2, fontsize='small')
plt.tight_layout()
plt.show()

# Problem 2b
gk_final = []
# Solve for n and gK
for v in V:
    r = rates(v)
    n1 = n(t, r[0], r[1])
    gK1 = gK(n1)
    gk_final.append(gK1[-1])
    
gK_final = np.array(gk_final)
# Normalize gK for Plotting
gK_norm = gk_final/np.max(gk_final)
# Plot Potassium Activation Curve
plt.figure(figsize=(6,5))
plt.plot(V, gK_norm, 'o-', color='purple', lw=2)
plt.xlabel('Voltage (mV)',fontsize = 14)
plt.ylabel('Normalized gK (Activation)',fontsize = 14)
plt.title('Potassium Channel Activation Curve',fontsize = 16)
plt.grid(True)
plt.tight_layout()
plt.show()

# Problem 2c
tNa = np.linspace(0,100,2000)
# Solve for dm/dt
def odefun2 (m,tNa,a_m,b_m):
    dmdt = a_m*(1-m)-b_m*m
    return dmdt
# Solve for dh/dt
def odefun3(h,tNa,a_h,b_h):
    dhdt = a_h*(1-h)-b_h*h
    return dhdt
#solve for m
def m (tNa, a_m, b_m):
    a_m_0,b_m_0 = rates(0)[2:4]
    m0 = a_m_0/(a_m_0+b_m_0)
    soln = odeint(odefun2,m0,tNa, args = (a_m,b_m))
    return soln[:,0]
# Solve for h
def h (tNa,a_h,b_h):
    a_h_0,b_h_0 = rates(0)[4:6]
    h0 = a_h_0/(a_h_0+b_h_0)
    soln = odeint(odefun3,h0,tNa, args = (a_h,b_h))
    return soln[:,0]
# Solve for gNa
def gNa (m,h):
    m = np.array(m)
    h = np.array(h)
    gNa_bar = 70.7 
    gNa1 = gNa_bar*(m**3)*h
    return gNa1
# Solve for gNa and plot for each voltage step
plt.figure(figsize=(8,6))
for v in V:
    r = rates(v)
    m1 = m(tNa, r[2], r[3])
    h1 = h(tNa, r[4], r[5])
    gNa1 = gNa(m1, h1)
    plt.plot(tNa, gNa1, label=f'{v:.0f} mV')

plt.xlabel('Time (ms)',fontsize = 14)
plt.ylabel('gNa (mmol/cm²)',fontsize = 14)
plt.title('Sodium Conductance vs Time for Voltage Steps',fontsize = 16)
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
# Normalize gNa for plotting
gNa_norm = gNa_final/np.max(gNa_final)

# Plot Na activation curve
plt.figure()
plt.plot(V, gNa_norm, 'o-', color='purple', lw=2)
plt.xlabel('Voltage (mV)',fontsize = 14)
plt.ylabel('Normalized gNa (Activation)',fontsize = 14)
plt.title('Sodium Channel Activation Curve',fontsize = 16)
plt.grid(True)
plt.tight_layout()
plt.show()


# Problem 2d

V_inac = np.linspace(-40, 80, 13)  # prepulse voltages

# Protocol voltages
V_hold = -80.0    # holding potential (mV) before prepulse
V_test = 60.0     # test pulse (mV)

# Durations and time vectors (ms)
t_hold = np.linspace(0, 10, 200)    # 10 ms hold at V_hold to produce a resting onset
t_pulse = np.linspace(0, 200, 2000) # 200 ms prepulse
t_test  = np.linspace(0, 20, 400)   # 20 ms test pulse

E_Na = 110.0
gNa_bar = 70.7

INa_Peak = []

plt.figure(figsize=(9,6))

for v in V_inac:
    # --- 1) Start from hold at V_hold to ensure channels are at resting state and produce an onset at prepulse ---
    a_m_hold, b_m_hold, a_h_hold, b_h_hold = rates(V_hold)[2:6]
    m0_hold = a_m_hold / (a_m_hold + b_m_hold)
    h0_hold = a_h_hold / (a_h_hold + b_h_hold)

    # integrate over hold to get initial condition at the moment of prepulse
    m_hold = odeint(odefun2, m0_hold, t_hold, args=(a_m_hold, b_m_hold))[:, 0]
    h_hold = odeint(odefun3, h0_hold, t_hold, args=(a_h_hold, b_h_hold))[:, 0]
    INa_hold = gNa_bar * (m_hold**3) * h_hold * (V_hold - E_Na)

    # --- 2) Prepulse to voltage v (this will produce the first transient peak at prepulse onset) ---
    a_m_pre, b_m_pre, a_h_pre, b_h_pre = rates(v)[2:6]
    m0_prep = m_hold[-1]
    h0_prep = h_hold[-1]
    m_prep = odeint(odefun2, m0_prep, t_pulse, args=(a_m_pre, b_m_pre))[:, 0]
    h_prep = odeint(odefun3, h0_prep, t_pulse, args=(a_h_pre, b_h_pre))[:, 0]
    INa_prep = gNa_bar * (m_prep**3) * h_prep * (v - E_Na)

    # --- 3) Step to test pulse (V_test) starting from end of prepulse (this gives the second transient) ---
    a_m_test, b_m_test, a_h_test, b_h_test = rates(V_test)[2:6]
    m0_test = m_prep[-1]
    h0_test = h_prep[-1]
    m_test = odeint(odefun2, m0_test, t_test, args=(a_m_test, b_m_test))[:, 0]
    h_test = odeint(odefun3, h0_test, t_test, args=(a_h_test, b_h_test))[:, 0]
    INa_test = gNa_bar * (m_test**3) * h_test * (V_test - E_Na)

    #  combine time and currents  for plotting
    t_combined = np.concatenate([
        t_hold,
        t_hold[-1] + t_pulse,
        t_hold[-1] + t_pulse[-1] + t_test
    ])
    INa_combined = np.concatenate([INa_hold, INa_prep, INa_test])

    # record the test-peak magnitude for the inactivation curve 
    INa_Peak.append(np.max(np.abs(INa_test)))

    # plot the full trace with label
    plt.plot(t_combined, INa_combined, label=f'{v:.0f} mV')


# plot layout
plt.xlim(-1, t_combined[-1] + 1)
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('I_Na (nA)', fontsize=14)
plt.title('Sodium Current Traces for Innactivation ', fontsize=16)
plt.legend(ncol=2, fontsize='x-small', loc='center')
plt.grid(True)
plt.tight_layout()
plt.show()

# Inactivation curve
INa_Peak = np.array(INa_Peak)
INa_norm = INa_Peak / np.max(INa_Peak)

plt.figure(figsize=(5,4))
plt.plot(V_inac, INa_norm, 'o-', linewidth=2)
plt.xlabel("Prepulse Voltage (mV)", fontsize=14)
plt.ylabel("Normalized Peak I_Na (test)", fontsize=14)
plt.title("Sodium Inactivation Curve", fontsize=16)
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
# Function to solve for m and h under the pulse protocol
def run_pulse(Vm, tspan, m0, h0):
    a_m, b_m, a_h, b_h = rates(Vm)[2:6]
    m_sol = odeint(odefun2, m0, tspan, args=(a_m, b_m))[:, 0]
    h_sol = odeint(odefun3, h0, tspan, args=(a_h, b_h))[:, 0]
    return m_sol, h_sol

INa_peaks = []

# Loop over durations and solve for needed values of n,m and h
for t2 in t2_values:
    # First pulse (to inactivate)
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

# Plot single trace for t2 = 100 ms
t2_example = 100

#First pulse (V1 = 60 mV, 200 ms)
tspan1 = np.linspace(0, t1, 1000)
a_m1, b_m1, a_h1, b_h1 = rates(V1)[2:6]
m_inf1 = a_m1 / (a_m1 + b_m1)
h_inf1 = a_h1 / (a_h1 + b_h1)
m1, h1 = run_pulse(V1, tspan1, m_inf1, h_inf1)
INa1 = gNa_bar * (m1**3) * h1 * (V1 - E_Na)

# Interim (V2 = 0 mV, variable duration
tspan2 = np.linspace(0, t2_example, max(2, int(t2_example * 10 + 2)))
m2, h2 = run_pulse(V2, tspan2, m1[-1], h1[-1])
INa2 = gNa_bar * (m2**3) * h2 * (V2 - E_Na)   # basically ~0 nA

#Second pulse (V3 = 60 mV, 20 ms)
tspan3 = np.linspace(0, t3, 1000)
m3, h3 = run_pulse(V3, tspan3, m2[-1], h2[-1])
INa3 = gNa_bar * (m3**3) * h3 * (V3 - E_Na)

# Concatenate times and currents
t_full = np.concatenate([
    tspan1,
    t1 + tspan2,
    t1 + t2_example + tspan3
])
INa_full = np.concatenate([
    INa1,
    INa2,
    INa3
])


#Plot sodium current Recovery
plt.figure(figsize=(8,5))
plt.plot(t_full, INa_full, color='navy')
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('I_Na (nA)', fontsize=14)
plt.title(f'Sodium Current Recovery (Interim = {t2_example} ms)', fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.show()
# Plot recovery curve 
plt.figure(figsize=(6,5))
plt.plot(t2_values, INa_norm, 'o', label='Data', markersize = 4)
plt.plot(t2_values, exp_recovery(t2_values, *popt), '-', label=f'Fit: τ = {tau_fit:.2f} ms')
plt.xlabel('Interim Duration (ms)',fontsize = 14)
plt.ylabel('Normalized Peak I_Na',fontsize = 14)
plt.title('Sodium Channel Recovery from Inactivation',fontsize = 16)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()



























