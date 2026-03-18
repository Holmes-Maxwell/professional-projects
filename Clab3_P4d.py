# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 2025
@author: holme
"""

# Problem 4d — Effect of gK Reduction on Action Potential and Refractory Period
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# Parameters
gNa_bar = 70.7
gK_Norm = 24.31
gK_bar50 = gK_Norm * 0.5
gK_bar25 = gK_Norm * 0.75
ENa = 110
EK = -12
Cm = 1
gL_bar = 0.3
EL = 0
t = np.linspace(0, 50, 5000)


# Rate Functions
def rates(Vm):
    if Vm == -10:
        a_n = 0.1
    else:
        a_n = (10 - Vm) / (100 * (np.exp((10 - Vm) / 10) - 1))
    b_n = 0.125 * np.exp(-Vm / 80)

    if Vm == -25:
        a_m = 1
    else:
        a_m = (25 - Vm) / (10 * (np.exp((25 - Vm) / 10) - 1))
    b_m = 4 * np.exp(-Vm / 18)

    a_h = 0.07 * np.exp(-Vm / 20)
    b_h = 1 / (np.exp((30 - Vm) / 10) + 1)
    return [a_n, b_n, a_m, b_m, a_h, b_h]


# Stimulus
def IStim(t):
    return 10 if 0.5 <= t <= 2.5 else 0

# HH Equations
def derivatives(y, t, gK_bar, IStim):
    V, n, m, h = y
    a_n, b_n, a_m, b_m, a_h, b_h = rates(V)
    dndt = a_n * (1 - n) - b_n * n
    dmdt = a_m * (1 - m) - b_m * m
    dhdt = a_h * (1 - h) - b_h * h

    IK  = gK_bar  * n**4 * (V - EK)
    INa = gNa_bar * m**3 * h * (V - ENa)
    IL  = gL_bar  * (V - EL)

    I_ext = IStim(t)
    dVdt = (I_ext - (IK + INa + IL)) / Cm
    return [dVdt, dndt, dmdt, dhdt]


# Initial Conditions
a_n, b_n, a_m, b_m, a_h, b_h = rates(0)
n0 = a_n / (a_n + b_n)
m0 = a_m / (a_m + b_m)
h0 = a_h / (a_h + b_h)
y0 = [0, n0, m0, h0]


# Solve for Normal, 50%, and 25% gK
sol_norm = odeint(derivatives, y0, t, args=(gK_Norm, IStim))
sol50 = odeint(derivatives, y0, t, args=(gK_bar50, IStim))
sol25 = odeint(derivatives, y0, t, args=(gK_bar25, IStim))

Vn, nn, mn, hn = sol_norm.T
V50, n50, m50, h50 = sol50.T
V25, n25, m25, h25 = sol25.T


# Compute Ionic Currents
def compute_currents(V, n, m, h, gK_bar):
    IK  = gK_bar  * n**4 * (V - EK)
    INa = gNa_bar * m**3 * h * (V - ENa)
    IL  = gL_bar  * (V - EL)
    return IK, INa, IL

IKn, INan, ILn = compute_currents(Vn, nn, mn, hn, gK_Norm)
IK50, INa50, IL50 = compute_currents(V50, n50, m50, h50, gK_bar50)
IK25, INa25, IL25 = compute_currents(V25, n25, m25, h25, gK_bar25)


# Plot: Stimulus
I_values = np.array([IStim(tt) for tt in t])
plt.figure(figsize=(6,3))
plt.plot(t, I_values, color='black', lw=2)
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('I_stim (nA)', fontsize=14)
plt.title('Stimulus Current Pulse', fontsize=16)
plt.xlim(0, 20)
plt.grid(True)
plt.show()

# Plot: Ionic Currents
plt.figure(figsize=(8,4))
plt.plot(t, IKn, label='I_K (Normal)')
plt.plot(t, IK50, label='I_K (50%)')
plt.plot(t, IK25,'--',  label='I_K (25%)')
plt.plot(t, INa25,'--', label='I_Na (25%)')
plt.plot(t, INa50, label='I_Na (50%)')
plt.plot(t,INan, label = 'I_Na (Normal)')
plt.plot(t, IL25, '--', label='I_L (25%)')
plt.plot(t, IL50, label='I_L (50%)')
plt.plot(t, ILn, label='I_L (Normal)')
plt.title('Ionic Currents with Reduced gK', fontsize=16)
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('Current (nA)', fontsize=14)
plt.legend(loc = 'upper right')
plt.xlim(0, 20)
plt.grid(True)
plt.show()


# Plot: Action Potentials
plt.figure(figsize=(8,5))
plt.plot(t, Vn, label='Normal gK', color='blue', lw=2)
plt.plot(t, V50, label='50% gK', color='red', lw=2)
plt.plot(t, V25, label='25% gK', color='green', lw=2)
plt.title('Action Potentials with gK Reduction', fontsize=16)
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('Membrane Voltage (mV)', fontsize=14)
plt.legend()
plt.xlim(0, 20)
plt.grid(True)
plt.show()


# Plot: Gating Variables
plt.figure(figsize=(8,4))
plt.plot(t, n50, label='n (50%)')
plt.plot(t, n25, '--', label='n (25%)')
plt.plot(t,nn, label = 'n (Normal)')
plt.plot(t, m50, label='m (50%)')
plt.plot(t, m25, '--', label='m (25%)')
plt.plot(t,mn, label = 'm (Normal)')
plt.plot(t, h50, label='h (50%)')
plt.plot(t, h25, '--', label='h (25%)')
plt.plot(t,hn, label = 'h (Normal)')
plt.title('Gating Variables for gK Reduction', fontsize=16)
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('Gating Variable', fontsize=14)
plt.legend()
plt.xlim(0, 20)
plt.grid(True)
plt.show()


# Refractory Period Analysis
def IStim_double(t, dt, amp=10):
    if (0.5 <= t <= 2.5) or (0.5 + dt <= t <= 2.5 + dt):
        return amp
    else:
        return 0

def simulate_dt(dt, amp, gK_bar):
    t = np.linspace(0, 50 + dt, 6000)
    Is_func = lambda tt: IStim_double(tt, dt, amp)
    sol = odeint(derivatives, y0, t, args=(gK_bar, Is_func))
    V = sol[:, 0]
    spike_inds = np.where((V[1:] > 0) & (V[:-1] <= 0))[0]
    num_spikes = len(spike_inds)
    return num_spikes, V, t

dt_values = np.arange(1, 31, 1)
I1 = 10

#  50% gK 
gK_bar = gK_bar50
I2_ratio_reduced = []
for dt in dt_values:
    amp = I1
    spikes, _, _ = simulate_dt(dt, amp, gK_bar)
    if spikes < 2:
        amp_test = I1
        while amp_test < 400:
            spikes, _, _ = simulate_dt(dt, amp_test, gK_bar)
            if spikes >= 2:
                break
            amp_test += 2
        I2_ratio_reduced.append(amp_test / I1)
    else:
        I2_ratio_reduced.append(1.0)

#  25% gK 
gK_bar = gK_bar25
I2_ratio_25 = []
for dt in dt_values:
    amp = I1
    spikes, _, _ = simulate_dt(dt, amp, gK_bar)
    if spikes < 2:
        amp_test = I1
        while amp_test < 400:
            spikes, _, _ = simulate_dt(dt, amp_test, gK_bar)
            if spikes >= 2:
                break
            amp_test += 2
        I2_ratio_25.append(amp_test / I1)
    else:
        I2_ratio_25.append(1.0)

# Normal gK 
gK_bar = gK_Norm
I2_ratio_norm = []
for dt in dt_values:
    amp = I1
    spikes, _, _ = simulate_dt(dt, amp, gK_bar)
    if spikes < 2:
        amp_test = I1
        while amp_test < 400:
            spikes, _, _ = simulate_dt(dt, amp_test, gK_bar)
            if spikes >= 2:
                break
            amp_test += 2
        I2_ratio_norm.append(amp_test / I1)
    else:
        I2_ratio_norm.append(1.0)


# Plot: Refractory Period Curves
plt.figure(figsize=(7,5))
plt.plot(dt_values, I2_ratio_norm, label='Normal $g_K$', color='blue', lw=2)
plt.plot(dt_values, I2_ratio_reduced, label='50% $g_K$', color='red', lw=2)
plt.plot(dt_values, I2_ratio_25, label='25% $g_K$', color='green', lw=2)
plt.axhline(y=1, color='k', linestyle='--', linewidth=1)
plt.xlim(0, 30)
plt.ylim(0, 20)
plt.xlabel(r'$\Delta t$ (ms)', fontsize=14)
plt.ylabel(r'$\frac{I_2}{I_1}$', fontsize=14)
plt.title('Refractory Period Comparison: Normal vs Reduced $g_K$', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True)
plt.show()
