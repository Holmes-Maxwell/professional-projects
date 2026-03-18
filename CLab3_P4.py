# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 19:17:04 2025

@author: holme
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 18:42:13 2025

@author: holme
"""

# Problem 4a
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
# Parameters
gNa_Norm = 70.7
gNa_bar50 = 70.7 * 0.5   # 50% reduction
gNa_bar25 = 70.7 * 0.75  # 25% reduction
gK_bar = 24.31
ENa = 110
EK = -12
Cm = 1
gL_bar = 0.3
EL = 0
t = np.linspace(0, 50, 5000)

# Rate function
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
    return 25 if 0.5 <= t <= 2.5 else 0

# Hodgkin-Huxley equations 
def derivatives(y, t, gNa_bar, IStim):
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

# Initial conditions
a_n, b_n, a_m, b_m, a_h, b_h = rates(0)
n0 = a_n / (a_n + b_n)
m0 = a_m / (a_m + b_m)
h0 = a_h / (a_h + b_h)
y0 = [0, n0, m0, h0]

# Solve for both cases
soln = odeint(derivatives, y0, t, args = (gNa_Norm, IStim))
sol50 = odeint(derivatives, y0, t, args=(gNa_bar50, IStim))
sol25 = odeint(derivatives, y0, t, args=(gNa_bar25, IStim))

V50, n50, m50, h50 = sol50.T
V25, n25, m25, h25 = sol25.T
Vn, nn, mn, hn = soln.T

#Compute currents 
IK50  = gK_bar * n50**4 * (V50 - EK)
INa50 = gNa_bar50 * m50**3 * h50 * (V50 - ENa)
IL50  = gL_bar * (V50 - EL)

IK25  = gK_bar * n25**4 * (V25 - EK)
INa25 = gNa_bar25 * m25**3 * h25 * (V25 - ENa)
IL25  = gL_bar * (V25 - EL)

IKNorm = gK_bar * n25**4 * (V25 - EK)
INaNorm = gNa_Norm * m25**3 * h25 * (V25 - ENa)
ILNorm = gL_bar * (V25 - EL)

# Plot stimulus 
I_values = np.array([IStim(tt) for tt in t])
plt.figure(figsize=(6,3))
plt.plot(t, I_values, color='black', lw=2)
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('I_stim (nA)', fontsize=14)
plt.title('Stimulus Current Pulse', fontsize=16)
plt.xlim(0, 20)
plt.grid(True)
plt.show()

#Plot ionic currents 
plt.figure(figsize=(8,4))
plt.plot(t, INa50, label='I_Na (50% GNa)')
plt.plot(t, INa25, '--',label='I_Na (25% GNa)')
plt.plot(t, INaNorm, label='I_Na (Normal GNa)')
plt.plot(t, IK50, label='I_K (50% GNa)')
plt.plot(t, IK25,'--', label='I_K (25% GNa)')
plt.plot(t, IKNorm, label='I_K (Normal GNa)')
plt.plot(t, IL50, label='I_L (50% GNa)')
plt.plot(t, IL25,'--', label='I_L (25% GNa)')
plt.plot(t, ILNorm, label='I_L (Normal GNa)')
plt.title('Ionic Currents with 25% and 50% Reduction of GNa', fontsize=16)
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('Ionic Currents (nA)', fontsize=14)
plt.legend(loc = 'upper right')
plt.xlim(0, 20)
plt.grid(True)
plt.show()

# Plot action potentials 
plt.figure(figsize=(8,5))
plt.plot(t, Vn, label='Normal GNa', color='purple', lw=2)
plt.plot(t, V50, label='50% GNa', color='orange', lw=2)
plt.plot(t, V25, label='25% GNa', color='green', lw=2)
plt.title('Action Potentials with 25% and 50% Reduction of GNa', fontsize=16)
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('Membrane Voltage (mV)', fontsize=14)
plt.legend()
plt.xlim(0, 20)
plt.grid(True)
plt.show()

# Plot gating variables 
plt.figure(figsize=(8,4))
plt.plot(t, n50, label='n (50%)')
plt.plot(t, n25, '--',label='n (25%)')
plt.plot(t, nn, label='n (Normal)')
plt.plot(t, m50, label='m (50%)')
plt.plot(t, m25,'--', label='m (25%)')
plt.plot(t, mn, label='m (Normal)')
plt.plot(t, h50,  label='h (50%)')
plt.plot(t, h25, '--', label='h (25%)')
plt.plot(t, hn,label='h (Normal)')
plt.title('Gating Variables for 25% and 50% GNa Reduction', fontsize=16)
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('Gating Variable', fontsize=14)
plt.legend(loc = 'upper right')
plt.xlim(0, 20)
plt.grid(True)
plt.show()

# Problem 4b

# Function to Create two pulses of the same magnitude a given delta t apart
def IStim_double(t, dt, amp=25):
    if (0.5 <= t <= 2.5) or (0.5 + dt <= t <= 2.5 + dt):
        return amp
    else:
        return 0

# Function to stimulate and measure excitability
def simulate_dt(dt, amp):
    # Use longer time to ensure we capture both pulses
    t = np.linspace(0, 50 + dt, 6000)   
    Is_func = lambda tt: IStim_double(tt, dt, amp)
    
    # Use fresh initial conditions for each simulation
    a_n, b_n, a_m, b_m, a_h, b_h = rates(0)
    n0 = a_n / (a_n + b_n)
    m0 = a_m / (a_m + b_m)
    h0 = a_h / (a_h + b_h)
    V0 = 0
    y0_local = [V0, n0, m0, h0]  # Local initial conditions
    
    sol = odeint(derivatives, y0_local, t, args=(gNa_bar,Is_func))
    V = sol[:, 0]
    
    # Count spikes using simple threshold detection
    spike_inds = np.where((V[1:] > 0) & (V[:-1] <= 0))[0]
    num_spikes = len(spike_inds)
    return num_spikes, V, t

# Reset gNa_bar to reduced value for the first part (50%)
gNa_bar = gNa_Norm * 0.5

dt_values = np.arange(1, 31, 1)  # delta t values 
I1 = 25
I2_ratio_reduced = []

for dt in dt_values:
    amp = I1
    spikes, _, _ = simulate_dt(dt, amp)
    
    if spikes < 2:
        amp_test = I1
        while amp_test < 600:  # upper cap
            spikes, _, _ = simulate_dt(dt, amp_test)
            if spikes >= 2:
                break
            amp_test += 2
        I2_ratio_reduced.append(amp_test / I1)
    else:
        I2_ratio_reduced.append(1.0)

I_ratio_reduced = np.array(I2_ratio_reduced)

# Rerun for 25% GNa
gNa_bar = gNa_Norm * 0.75

I1 = 25
I2_ratio_25 = []

for dt in dt_values:
    amp = I1
    spikes, _, _ = simulate_dt(dt, amp)
    
    if spikes < 2:
        amp_test = I1
        while amp_test < 600:  # upper cap
            spikes, _, _ = simulate_dt(dt, amp_test)
            if spikes >= 2:
                break
            amp_test += 2
        I2_ratio_25.append(amp_test / I1)
    else:
        I2_ratio_25.append(1.0)

I_ratio_25 = np.array(I2_ratio_25)
# Rerun for normal Gna
gNa_bar = gNa_Norm  # Reset to normal value
I2_ratio_norm = []

for dt in dt_values:
    amp = I1
    spikes, _, _ = simulate_dt(dt, amp)

    if spikes < 2:
        amp_test = I1
        while amp_test <600:
            spikes, _, _ = simulate_dt(dt, amp_test)
            if spikes >= 2:
                break
            amp_test += 2
        I2_ratio_norm.append(amp_test / I1)
    else:
        I2_ratio_norm.append(1.0)

I_ratio_norm = np.array(I2_ratio_norm)

# Plot Refractory Period Curve 
plt.figure(figsize=(7,5))
plt.plot(dt_values, I_ratio_norm, label='Normal $g_{Na}$', color='blue', lw=2)
plt.plot(dt_values, I_ratio_reduced, label='50% $g_{Na}$', color='red', lw=2)
plt.plot(dt_values, I_ratio_25, label='25% $g_{Na}$', color='green', lw=2)
plt.axhline(y=1, color='k', linestyle='--', linewidth=1)
plt.xlim(0, 30)
plt.ylim(0, 20)
plt.xlabel(r'$\Delta t$ (ms)', fontsize=14)
plt.ylabel(r'$\frac{I_2}{I_1}$', fontsize=14)
plt.title('Refractory Period Comparison: Normal vs 50% $g_{Na}$ vs 25% $g_{Na}$ ', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True)
plt.show()