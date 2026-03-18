import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Constants
# Problem 3a
gNa_bar = 70.7
gK_bar  = 24.31
ENa = 110
EK = -12
Cm = 1
gL_bar = 0.3
EL = 0
t = np.linspace(0, 20, 2000)
# Function to solve for rates
def rates(Vm):
    if abs(Vm + 10) < 1e-6:
        a_n = 0.1
    else:
        a_n = (10 - Vm) / (100 * (np.exp((10 - Vm) / 10) - 1))
    b_n = 0.125 * np.exp(-Vm / 80)

    if abs(Vm + 25) < 1e-6:
        a_m = 1
    else:
        a_m = (25 - Vm) / (10 * (np.exp((25 - Vm) / 10) - 1))
    b_m = 4 * np.exp(-Vm / 18)

    a_h = 0.07 * np.exp(-Vm / 20)
    b_h = 1 / (np.exp((30 - Vm) / 10) + 1)

    return a_n, b_n, a_m, b_m, a_h, b_h
# Function to Output the stimulus current
def IStim(t):
    return 20 if 0.5 <= t <= 2.5 else 0
# Function to Solve the HH Equations
def derivatives(y, t, I_func):
    V, n, m, h = y
    a_n, b_n, a_m, b_m, a_h, b_h = rates(V)

    dndt = a_n*(1-n) - b_n*n
    dmdt = a_m*(1-m) - b_m*m
    dhdt = a_h*(1-h) - b_h*h

    IK  = gK_bar * n**4 * (V - EK)
    INa = gNa_bar * m**3 * h * (V - ENa)
    IL  = gL_bar * (V - EL)

    dVdt = (I_func(t) - (IK+INa+IL)) / Cm
    return [dVdt, dndt, dmdt, dhdt]

# Initial gating variables
a_n, b_n, a_m, b_m, a_h, b_h = rates(0)
# Compute initial gating variables at V = 0
a_n, b_n, a_m, b_m, a_h, b_h = rates(0)
n0 = a_n / (a_n + b_n)
m0 = a_m / (a_m + b_m)
h0 = a_h / (a_h + b_h)
y0 = [0, a_n/(a_n+b_n), a_m/(a_m+b_m), a_h/(a_h+b_h)]

sol = odeint(derivatives, y0, t, args=(IStim,))
V, n, m, h = sol.T

# Compute currents 
IK  = gK_bar * n**4 * (V - EK)
INa = gNa_bar * m**3 * h * (V - ENa)
IL  = gL_bar * (V - EL)

# Plot Stimulus Current
I_vals = np.array([IStim(T) for T in t])  # Evaluate IStim over time vector

plt.figure(figsize=(7,4))
plt.step(t, I_vals, where='post', linewidth=2)
plt.title("Applied Stimulus Current",fontsize = 16)
plt.xlabel("Time (ms)",fontsize = 14)
plt.ylabel("I_stim (nA)",fontsize = 14)
plt.ylim(-1, max(I_vals) + 5)  
plt.xlim(0,20)
plt.grid(True)
plt.show()

# Plot Ionic Currents  
plt.figure(figsize=(7,4))
plt.plot(t, IK, label="I_K")
plt.plot(t, INa, label="I_Na")
plt.plot(t, IL, label="I_L")
plt.title("Simulated Hodgkin-Huxley Ionic Currents",fontsize = 16)
plt.xlabel("Time (ms)",fontsize = 14)
plt.ylabel("Ionic Current (nA)",fontsize = 14)
plt.xlim(0,20)
plt.legend()
plt.grid(True)
plt.show()
# Plot Action Potential
plt.figure(figsize=(8,5))
plt.plot(t, V, lw=2, color='purple')
plt.title("Simulated Hodgkin–Huxley Action Potential",fontsize = 16)
plt.xlabel("Time (ms)",fontsize = 14)
plt.ylabel("Membrane Voltage (mV)",fontsize = 14)
plt.grid(True)
plt.xlim(0,20)
plt.show()
# Plot Gating Variables
plt.figure(figsize=(7,4))
plt.plot(t, n, label="n")
plt.plot(t, m, label="m")
plt.plot(t, h, label="h")
plt.title("Simulated Hodgkin-Huxley Gating Variables",fontsize = 16)
plt.xlabel("Time (ms)",fontsize = 14)
plt.ylabel("Gating Value",fontsize = 14)
plt.xlim(0,20)
plt.legend()
plt.grid(True)
plt.show()

 # Problem 3b : Forward Euler 
dt = 0.0115
y = np.array([0, n0, m0, h0], dtype=float)

V_FE = np.zeros_like(t)
n_FE = np.zeros_like(t)
m_FE = np.zeros_like(t)
h_FE = np.zeros_like(t)

IK_FE = np.zeros_like(t)
INa_FE = np.zeros_like(t)
IL_FE = np.zeros_like(t)

for i in range(len(t)):
    V_FE[i], n_FE[i], m_FE[i], h_FE[i] = y

    # Currents
    IK_FE[i]  = gK_bar * n_FE[i]**4 * (V_FE[i] - EK)
    INa_FE[i] = gNa_bar * m_FE[i]**3 * h_FE[i] * (V_FE[i] - ENa)
    IL_FE[i]  = gL_bar * (V_FE[i] - EL)

    # Forward Euler update
    y = y + dt * np.array(derivatives(y, t[i], IStim))


# Plot AP
plt.figure(figsize=(8,5))
plt.plot(t, V_FE, color='purple')
plt.title("Simulated HH AP using the Forward Euler Approximation",fontsize = 16)
plt.xlabel("Time (ms)",fontsize = 14)
plt.ylabel("Membrane Voltage (mV)",fontsize = 14)
plt.xlim(0,20)
plt.grid(True)
plt.show()

# Plot Ionic Currents
plt.figure(figsize=(7,4))
plt.plot(t, IK_FE, label="IK")
plt.plot(t, INa_FE, label="INa")
plt.plot(t, IL_FE, label="IL")
plt.title("Simulated HH Ionic Currents using the Forward Euler Approximation",fontsize = 16)
plt.xlabel("Time (ms)",fontsize = 14)
plt.ylabel("Current (nA)",fontsize = 14)
plt.legend()
plt.xlim(0,20)
plt.grid(True)
plt.show()

# Plot Gates
plt.figure(figsize=(7,4))
plt.plot(t,n_FE, label = 'n')
plt.plot(t,m_FE, label = 'm')
plt.plot(t,h_FE, label = 'h')
plt.title('Simulated HH Gating Variables using the Forward Euler Approximation',fontsize = 16)
plt.xlabel('Time(ms)',fontsize = 14)
plt.ylabel('Gating Variable',fontsize = 14)
plt.legend()
plt.xlim(0,20)
plt.grid(True)
plt.show()

#  Problem 3c: Midpoint Method 
y_MP = np.array([0, n0, m0, h0], dtype=float) 
dtm = 0.05
t = np.arange(0, 20, dtm)

V_MP = np.zeros_like(t)
n_MP = np.zeros_like(t)
m_MP = np.zeros_like(t)
h_MP = np.zeros_like(t)

IK_MP = np.zeros_like(t)
INa_MP = np.zeros_like(t)
IL_MP = np.zeros_like(t)

for i in range(len(t)):
    V_MP[i], n_MP[i], m_MP[i], h_MP[i] = y_MP

    IK_MP[i]  = gK_bar * n_MP[i]**4 * (V_MP[i] - EK)
    INa_MP[i] = gNa_bar * m_MP[i]**3 * h_MP[i] * (V_MP[i] - ENa)
    IL_MP[i]  = gL_bar * (V_MP[i] - EL)

    k1 = np.array(derivatives(y_MP, t[i], IStim))
    k2 = np.array(derivatives(y_MP + dtm/2 * k1, t[i] + dtm/2, IStim))
    y_MP = y_MP + dtm * k2


# Plot AP
plt.figure(figsize=(8,5))
plt.plot(t, V_MP, color='purple')
plt.title("Simulated HH AP using the Midpoint Method",fontsize = 16)
plt.xlabel("Time (ms)",fontsize = 14)
plt.ylabel("Vm (mV)",fontsize = 14)
plt.grid(True)
plt.xlim(0,20)
plt.show()

# Plot Ionic Currents
plt.figure(figsize=(7,4))
plt.plot(t, IK_MP, label="IK")
plt.plot(t, INa_MP, label="INa")
plt.plot(t, IL_MP, label="IL")
plt.title("Simulated HH Ionic Currents using the Midpoint Method",fontsize = 16)
plt.xlabel("Time (ms)",fontsize = 14)
plt.ylabel("Current (nA)",fontsize = 14)
plt.legend()
plt.grid(True)
plt.xlim(0,20)
plt.show()

# Plot Gates
plt.figure(figsize=(7,4))
plt.plot(t,n_MP, label = 'n')
plt.plot(t,m_MP, label = 'm')
plt.plot(t,h_MP, label = 'h')
plt.title('Simulated HH Gating Variables using the Midpoint Method',fontsize = 16)
plt.xlabel('Time(ms)',fontsize = 14)
plt.ylabel('Gating Variable',fontsize = 14)
plt.legend()
plt.grid(True)
plt.xlim(0,20)
plt.show()












































