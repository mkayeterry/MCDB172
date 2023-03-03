import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

# E = golden eagle (predator)
E0 = 28 # placeholder

# F = fox (prey)
F0 = 1312

# S = skunk (prey)
S0 = 1000

# P = piglet (prey)
P0 = 13827

# r_i = instrinsic growth rate
r_f = 0.32
r_p = 0.78
r_s = r_f

# K_i = carrying capacity
K_f = 1544
K_s = 2490
K_p = 15189

# beta_ij = measure of energetic resource competition between foxes and skunks only
beta_sf = 500.48 / 181.58
beta_fs = 181.58 / 500.48

# mu_i = predation rate by eagles
mu_f = 0.086
mu_p = 0.019
mu_s = 0.159

# phi = term of eagle preference for foxes
phi = 8.1

# sigma = term of eagle preference for skunks relative to piglets
sigma = 3.1

# nu = eagle mortality rate
nu = 0.09

# lambda_i = the rate at which prey i are turned into new predators
lambda_f = 7.7e-4
lambda_p = 7.7e-4
lambda_s = 2.5e-4

# ------------------------------------------------------------------------------------------------
# Without Pigs
def LV(state, t):
    F, S, E = state

    dFdt = r_f * F * (1 - ((F + (beta_fs * S)) / K_f)) - mu_f * ((phi * F) / ((phi * F) + (sigma * S))) * E * F
    dSdt = r_s * S * (1 - ((S + (beta_sf * F)) / K_s)) - mu_s * ((sigma * S) / ((phi * F) + (sigma * S))) * E * S
    dEdt = (((lambda_f * mu_f * phi * (F**2)) + (lambda_s * mu_s * sigma * (S**2)) * E) / ((phi * F) + (sigma * S))) - (nu * E)
    return np.array([dFdt, dSdt, dEdt])

t = np.linspace(0, 100, 15)

state0 = [F0, S0, E0]

state = odeint(LV, state0, t)

F, S, E = state.T

def y1_y2(x):
    return ((2*x)/100)

def y2_y1(x):
    return x*50

fig, ax = plt.subplots()
ax.plot(t, F, '-', linewidth = 3, label = 'Foxes', color = 'blue')
ax.plot(t, S, '-', linewidth = 3, label = 'Skunks', color = 'lightgreen')
ax.plot(t, E, '-', linewidth = 3, label = 'Eagles', color = '#ee5727')
ax.set_xlabel('Time')
ax.set_ylabel('Foxes, Skunks')
secax = ax.secondary_yaxis('right', functions = (y1_y2, y2_y1))
secax.set_ylabel('Eagles')
secax.set_ylim(0, 30)
ax.set_xlim(0, 100)
ax.set_ylim(0, 1400)
ax.legend()
plt.show()
# -------------------------------------------------------------------------------

# With Pigs
def LV_pigs(state_2, t):
    F_2, S_2, P_2, E_2 = state_2

    dFdt = (r_f * F_2 * (1 - ((F_2 + (beta_fs * S_2)) / K_f))) - mu_f * ((phi * F_2) / ((phi * F_2) + (sigma * S_2) + (P_2))) * E_2 * F_2
    dSdt = (r_s * S_2 * (1 - ((S_2 + (beta_sf * F_2)) / K_s))) - mu_s * ((sigma * S_2) / ((phi * F_2) + (sigma * S_2) + (P_2))) * E_2 * S_2
    dPdt = (r_p * P_2 * (1 - (P_2 / K_p))) - mu_p * (P_2 / ((phi * F_2) + (sigma * S_2) + (P_2))) * E_2 * P_2
    dEdt = (((lambda_f * mu_f * phi * (F_2**2)) + (lambda_s * mu_s * sigma * (S_2**2)) + (lambda_p * mu_p * (P_2**2)) * E_2) / ((phi * F_2) + (sigma * S_2) + (P_2))) - (nu * E_2)
    return np.array([dFdt, dSdt, dPdt, dEdt])

state0_2 = [F0, S0, P0, E0]

state_2 = odeint(LV_pigs, state0_2, t)

F_2, S_2, P_2, E_2 = state_2.T

def y1_y2(x):
    return ((2*x)/100)

def y2_y1(x):
    return x*50

fig, ax = plt.subplots()
ax.plot(t, F_2, '-', linewidth = 3, label = 'Foxes', color = 'blue')
ax.plot(t, S_2, '-', linewidth = 3, label = 'Skunks', color = 'lightgreen')
ax.plot(t, E_2, '-', linewidth = 3, label = 'Eagles', color = '#ee5727')
ax.plot(t, P_2/1000, '-', linewidth = 3, label = 'Pigs')
ax.set_xlabel('Time')
ax.set_ylabel('Foxes, Skunks')
secax = ax.secondary_yaxis('right', functions = (y1_y2, y2_y1))
secax.set_ylabel('Eagles')
secax.set_ylim(0, 30)
ax.set_xlim(0, 100)
ax.set_ylim(0, 1400)
ax.legend()
plt.show()


