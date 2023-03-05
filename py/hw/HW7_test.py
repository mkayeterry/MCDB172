import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import sin

def ImperfectImmunity_latent(state, t):
    I = state[0]
    X = state[1]
    S = state[2]


    # define system with imperfect immunity and a latent period
    dIdt = a * I * S - I/tau_inf
    dXdt =  a * (I + X) * S - X/tau_lat
    R = N - I - S 
    dSdt = -a * I * S + R/tau_imm

    return(dIdt, dXdt, dSdt)

# define parameters
a = .4
tau_inf = .1
tau_imm = 3
tau_lat = tau_inf - tau_imm
state0 = [1, 0, 99]
N = np.sum(state0)
t = np.linspace(0, 20, 1000)

state = odeint(ImperfectImmunity_latent, state0, t)

plt.plot(t, state[:, 0], c='#ff6961', linewidth=3, label='Infected')
plt.plot(t, state[:, 1], c='#FFA500', linewidth=3, label='Exposed')
plt.plot(t, state[:, 2], c='#7b519d', linewidth=3, label='Susceptible')
plt.legend(fontsize=10, loc='upper right')
plt.xlabel('Time', fontsize=15)
plt.ylabel('Individuals', fontsize=15)
plt.xlim(-.2, 10)
plt.ylim(-1, 100)
plt.show()

