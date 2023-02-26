import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import sin

def ImperfectImmunity(state, t):
    # alpha = a * (np.sin(t) + 1)
    I = state[0]
    S = state[1]

    dIdt = np.zeros(t.size)
    dSdt = np.zeros(t.size)

    for I_idx, i in enumerate(I):
            for S_idx, s in enumerate(S):
                dIdt[I_idx] = alpha * i * s - i/tau
                R = N - i - s
                dSdt[S_idx] = -alpha * i * s + R/tau2

    return(dIdt, dSdt)

def Rt(N, It, St):
    Nvec = N * np.ones(St.shape)
    return (Nvec - St - It)

# nullcline function
def get_nullclines():
    # alpha = a * (np.sin(t) + 1)
    I_nullcline = 1 / (alpha * tau)
    S_nullcline = tau * (N - 1 / (alpha * tau)) / (tau2 + tau)
    return (S_nullcline, I_nullcline)

# define parameters
# alpha = 0.4
tau = .1
tau2 = 3
state0 = [1,99]
N = np.sum(state0)
t = np.linspace(0, 10, 1000)
alpha = np.zeros(t.size)
for i, v in enumerate(t):
    alpha[i] = .4 * (sin(v) + 1)


state = odeint(ImperfectImmunity, state0, t)

I_nullcline = get_nullclines()[1]
S_nullcline = get_nullclines()[0]

plt.plot(t, Rt(N, state[:,0], state[:,1]), c='#77dd77', linewidth = 3, label = 'Recovered')
plt.plot(t, state[:, 0], c='#ff6961', linewidth = 3, label = 'Infected')
plt.plot(t, state[:, 1], c='#7b519d', linewidth = 3, label = 'Susceptible')
plt.plot(t, I_nullcline*np.ones(1000), c='#5B2C6F', linestyle = '--', linewidth = 2, label='I nullcline')
plt.plot(t, S_nullcline*np.ones(1000), c='#953732', linestyle = '--', linewidth = 2, label='S nullcline')
plt.legend(fontsize = 10, loc = 'upper right')
plt.xlabel('Time', fontsize =15)
plt.ylabel('Individuals', fontsize =15)
plt.show()





