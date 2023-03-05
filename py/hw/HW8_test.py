import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint



def fitzHughNagumo_ForceFxn(state, t, Ton, q):
    
    v = state[0]
    w = state[1]

    dv = -v * (v - a) * (v - 1) - w + force(t, Ton, q)
    dw = eps * (v - (gamma * w))

    return(dv, dw)

# forcing function
def force(t, Ton, q):

    if t < Ton:
        I = 0

    if t >= Ton:
        I = q
        
    return I

# for plotting
def computeI(t, Ton, q):

    I = np.zeros(len(t))

    for idx in np.arange(0, len(t), 1):
        I[idx] = force(t[idx], Ton, q)

    return I

# minimum input amplitude for full spiking
def min_full_spike(state0, t, Ton):

    q_lst = np.linspace(0, 0.2, 100)
    max_v = []

    for q in q_lst:
        state = odeint(fitzHughNagumo_ForceFxn, state0, t, args=(Ton, q))
        max_v.append(np.max(state[:, 0]))

    idx = np.argmax(np.array(max_v) >= 0.9)
    q_thresh = q_lst[idx]

    return q_thresh


# minimum input amplitude for repetitive firing
def min_rep_fire(state0, t, Ton, max_iter=1000):
    n_spikes = 0
    q_thresh = min_full_spike(state0, t, Ton)
    q_rep_thresh = q_thresh
    
    state = odeint(fitzHughNagumo_ForceFxn, state0, t, args=(Ton, q_rep_thresh))
    
    for i in range(1, t.size):
        if state[:, 0][i] > q_thresh and state[:, 0][i-1] < q_thresh:
            n_spikes += 1
    
    while n_spikes < 3:
        q_rep_thresh += 0.001
        state = odeint(fitzHughNagumo_ForceFxn, state0, t, args=(Ton, q_rep_thresh))
        n_spikes = 0
        for i in range(1, len(t)):
            if state[:, 0][i] > q_thresh and state[:, 0][i-1] < q_thresh:
                n_spikes += 1
        max_iter -= 1
    
    return q_rep_thresh



# define parameters
state0 = [0, 0]
t = np.linspace(0, 300, 1500)

eps = 0.01
a = 0.1
gamma = 5
Ton = 50

state = odeint(fitzHughNagumo_ForceFxn, state0, t, args=(Ton, min_full_spike(state0, t, Ton)))

print(f"Minimum input amplitude for full spiking: {min_full_spike(state0, t, Ton):.3f}")
print(f"Minimum input amplitude for repetitive firing: {min_rep_fire(state0, t, Ton):.3f}")

plt.plot(t, state[:, 0], label = 'Voltage')
plt.plot(t, computeI(t, Ton, min_full_spike(state0, t, Ton)), label = 'Conductance')
plt.legend()
plt.show()

