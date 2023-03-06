import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

gammaV = 510
gammaVA = 619.2
gammaVH = 1.02
alphaV = 1.7
aV1 = 100
aV2 = 23000
bHD = 4
aR = 1
gammaHV = 0.34
bHF = 0.01
bIE = 0.066
aI = 1.5
bMD = 1
bMV = 0.0037
aM = 1
bF = 250000
cF = 2000
bFH = 17
aF = 8
bEM = 8.8
bEI = 2.72
aE = 0.4
bPM = 11.5
aP = 0.4
bA = 0.043
gammaAV = 146.2
aA = 0.043
r = 3e-5

# standard behavior
V0 = 0.01
H0 = 1
I0 = 0
M0 = 0
F0 = 0
R0 = 0
E0 = 1
P0 = 1
A0 = 1
S0 = 0.1


def deriv_std(state, t):
    V, H, I, M, F, R, E, P, A, S = state

    dVdt = (gammaV * I) - (gammaVA * S * A * V) - (gammaVH * H * V) - (alphaV * V) - ((aV1 * V)/(1 + (aV2 * V)))
    dHdt = ((bHD * D) * (H + R)) + (R * R) - (gammaHV * V * H) - (bHF * F * H)
    dIdt = (gammaHV * V * H) - (bIE * E * I) - (aI * I)
    dMdt = (((bMD * D) + (bMV * V)) * (1 - M)) - (aM * M)
    dFdt = (bF * M) + (cF * I) - (bFH * H * F) - (aF * F)
    dRdt = (bHF * F * H) - (aR * R)
    dEdt = (bEM * M * E) - (bEI * I * E) + (aE * (1 - E))
    dPdt = (bPM * M * P) + (aP * (1 - P))
    dAdt = (bA * P) - (gammaAV * S * A * V) - (aA * A)
    dSdt = (r * P) * (1 - S)

    return np.array([dVdt, dHdt, dIdt, dMdt, dFdt, dRdt, dEdt, dPdt, dAdt, dSdt])

t = np.linspace(0, 15, 100)

state0 = [V0, H0, I0, M0, F0, R0, E0, P0, A0, S0]

state = odeint(deriv_std, state0, t)

V, H, I, M, F, R, E, P, A, S = state.T
D = 1 - H - R - I


# Create a figure with 5 rows and 2 columns of subplots
fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(5, 10))
fig.set_figwidth(17)

# Plot each variable on a separate subplot
axes[0, 0].plot(t, V)
axes[0, 0].set_xlabel('Days', fontsize = 13)
axes[0, 0].set_ylabel('V', fontsize = 13)

axes[0, 1].plot(t, H)
axes[0, 1].set_xlabel('Days', fontsize = 13)
axes[0, 1].set_ylabel('H', fontsize = 13)

axes[0, 2].plot(t, I)
axes[0, 2].set_xlabel('Days', fontsize = 13)
axes[0, 2].set_ylabel('I', fontsize = 13)

axes[1, 0].plot(t, M)
axes[1, 0].set_xlabel('Days', fontsize = 13)
axes[1, 0].set_ylabel('M', fontsize = 13)

axes[1, 1].plot(t, F)
axes[1, 1].set_xlabel('Days', fontsize = 13)
axes[1, 1].set_ylabel('F', fontsize = 13)

axes[1, 2].plot(t, R)
axes[1, 2].set_xlabel('Days', fontsize = 13)
axes[1, 2].set_ylabel('R', fontsize = 13)

axes[2, 0].plot(t, E)
axes[2, 0].set_xlabel('Days', fontsize = 13)
axes[2, 0].set_ylabel('E', fontsize = 13)

axes[2, 1].plot(t, P)
axes[2, 1].set_xlabel('Days', fontsize = 13)
axes[2, 1].set_ylabel('P', fontsize = 13)

axes[2, 2].plot(t, A)
axes[2, 2].set_xlabel('Days', fontsize = 13)
axes[2, 2].set_ylabel('A', fontsize = 13)

axes[3, 0].plot(t, D)
axes[3, 0].set_xlabel('Days', fontsize = 13)
axes[3, 0].set_ylabel('D', fontsize = 13)

axes[3, 1].plot(t, S)
axes[3, 1].set_xlabel('Days', fontsize = 13)
axes[3, 1].set_ylabel('S', fontsize = 13)



t_pad = np.pad(t, (5, 0), 'constant')
H_pad = np.pad(H, (5, 0), 'constant')
R_pad = np.pad(R, (5, 0), 'constant')
I_pad = np.pad(I, (5, 0), 'constant')

axes[3, 2].plot(t_pad, H_pad, color = 'blue')
axes[3, 2].plot(t_pad, R_pad, color = 'green')
axes[3, 2].plot(t_pad, I_pad*500, color = 'red')
axes[3, 2].plot(t_pad, np.zeros(105), color = 'black')
axes[3, 2].set_xlabel('Days', fontsize = 13)
axes[3, 2].set_ylabel('Cell Proportions', fontsize = 13)

# Add some space between the subplots
fig.tight_layout()

# Show the figure
plt.show()