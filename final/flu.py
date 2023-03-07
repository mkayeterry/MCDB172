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
    D = 1 - H - R - I

    dVdt = (gammaV * I) - (gammaVA * S * A * V) - (gammaVH * H * V) - (alphaV * V) - ((aV1 * V)/(1 + (aV2 * V)))
    dHdt = ((bHD * D) * (H + R)) + (aR * R) - (gammaHV * V * H) - (bHF * F * H)
    dIdt = (gammaHV * V * H) - (bIE * E * I) - (aI * I)
    dMdt = (((bMD * D) + (bMV * V)) * (1 - M)) - (aM * M)
    dFdt = (bF * M) + (cF * I) - (bFH * H * F) - (aF * F)
    dRdt = (bHF * F * H) - (aR * R)
    dEdt = (bEM * M * E) - (bEI * I * E) + (aE * (1 - E))
    dPdt = (bPM * M * P) + (aP * (1 - P))
    dAdt = (bA * P) - (gammaAV * S * A * V) - (aA * A)
    dSdt = (r * P) * (1 - S)

    return np.array([dVdt, dHdt, dIdt, dMdt, dFdt, dRdt, dEdt, dPdt, dAdt, dSdt])


t = np.linspace(0, 15, 1000)

state0 = [V0, H0, I0, M0, F0, R0, E0, P0, A0, S0]

state = odeint(deriv_std, state0, t)

V, H, I, M, F, R, E, P, A, S = state.T

D = 1 - H - R - I

# Figure 3 a and b
V0_lst = [0.001, 0.00227, 0.01, 0.1, 200]
V_lst = []
D_lst = []

for i in V0_lst:

    state0 = [i, H0, I0, M0, F0, R0, E0, P0, A0, S0]
    state = odeint(deriv_std, state0, t)

    V_lst.append(state.T[0])

    # need help here
    Dt = D * (V0 ** i)
    D_lst.append(Dt)


plt.plot(t, V_lst[0], label = 'V0 = 0.001')
plt.plot(t, V_lst[1], label = 'V0 = 0.00227')
plt.plot(t, V_lst[2], label = 'V0 = 0.01')
plt.plot(t, V_lst[3], label = 'V0 = 0.1')
plt.plot(t, V_lst[4], label = 'V0 = 200')
plt.xlabel('Days', fontsize = 15)
plt.ylabel('V', fontsize = 15)
plt.legend()
plt.show()

plt.plot(t, D_lst[2], label = 'V0 = 0.001')
plt.plot(t, D_lst[3], label = 'V0 = 0.00227')
plt.plot(t, D_lst[4], label = 'V0 = 0.01')
plt.xlabel('Days', fontsize = 15)
plt.ylabel('D', fontsize = 15)
plt.legend()
plt.show()