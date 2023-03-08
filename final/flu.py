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


# Figure 3 b
V0_lst = [0.001, 0.00227, 0.01, 0.1, 200]
V_lst = []
D_lst = []

for i in range(len(V0_lst)):
    state0 = [V0_lst[i], H0, I0, M0, F0, R0, E0, P0, A0, S0]
    state = odeint(deriv_std, state0, t)

    V = state[:, 0]
    H = state[:, 1]
    I = state[:, 2]
    R = state[:, 5]

    D = 1 - H - R - I
    D_V = D * V

    V_lst.append(V)
    D_lst.append(D_V)


plt.figure(figsize=(8,6))
for i in range(len(V0_lst)):
    plt.plot(t, V_lst[i], label='V0 = {}'.format(V0_lst[i]))
plt.xlabel('Days', fontsize=15)
plt.ylabel('V', fontsize=15)
plt.legend()
plt.show()


plt.figure(figsize=(8,6))
for i in range(2, len(V0_lst)):
    plt.plot(t, D_lst[i], label='V0 = {}'.format(V0_lst[i]))
plt.xlabel('Days', fontsize=15)
plt.ylabel('D', fontsize=15)
plt.legend()
plt.show()

# Phase diagram
V_max = 1.3e10
D_max = 36
V1 = V_lst[1]
V2 = V_lst[2]
V3 = V_lst[3]
D1 = D_lst[1]
D2 = D_lst[2]
D3 = D_lst[3]


xlim = [10e-7, 10e2] # log D range
ylim = [0, 10] # log V range

npoints = 20

s1 = np.logspace(xlim[0], xlim[1], npoints)
s2 = np.logspace(ylim[0], ylim[1], npoints)

S1, S2 = np.meshgrid(s1, s2)

q1, p1 = np.zeros(S1.shape), np.zeros(S2.shape)
# q2, p2 = np.zeros(S1.shape), np.zeros(S2.shape)
# q3, p3 = np.zeros(S1.shape), np.zeros(S2.shape)

NI, NJ = S1.shape

# plotting phase diagrams
for i in range(NI):
    for j in range(NJ):
        V = S2[i, j]
        D = S1[i, j]

???

# plot quivers for each section of the phase diagram
fig, ax = plt.subplots(figsize=(8, 8))

# extreme disease section
ax.quiver(np.log10(S1[:, :]), np.log10(S2[:, :]), p1[:, :], q1[:, :], color='red', alpha=0.5, scale=5)
# # typical disease section
# ax.quiver(np.log10(S1[:, :]), np.log10(S2[:, :]), p2[:, :], q2[:, :], color='blue', alpha=0.5, scale=5)
# # healthy section
# ax.quiver(np.log10(S1[:, :]), np.log10(S2[:, :]), p3[:, :], q3[:, :], color='green', alpha=0.5, scale=5)

# plot trajectory lines
plt.plot(np.log10(D_lst[1]), np.log10(V_lst[1]), color='black')
# plt.plot(np.log10(D_lst[2]), np.log10(V_lst[2]), color='black')
# plt.plot(np.log10(D_lst[3]), np.log10(V_lst[3]), color='black')

# add labels and legend
plt.xlabel('log(D)', fontsize=15)
plt.ylabel('log(V)', fontsize=15)
plt.legend(['V0 = 0.00227', 'V0 = 0.01', 'V0 = 0.1', 'Thresholds: V1=0.00227, V2=0.1'], fontsize=12)

plt.show()

