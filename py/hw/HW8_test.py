import numpy as np
import matplotlib.pyplot as plt

# define parameters
alpha = 0.4
tau_inf = 0.1
tau_imm = 3
N = 100

# nullcline function
def get_nullclines(S, I):

    for i in np.where(tau_inf * alpha * S != 0):
        I_nullcline = 1 / (tau_inf * alpha * S)
        S_nullcline = (N - I) / (tau_imm * alpha)
        return (S_nullcline, I_nullcline)

# define axis ranges
I_range = np.arange(0, 100)
S_range = np.arange(0, 100)

# meshgrid for phase diagram
I_grid, S_grid = np.meshgrid(I_range, S_range)

# calculate derivatives
dI_dt = alpha * I_grid * S_grid - (I_grid / tau_inf)
R = (N - I_grid - S_grid)
dS_dt = -alpha * I_grid * S_grid + R / tau_imm

# plot phase space diagram
fig, ax = plt.subplots(figsize=(8, 8))
ax.quiver(I_grid, S_grid, dI_dt, dS_dt, color='#440154')
ax.plot(I_range, get_nullclines(S_range, I_range)[1], linewidth = 5, label='I-nullcline', color='#5ec962')
ax.plot(I_range, get_nullclines(S_range, I_range)[0], linewidth = 5, label='S-nullcline', color='#3b528b')
ax.set_xlabel('dIdt', fontsize =15)
ax.set_ylabel('dSdt', fontsize =15)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.legend()
plt.show()















import numpy as np
import matplotlib.pyplot as plt

# define parameters
alpha = 0.4
tau_inf = 0.1
tau_imm = 3
N = 100

# nullcline functions
def get_nullclines(S, I):
    I_nullcline = 1 / (alpha * tau_inf)
    S_nullcline = tau_inf * (N - 1 / (alpha * tau_inf)) / (tau_imm + tau_inf)
    return (S_nullcline, I_nullcline)

# define axis ranges
I_range = np.linspace(0, 100, 100)
S_range = np.linspace(0, 100, 100)

# meshgrid for phase diagram
I_grid, S_grid = np.meshgrid(I_range, S_range)

# calculate derivatives
dI_dt = alpha * I_grid * S_grid - (I_grid / tau_inf)
R = (N - I_grid - S_grid)
dS_dt = -alpha * I_grid * S_grid + R / tau_imm

# normalize vector lengths by the hypoteneuse
hypo_len = (np.hypot(dI_dt,dS_dt))
hypo_len[hypo_len==0] = 1.0
dI_dt /= hypo_len
dS_dt /= hypo_len

# plot phase space diagram
fig, ax = plt.subplots(figsize=(8, 8))
ax.quiver(I_grid, S_grid, dI_dt, dS_dt, scale = 100, color='#440154')
ax.plot(I_range, get_nullclines(S_range, I_range)[1]*np.ones(100), linewidth = 3, label='I-nullcline', color='#5ec962')
ax.plot(get_nullclines(S_range, I_range)[0]*np.ones(100), S_range, linewidth = 3, label='S-nullcline', color='#3b528b')
ax.set_xlabel('dIdt', fontsize =15)
ax.set_ylabel('dSdt', fontsize =15)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.legend()
plt.show()






import sympy as sym

# Define symbols
I, S, alpha, tau_inf, tau_imm, N = sym.symbols('I S alpha tau_inf tau_imm N')

# Define the equations
eq1 = alpha * I * S - I / tau_inf
eq2 = -alpha * I * S + (N - I - S) / tau_imm

# Solve for I_dot = 0 and S_dot = 0
sol = sym.solve((eq1, eq2), (I, S))

# Print the solutions
print(sol)