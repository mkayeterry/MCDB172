
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def ImperfectImmunity(state, t):
    I = state[0]
    S = state[1]

    # define system with imperfect immunity
    dIdt = alpha * I * S - I/tau_inf
    R = N - I - S 
    dSdt = -alpha * I * S + R/tau_imm 

    return(dIdt, dSdt)

def Rt(N, It, St):
    Nvec = N * np.ones(St.shape)
    return (Nvec - St - It)

# nullcline functions
def I_nullcline(S):
    return ( 1 / (tau_inf * alpha * S) )

def S_nullcline(I):
    return ( (N - I) / (tau_imm * alpha) )

# define parameters
alpha = 0.4
tau_inf = 0.1
tau_imm = 3

state0 = [1,99]

N = np.sum(state0)
t = np.linspace(0, 10, 1000)

state = odeint(ImperfectImmunity, state0, t)

# plotting a phase diagram
Ilim = [0, 100]
Slim = [0, 100]

npoints = 30 # pick the resolution of our phase portrait

# set up meshgrid function for evaluating the dynamics at many states
h1 = np.linspace(Ilim[0], Ilim[1], npoints)
h2 = np.linspace(Ilim[0], Ilim[1], npoints)

# generate meshgrid
H1, H2 = np.meshgrid(h1, h2)

u, v = np.zeros(H1.shape), np.zeros(H2.shape)

NJ, NK = H1.shape

# two nested for loops to evaluate the dynamics on the meshgrid
for j in range(NJ):
    for k in range(NK):

        Istate = H1[j, k]
        Sstate = H2[j, k]

        statejk = [Istate, Sstate]

        # call differential equation to calc system velocity
        ISdot = ImperfectImmunity(statejk, [])

        u[j,k] = ISdot[0]
        v[j,k] = ISdot[1]


# normalize vector lengths by the hypoteneuse
M = (np.hypot(u,v))
M[M==0] = 1.
u /= M
v /= M

plt.quiver(H1, H2, u, v, color = 'red')
plt.xlabel('Infected Indiviuals', fontsize =15)
plt.ylabel('Susceptible Indiviuals', fontsize =15)
plt.plot(state[:,0], state[:,1], linewidth = 5)
plt.plot(np.linspace(0, 100, 100), I_nullcline(state[0]), linewidth = 5, label = 'I Nullcline')
plt.plot(np.linspace(0, 100, 100), S_nullcline(state[1]), linewidth = 5, label = 'S Nullcline')
# plt.legend(loc = 'upper right')
plt.ylim((0,100))
plt.xlim((0,100))
plt.show()

















import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def ImperfectImmunity(state, t):
    I = state[0]
    S = state[1]

    # define system with imperfect immunity
    dIdt = alpha * I * S - I/tau_inf
    R = N - I - S 
    dSdt = -alpha * I * S + R/tau_imm 

    return(dIdt, dSdt)

def Rt(N, It, St):
    Nvec = N * np.ones(St.shape)
    return (Nvec - St - It)

# nullcline function
def get_nullclines(S, I):
    I_nullcline = 1 / (tau_inf * alpha * S)
    S_nullcline = (N - I) / (tau_imm * alpha)
    return (S_nullcline, I_nullcline)

# define parameters
alpha = 0.4
tau_inf = 0.1
tau_imm = 3
state0 = [1,99]
N = np.sum(state0)
t = np.linspace(0, 10, 1000)

state = odeint(ImperfectImmunity, state0, t)

# plotting a phase diagram
I_lim = [0, 100]
S_lim = [0, 100]

npoints = 30 # pick the resolution of our phase portrait

# set up meshgrid function for evaluating the dynamics at many states
I_range = np.linspace(I_lim[0], I_lim[1], npoints)
S_range = np.linspace(I_lim[0], I_lim[1], npoints)

# generate meshgrid
I_grid, S_grid = np.meshgrid(I_range, S_range)

I_vel, S_vel = np.zeros(I_grid.shape), np.zeros(S_grid.shape)

I_points, S_points = I_grid.shape


# two nested for loops to evaluate the dynamics on the meshgrid
for i in range(I_points):
    for j in range(S_points):

        I_state = I_grid[i, j]
        S_state = S_grid[i, j]

        current_state = [I_state, S_state]

        # call differential equation to calc system velocity
        ISdot = ImperfectImmunity(current_state, [])

        I_vel[i,j] = ISdot[0]
        S_vel[i,j] = ISdot[1]

# normalize vector lengths by the hypoteneuse
hypo_len = (np.hypot(I_vel,S_vel))
hypo_len[hypo_len==0] = 1.0
I_vel /= hypo_len
S_vel /= hypo_len

# initialize arrays to store the nullclines
S_nullclines = np.zeros_like(I_grid)
I_nullclines = np.zeros_like(S_grid)

# two nested for loops to evaluate the dynamics on the meshgrid
for i in range(I_points):
    for j in range(S_points):

        I_state = I_grid[i, j]
        S_state = S_grid[i, j]

        # evaluate the nullclines at the current state
        S_nullcline, I_nullcline = get_nullclines(S_state, I_state)

        # store the nullclines in arrays
        S_nullclines[i, j] = S_nullcline
        I_nullclines[i, j] = I_nullcline

plt.quiver(I_grid, S_grid, I_vel, S_vel, color = 'red')
plt.xlabel('Infected Indiviuals', fontsize =15)
plt.ylabel('Susceptible Indiviuals', fontsize =15)
plt.plot(state[:,0], state[:,1], linewidth = 5)
plt.plot(np.linspace(0, 100, 100), I_nullclines, linewidth = 5, label = 'I Nullcline')
# plt.plot(np.linspace(0, 100, 100), S_nullclines, linewidth = 5, label = 'S Nullcline')
# plt.legend(loc = 'upper right')
# plt.ylim((0,100))
# plt.xlim((0,100))
plt.show()