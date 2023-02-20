import numpy as np
import matplotlib.pyplot as plt


import scipy.integrate as sci
import scipy.optimize as spo                

def sharks_eat_tuna(x, t):
    # unpack variables, write ODEs
    S, T = x
    dSdt = 0.01*S*T - 0.2*S
    dTdt = 0.05*T - 0.01*S*T
    return np.array([dSdt, dTdt])

# calculate the derivative at any given point
def deriv(point):
    return sharks_eat_tuna(point, 0)

# set grid for phase plane
S, T = np.meshgrid(np.linspace(0, 50, 20), np.linspace(0, 50, 20))


# calculate speed and direction of vectors
S_mag, T_mag = sharks_eat_tuna([S, T], 0)
speed = np.sqrt(S_mag**2 + T_mag**2)
S_mag[:,1:]
# calculate normalized vectors, avoiding division by zero
if speed.any() > .0001:
    S_norm = S_mag/speed
    T_norm = T_mag/speed
else:
    S_norm = 0
    T_norm = 0

# plot the vector field
plt.quiver(S, T, S_norm, T_norm, speed, cmap='viridis')

# plot the vector field, omitting the first column as normalizing magnitudes would divide by zero
plt.quiver(S[:,1:], T[:,1:], S_mag[:,1:]/speed[:,1:], T_mag[:,1:]/speed[:,1:], speed[:,1:], cmap='viridis')

# find the equilibrium points
equilibrium_points = []
for i in range(S.shape[0]):
    for j in range(S.shape[1]):
        point = spo.fsolve(deriv, [S[i,j], T[i,j]])
        if np.all(point >= 0):
            equilibrium_points.append(point)

# plot the equilibrium points
for point in equilibrium_points:
    plt.plot(point[0], point[1], 'ko')

# set axis limits
plt.xlim([-1, 50])
plt.ylim([-1, 50])
plt.xlabel('Sharks')
plt.ylabel('Tuna')
plt.show()
print(equilibrium_points)

#butthole
