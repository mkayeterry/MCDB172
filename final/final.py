## Variables:

# E = golden eagle (predator)
E = 1 # placeholder

# F = fox (prey)
F = 1 # placeholder
F0 = 1312

# S = skunk (prey)
S = 1 # placeholder
S0 = 1000

# P = piglet (prey)
P = 1 # placeholder
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

## Differential Equations
dFdt = ( r_f * F * (1 - ((F + (beta_fs * S)) / K_f)) ) - ( mu_f * ((phi * F) / ((phi * F) + (sigma * S) + P)) * E * F )
dSdt = ( r_s * S * (1 - ((S + (beta_sf * F)) / K_s)) ) - ( mu_s * ((sigma * S) / ((phi * F) + (sigma * S) + P)) * E * S )
dPdt = ( r_p * P * (1 - (P / K_p)) ) - ( mu_p * (P / ((phi * F) + (sigma * S) + P)) * E * P )
dEdt = (((lambda_f * mu_f * phi * (F**2)) + (lambda_s * mu_s * sigma * (S**2)) + (lambda_p * mu_p * (P**2)) * E) / ((phi * F) + (sigma * S) + P)) - (nu * E)
