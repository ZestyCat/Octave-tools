# SAE ARP1533 Matrix Solver

# Fuel composition
m = 7.1576 # C
n = 13.9187 # H 
p = 0.00004 # O
q = 0 # N 
r = 0.0012 # S 

# Air composition
h = 0.01261 # Water
R = 0.20687 # Oxygen
S = 0.78036 # Nitrogen
T = 0.00032 # Carbon dioxide
U = 0.0000037 # Methane

# Measured concentrations (mole fraction)
CO2  = 1.77
CO   = 1
CxHy = 1
NOx  = 1
NO   = 1
x    = 1 # Carbon
y    = 1 # Hydrogen
hsd  = 0 # Semi-dry water concentration

# Correction factors
L     = -1.3e-4 # Mole CO per mole CO2 (zero shift)
M     = -4.5e-4 # Mole CO per mole H2O (zero shift)
L_prm = 0.14 # Percent reading of NO per percent of CO2 (concentration factor effect)
M_prm = 0.28 # Percent reading of NO per percent of H2O (concentration factor effect)
J     = 0.09 # Percent reading of CO2 per percent of O2 (concentration factor effect)
eta   = 0.95 # NOx converter efficiency

A = [0,                         1,                          0, 0,                          0,                         1,  x,   0,  0,  0, -T - U;
     0,                         0,                          0, 0,                          2,                         0,  y,   0,  0,  0, -2 * h - 4 * U;
     0,                         2,                          0, 2,                          1,                         1,  0,   2,  1,  2, -2 * R - 2 * T - h;
     0,                         0,                          2, 0,                          0,                         0,  0,   1,  1,  0, -2 * S;
     0,                         0,                          0, 0,                          0,                         0,  0,   0,  0,  1, 0;
     CO2/(1-hsd),               -1,                         0, J*CO2,                      -CO2/(1-hsd),              0,  0,   0,  0,  0, 0;
     (CO+hsd*M)/(1-hsd),        L,                          0, 0,                          -(CO+hsd*M)/(1-hsd),       -1, 0,   0,  0,  0, 0;
     CxHy,                      0,                          0, 0,                          0,                         0,  -x,  0,  0,  0, 0;
     (1+hsd*M_prm)*NOx/(1-hsd), L_prm*NOx,                  0, -(1+hsd*M_prm)*NOx/(1-hsd), 0,                         0,  eta, -1, -1, 0, 0;
     (1+hsd*M_prm)*NO/(1-hsd),  L_prm*NO,                   0, 0,                          -(1+hsd*M_prm)*NO/(1-hsd), 0,  0,   0,  -1, 0, 0;
     -1,                        1,                          1, 1,                          1,                         1,  1,   1,  1,  1, 0]

B = [m; n; p; q; r; 0; 0; 0; 0; 0; 0]

C = inv(A) * B
