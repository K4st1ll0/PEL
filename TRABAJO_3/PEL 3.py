from numpy import sqrt

# DATOS GENERALES
Mcp = 1500
Mf_N2H4 = 118

rho_N2H4 = 1032

E1 = 8      # N
E2 = 0.7    # N

x = 0.4
Pi_tobera = 80     # Pc/Ps

Cp_H2 = 29.9
Cp_N2 = 31.1
Cp_NH3 = 54
Cp_N2H4 = 98.8

deltaH_N2H4 = 50.6
deltaH_NH3 = -45.9

Tref = 298.15
Ru = 8314.46

Pc_bar = 10
Pc = Pc_bar * 1e5
Ps = Pc / Pi_tobera

### APARTADO A)
M_N = 14
M_H = 1

B = (4/3) * (1 - x)       # NH3
C = (1/3) * (1 + 2*x)     # N2
D = 2 * x                 # H2

n_totales = B + C + D

# MASAS MOLARES [g/mol]
M_NH3 = M_N + 3*M_H
M_N2  = 2*M_N
M_H2  = 2*M_H

M_prod = (B*M_NH3 + C*M_N2 + D*M_H2) / n_totales   # g/mol
M_prod_kg = M_prod / 1000                           # kg/mol

Cp_prod = (B*Cp_NH3 + C*Cp_N2 + D*Cp_H2) / n_totales   # J/mol/K

Cp_prod_m = Cp_prod / M_prod_kg     # J/kg/K

Rp = Ru / (M_prod)    

gamma = Cp_prod_m / (Cp_prod_m - Rp)
GAMMA = sqrt(gamma) * (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))

Hf_N2H4 = 50.6e3
Hf_NH3 = -45.9e3

Tc = Tref + (Hf_N2H4 - B*Hf_NH3)/(B*Cp_NH3 + C*Cp_N2 + D*Cp_H2)

c_star = sqrt(Rp * Tc) / GAMMA

epsilon = GAMMA / (((1 / Pi_tobera)**(1 / gamma)) * sqrt((2 * gamma / (gamma - 1)) * (1 - ((1 / Pi_tobera)**((gamma - 1) / gamma)))))

Ce = GAMMA * sqrt((2 * gamma / (gamma - 1)) * (1 - ((1 / Pi_tobera)**((gamma - 1) / gamma)))) + epsilon * (1 / Pi_tobera)

# AREAS GARGANTA
Ag1 = E1 / (Pc * Ce)
Ag2 = E2 / (Pc * Ce)

# GASTOS MASICOS 
m1 = Pc * Ag1 / c_star
m2 = Pc * Ag2 / c_star

print("Temperatura de cámara Tc = ", Tc)
print("Masa molar media (g/mol) = ", M_prod)
print("Cp", Cp_prod)
print("Cp_massico = ", Cp_prod_m)
print("gamma = ", gamma)
print("c* = ", c_star)
print()
print("Ag1 = ", Ag1, " m^2")
print("Ag2 = ", Ag2, " m^2")
print()
print("m1 = ", m1, " kg/s")
print("m2 = ", m2, " kg/s")

### APARTADO B)
P_tonk = 1.10 * Pc
V_tonk = Mf_N2H4 / rho_N2H4

print("V tonk", V_tonk)

### APARTADO C)
alpha = 0.15
