from numpy import log, sqrt, linspace, zeros

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

T_D0 = 294 

# GASES
gases = ["H2", "He", "Aire", "CO2", "Ar"]
R_gases = ([4124, 2076.9, 287, 188.9, 208.1]) #sacados del manual del overleaf
gamma_gases = ([1.41, 1.66, 1.40, 1.30, 1.67]) #nist chemistry webbook

# MATERIALES
materiales = ["Al7075T6", "Acero", "Titanio6AL4V"]
rho_materiales = [2810, 8190, 4506]  # mirar            
sigma_u_materiales = [0.45e9, 0.54e9, 1e9]       

P_D_min = (1 + alpha) * P_tonk

P_D_max = 28e5 #cambiar
paso = 100
P_D_vector = linspace(13e5 , P_D_max, paso) 

V_D_z = zeros((len(gases), len(materiales), paso)) # Volúmen del depósito
W_g_z = zeros((len(gases), len(materiales), paso))  # Masa del gas presurizado
W_m_z = zeros((len(gases), len(materiales), paso))  # Masa del depósito 
W_d_z = zeros((len(gases), len(materiales), paso)) # Masa total del sitema de presurizado


for i, gas in enumerate(gases):
    print("\n===================================")
    print(f"          GAS: {gas}")
    print("===================================\n")

    for j, P_D in enumerate(P_D_vector):

        # Volumen del depósito
        V_D = (gamma_gases[i] * P_tonk * V_tonk) / (P_D - (1 + alpha) * P_tonk)

        # Masa del gas
        W_g = (P_D * V_D) / (R_gases[i] * T_D0)

        print(f"\n---- Presión de diseño P_D = {P_D/1e5:.2f} bar ----")
        print(f"Volumen V_D = {V_D:.6f} m^3")
        print(f"Masa del gas W_g = {W_g:.6f} kg")

        for w, mat in enumerate(materiales):

            # Masa del depósito
            W_m = (3/2) * (P_D * V_D * rho_materiales[w]) / sigma_u_materiales[w]

            # Masa total
            W_d = W_g + W_m

            # Guardar valores
            V_D_z[i, w, j] = V_D
            W_g_z[i, w, j] = W_g
            W_m_z[i, w, j] = W_m
            W_d_z[i, w, j] = W_d

            Isp = (E1 + E2) / (m1 + m2)
            Delta_V = Isp * log((W_d+Mcp+Mf_N2H4)/(W_d+Mcp))

            print(f"\nMaterial: {mat}")
            print(f"  Masa del depósito W_m = {W_m:.6f} kg")
            print(f"  Masa total W_d = {W_d:.6f} kg")
            print(f"  Delta V = {Delta_V:.6f} m/s")


import pandas as pd

# Crear tabla con los valores calculados
rows = []

for i, gas in enumerate(gases):
    for j, P_D in enumerate(P_D_vector):

        V_D = (gamma_gases[i] * P_tonk * V_tonk) / (P_D - (1 + alpha) * P_tonk)
        W_g = (P_D * V_D) / (R_gases[i] * T_D0)

        for w, mat in enumerate(materiales):

            W_m = (3/2) * (P_D * V_D * rho_materiales[w]) / sigma_u_materiales[w]
            W_d = W_g + W_m

            Isp = (E1 + E2) / (m1 + m2)
            Delta_V = Isp * log((W_d + Mcp + Mf_N2H4) / (W_d + Mcp))

            rows.append({
                "Gas": gas,
                "P_D_bar": P_D/1e5,
                "Material": mat,
                "V_D_m3": V_D,
                "W_g_kg": W_g,
                "W_m_kg": W_m,
                "W_d_kg": W_d,
                "DeltaV_m_s": Delta_V
            })

# Convertir a DataFrame
df = pd.DataFrame(rows)

# Guardar el CSV SIEMPRE con este nombre
df.to_csv("tabla_valores.csv", index=False)

print("\nCSV generado: tabla_valores.csv")



import matplotlib.pyplot as plt

# Crear graficas: una por material
for w, mat in enumerate(materiales):

    plt.figure(figsize=(8,6))
    
    # Para cada gas, extraemos las curvas correspondientes al material w
    for i, gas in enumerate(gases):

        Delta_V_vector = []

        for j, P_D in enumerate(P_D_vector):
            V_D = (gamma_gases[i] * P_tonk * V_tonk) / (P_D - (1 + alpha) * P_tonk)
            W_g = (P_D * V_D) / (R_gases[i] * T_D0)
            W_m = (3/2) * (P_D * V_D * rho_materiales[w]) / sigma_u_materiales[w]
            W_d = W_g + W_m

            Isp = (E1 + E2) / (m1 + m2)
            Delta_V = Isp * log((W_d + Mcp + Mf_N2H4) / (W_d + Mcp))

            Delta_V_vector.append(Delta_V)

        # Curva sin puntos → solo línea
        plt.plot(P_D_vector/1e5, Delta_V_vector, label=gas)

    # Acabado de la gráfica
    plt.xlabel("Presión de diseño $P_D$ [bar]")
    plt.ylabel("Incremento de velocidad $\Delta V$ [m/s]")
    plt.title(f"ΔV vs P_D para material: {mat}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.savefig(f"grafica_{mat}.png", dpi=300)

plt.show()

