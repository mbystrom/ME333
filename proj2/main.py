from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plot
from math import sqrt

v_1 = 320          # m/s
P_1 = 30.15 * 1000 # Pa
T_1 = 226.4        # K

k_air = 1.4

N_c = 0.73
N_t = 0.76

EPR_max = 18
T_max = 2450 # K

HHV = 46.2 # MJ/kg

thrust = 150 # kN

T_dead = 300 # K
P_dead = 100 # kPa


c_1 = PropsSI('A','T',T_1,'P',P_1,'air') # m/s
M_1 = v_1 / c_1

# stagnation temp/pressure from isentropic flow equations
T0 = T_1 * (1 + ((k_air-1)/2)*(M_1**2))                     # K
P0 = P_1 * (1 + ((k_air-1)/2)*(M_1**2))**(k_air/(k_air-1))  # Pa

h_1 = PropsSI('H','T',T_1,'P',P_1,'air') # J/kg
s_1 = PropsSI('S','T',T_1,'P',P_1,'air') # J/kg-K

T_2 = T0 # K
P_2 = P0 # Pa

h_2 = PropsSI('H','T',T_2,'P',P_2,'air') # J/kg
s_2 = PropsSI('S','T',T_2,'P',P_2,'air') # J/kg-K

P_3 = P_2 * EPR_max                       # Pa
h_3s = PropsSI('H','P',P_3,'S',s_2,'air') # J/kg
T_3s = PropsSI('T','P',P_3,'S',s_2,'air') # K

h_3 = h_2 + (h_3s - h_2)/N_c # J/kg
s_3 = PropsSI('S','H',h_3,'P',P_3,'air') # J/kg-K
T_3 = PropsSI('T','H',h_3,'P',P_3,'air') # K

T_4 = T_max # K
P_4 = P_3   # Pa
h_4 = PropsSI('H','T',T_4,'P',P_4,'air') # J/kg
s_4 = PropsSI('S','T',T_4,'P',P_4,'air') # J/kg-K

# use known BWR of 1.0 to find h_5
BWR = 1
w_in = h_3 - h_2    # W/kg
w_out = w_in * BWR  # W/kg
h_5 = h_4 - w_out   # W/kg

# find h_5s from h_5 & N_t, use that and known entropy at h_5s to find P_5
h_5s = h_4 - (h_4-h_5)/N_t                  # J/kg
P_5s = PropsSI('P','H',h_5s,'S',s_4,'air')  # Pa

P_5 = P_5s # Pa

T_5  = PropsSI('T','H',h_5,'P',P_5,'air')  # K
s_5  = PropsSI('S','T',T_5,'P',P_5,'air')  # J/kg-K

P_6 = P_1 # Pa

# isentropic flow to output
s_6 = s_5
h_6 = PropsSI('H','S',s_6,'P',P_6,'air') # J/kg
T_6 = PropsSI('T','S',s_6,'P',P_6,'air') # K

v_6 = sqrt(2*((h_5 - h_6))) # m/s
dv = v_6 - v_1              # m/s
m_a = thrust*1000 / dv      # kg/s

m_f = (m_a * (h_4 - h_3)) / (HHV*(10**6)) # kg/s
TSFC = m_f / thrust                       # kg/s/kN

c_6 = PropsSI('A','T',T_6,'P',P_6,'air') # m/s
M_6 = v_6 / c_6
#print(f"mach @6: {M_6}")


# EXERGY

s_c   = abs(s_2 - s_3) * m_a # J/K
X_d_c = s_c * T_dead         # W

s_t   = abs(s_4 - s_5) * m_a # J/kg-K
X_d_t = s_t * T_dead         # W

# Area

# m_dot = rho * A * v -> A = m_dot / (rho * v)

rho_1 = PropsSI('D','T',T_1,'P',P_1,'air') # kg/m3
rho_2 = PropsSI('D','T',T_2,'P',P_2,'air') # kg/m3
rho_3 = PropsSI('D','T',T_3,'P',P_3,'air') # kg/m3
rho_4 = PropsSI('D','T',T_4,'P',P_4,'air') # kg/m3
rho_5 = PropsSI('D','T',T_5,'P',P_5,'air') # kg/m3
rho_6 = PropsSI('D','T',T_6,'P',P_6,'air') # kg/m3

A_1 = m_a / (rho_1 * v_1) # m2
A_6 = m_a / (rho_6 * v_6) # m2

# need specific volume, not velocity!!!!
print(f"state  1 - T: {T_1-273.15:.2f}\N{DEGREE SIGN}C,  P: {P_1/1000:.2f}kPa,   h: {h_1/1000:.2f}kJ/kg,  s: {s_1/1000:.3f}kJ/kg-K  v: {1/rho_1:.3f}m3/kg")
print(f"state  2 - T: {T_2-273.15:.2f}\N{DEGREE SIGN}C,    P: {P_2/1000:.2f}kPa,   h: {h_2/1000:.2f}kJ/kg,  s: {s_2/1000:.3f}kJ/kg-K  v: {1/rho_2:.3f}m3/kg")
print(f"state 3s - T: {T_3s-273.15:.2f}\N{DEGREE SIGN}C,  P: {P_3/1000:.2f}kPa, h: {h_3s/1000:.2f}kJ/kg,  s: {s_2/1000:.3f}kJ/kg-K  v: -----m3/kg")
print(f"state  3 - T: {T_3-273.15:.2f}\N{DEGREE SIGN}C,  P: {P_3/1000:.2f}kPa, h: {h_3/1000:.2f}kJ/kg,  s: {s_3/1000:.3f}kJ/kg-K  v: {1/rho_3:.3f}m3/kg")
print(f"state  4 - T: {T_4-273.15:.2f}\N{DEGREE SIGN}C, P: {P_4/1000:.2f}kPa, h: {h_4/1000:.2f}kJ/kg, s: {s_4/1000:.3f}kJ/kg-K  v: {1/rho_4:.3f}m3/kg")
print(f"state 5s - T: ------\N{DEGREE SIGN}C,  P: {P_5/1000:.2f}kPa,  h: {h_5s/1000:.2f}kJ/kg, s: {s_4/1000:.3f}kJ/kg-K  v: -----m3/kg")
print(f"state  5 - T: {T_5-273.15:.2f}\N{DEGREE SIGN}C, P: {P_5/1000:.2f}kPa,  h: {h_5/1000:.2f}kJ/kg, s: {s_5/1000:.3f}kJ/kg-K  v: {1/rho_5:.3f}m3/kg")
print(f"state  6 - T: {T_6-273.15:.2f}\N{DEGREE SIGN}C,  P: {P_6/1000:.2f}kPa,   h: {h_6/1000:.2f}kJ/kg, s: {s_6/1000:.3f}kJ/kg-K  v: {1/rho_6:.3f}m3/kg")

print()

print(f"mass flow rate (air): {m_a:.4f}kg/s")
print(f"mass flow rate (fuel): {m_f:.5f}kg/s")
print(f"nozzle velocity: {v_6:.2f}m/s")
print(f"turbine power produced: {(w_out*m_a)/1000:.2f}kW")
print(f"compressor exergy destruction: {X_d_c/1000:.2f}kW")
print(f"turbine exergy destruction: {X_d_t/1000:.2f}kW")
print(f"thrust-specific fuel consumption: {TSFC:.8f}kg/kN-s")

print()

print(f"diffuser inlet area: {A_1:.5f}m^2")
print(f"nozzle outlet area: {A_6:.5f}m^2")

print()

###########################################################################################################

cycleT = [T_1,T_2,T_3,T_4,T_5,T_6]
cycles = [s_1,s_2,s_3,s_4,s_5,s_6]
plot.subplot(2,4,1)
plot.plot(cycles,cycleT,label='Cycle T-s',color='r')
plot.plot(cycles,cycleT,marker='o',color='r')
plot.xlabel("s (J/kg-K)")
plot.ylabel("T (K)")
plot.legend()
plot.title("T-s Diagram")
plot.grid()

###########################################################################################################

# h_5 won't change since the compressor always needs the same power, but h_5s and P_5 will change, which changes output velocity and mass flow rate

TSFC_n = []
X_d_n = []
n = []

for i in range(50,101):

    try:
        
        h5sn = h_4 - (h_4-h_5)/(i/100)
        P_5n = PropsSI('P','H',h5sn,'S',s_4,'air')
        
        s_5n = PropsSI('S','T',T_5,'P',P_5n,'air')
        h_6n = PropsSI('H','S',s_5n,'P',P_6,'air')
        v_6n = sqrt(2*((h_5 - h_6n)))

        m_an = (thrust*1000) / (v_6n - v_1)
        m_fn = (m_an * (h_4 - h_3)) / (HHV*(10**6))
        X_dn = (T_dead * m_an * (s_4 - s_5n)) / 1000

        #print(f"{i}% efficiency: h5sn: {h5sn/1000:.3f}kJ/kg, P_5n: {P_5n/1000:.3f}kJ/kg-K, h6n: {h_6n/1000:.3f}kJ/kg, v_6n: {v_6n:.2f}m/s")
        #print(f"{i}% efficiency: m_an: {m_an:.5f}, TSFC: {m_fn / 1000}")

        n.append(i)
        TSFC_n.append(m_fn / thrust)
        X_d_n.append(X_dn)
    except:
        print(f"ERROR ({i}% efficiency): not compatible with coolprop library (Equations of State for h/s only valid to 2000K)")

plot.subplot(2,4,2)
plot.plot(n,TSFC_n,label='TSFC')
plot.xlabel("Turbine Isentropic Efficiency (%)")
plot.ylabel("Thrust-Specific Fuel Consumption (kg/s/kN)")
plot.title("TSFC vs. Turbine Efficiency")
plot.legend()
plot.grid()
plot.subplot(2,4,6)
plot.plot(n,X_d_n,label='X_d')
plot.xlabel("Turbine Isentropic Efficiency (%)")
plot.ylabel("Exergy Destroyed (kW)")
plot.title("Turbine Exergy Destruction vs. Turbine Efficiency")
plot.legend()
plot.grid()

###########################################################################################################

TSFC_P = []
m_aP = []
Pr = []

for i in range(10,26):
    try:
        P_3p = P_2 * i
        h3sp = PropsSI('H','P',P_3p,'S',s_2,'air') # J/kg
        T3sp = PropsSI('T','P',P_3p,'S',s_2,'air') # K

        h_3p = h_2 + (h3sp - h_2)/N_c # J/kg
        s_3p = PropsSI('S','H',h_3p,'P',P_3p,'air') # J/kg-K
        T_3p = PropsSI('T','H',h_3p,'P',P_3p,'air') # K

        h_4p = PropsSI('H','T',T_4,'P',P_3p,'air') # J/kg
        s_4p = PropsSI('S','T',T_4,'P',P_3p,'air') # J/kg-K

        w_inp = h_3p - h_2    # W/kg
        w_outp = w_in * BWR  # W/kg
        h_5p = h_4p - w_out   # W/kg

        h_5sp = h_4p - (h_4p-h_5p)/N_t                 # J/kg
        P_5sp = PropsSI('P','H',h_5sp,'S',s_4p,'air')  # Pa

        s_5p = PropsSI('S','H',h_5p,'P',P_5sp,'air')
        h_6p = PropsSI('H','S',s_5p,'P',P_6,'air')
        v_6p = sqrt(2*((h_5p - h_6p)))

        m_a_p = (thrust*1000) / (v_6p - v_1)
        m_fp = (m_a_p * (h_4p - h_3p)) / (HHV*(10**6))

        Pr.append(i)
        TSFC_P.append(m_fp / thrust)
        m_aP.append(m_a_p)

    except:
        print(f"ERROR (pressure ratio {i}): not compatible with coolprop library")

plot.subplot(2,4,3)
plot.plot(Pr,TSFC_P,label="TSFC")
plot.xlabel("Pressure Ratio")
plot.ylabel("Thrust-Specific Fuel Consumption (kg/s/kN)")
plot.legend()
plot.title("TSFC vs. Pressure Ratio")
plot.grid()
plot.subplot(2,4,7)
plot.plot(Pr,m_aP,label="m_a")
plot.xlabel("Pressure Ratio")
plot.ylabel("Mass Flow Rate, Air (kg/s)")
plot.legend()
plot.title("Air Flowrate vs. Pressure Ratio")
plot.grid()

###########################################################################################################

TSFC_T = []
X_d_T = []
maxT = []

for i in range(1500,3050,50):
    try:
        T_4T = i # K
        P_4T = P_3   # Pa
        h_4T = PropsSI('H','T',T_4T,'P',P_4T,'air') # J/kg
        s_4T = PropsSI('S','T',T_4T,'P',P_4T,'air') # J/kg-K

        w_inT = h_3 - h_2    # W/kg
        w_out = w_in * BWR  # W/kg
        h_5T = h_4T - w_out   # W/kg

        h5sT = h_4T - (h_4T-h_5T)/N_t
        P_5T = PropsSI('P','H',h5sT,'S',s_4T,'air')

        h5sT = h_4T - (h_4T-h_5T)/N_t
        P_5T = PropsSI('P','H',h5sT,'S',s_4T,'air')
        
        s_5T = PropsSI('S','H',h_5T,'P',P_5T,'air')
        h_6T = PropsSI('H','S',s_5T,'P',P_6,'air')
        #print(f"h_5T: {h_5T:.2f}J/kg, h_6T: {h_6T:.2f}J/kg")
        v_6T = sqrt(2*((h_5T - h_6T)))

        m_aT = (thrust*1000) / (v_6T - v_1)
        m_fT = (m_aT * (h_4T - h_3)) / (HHV*(10**6))
        X_dT = (T_dead * m_aT * (s_4T - s_5T)) / 1000

        maxT.append(i)
        TSFC_T.append(m_fT / thrust)
        X_d_T.append(X_dT)

    except:
        print(f"ERROR ({i}K): inlet temperature not compatible with coolprop library (Equations of State for h/s only valid to 2000K)")

plot.subplot(2,4,4)
plot.plot(maxT,TSFC_T,label='TSFC')
plot.xlabel("Turbine Inlet Temperature (K)")
plot.ylabel("Thrust-Specific Fuel Consumption (kg/s/kN)")
plot.title("TSFC vs. Turbine Inlet Temp")
plot.legend()
plot.grid()
plot.subplot(2,4,8)
plot.plot(maxT,X_d_T,label='X_d')
plot.xlabel("Turbine Inlet Temperature (K)")
plot.ylabel("Exergy Destroyed (kW)")
plot.title("Turbine Exergy Destruction vs. Turbine Inlet Temp")
plot.legend()
plot.grid()

###########################################################################################################

plot.show()

