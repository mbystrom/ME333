# to use coolprop, must run "pip install coolprop" from terminal

from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plot

N_t = 0.8
N_p = 0.85

p_a = 101325 #ambient pressure, Pa

T_1 = 60+273.15 #K
x_1 = 1
T_2 = 28+273.15 #K
x_3 = 0

T_5 = 70+273.15 #K
T_6 = 65+273.15 #K
T_7 = 20+273.15 #K
T_8 = 25+273.15 #K

h_1 = PropsSI('H','T',T_1,'Q',x_1,'propane') #J/kg
p_1 = PropsSI('P','T',T_1,'Q',x_1,'propane') #Pa
p_4 = p_1

s_1  = PropsSI('S','T',T_1,'Q',x_1,'propane')  #J/kg-K
h_2s = PropsSI('H','T',T_2,'S',s_1,'propane')  #J/kg
T_2s = PropsSI('T','H',h_2s,'S',s_1,'propane') #K
h_2  = h_1 - N_t*(h_1-h_2s)                    #J/kg


p_2t = PropsSI('P','T',T_2,'Q',0,'propane') #Pa
h_2f = PropsSI('H','T',T_2,'Q',0,'propane') #sat. fluid
h_2g = PropsSI('H','T',T_2,'Q',1,'propane') #sat. vapor

#p_2  = PropsSI('P','T',T_2,'H',h_2,'propane') #Pa #T/h call doesn't work
if (h_2 >= h_2f and h_2 <= h_2g): #PropsSI call with T/h call doesn't work; check that result is saturated, otherwise not enough info to solve
    p_2 = p_2t
else:
    print("ERROR: Superheated Turbine Exit - CoolProp Cannot retrieve")
    exit()

x_2 = (h_2 - h_2f) / (h_2g - h_2f)
s_2 = PropsSI('S','T',T_2,'Q',x_2,'propane')

p_3 = p_2
h_3 = PropsSI('H','P',p_3,'Q',x_3,'propane') #J/kg
s_3 = PropsSI('S','P',p_3,'Q',x_3,'propane') #J/kg-K
T_3 = PropsSI('T','P',p_3,'Q',x_3,'propane') #K

h_4s = PropsSI('H','P',p_4,'S',s_3,'propane')  #J/kg
h_4  = h_3 + (h_4s - h_3)/N_p                  #J/kg
s_4  = PropsSI('S','P',p_4,'H',h_4,'propane')  #J/kg-K
T_4  = PropsSI('T','P',p_4,'H',h_4,'propane')  #K
T_4s = PropsSI('T','P',p_4,'H',h_4s,'propane') #K

w_out = (h_1-h_2)/1000 #kW/kg
w_in  = (h_3-h_4)/1000 #kW/kg
w_net = w_out + w_in   #kW/kg
q_in  = (h_1-h_4)/1000 #kW/kg
q_out = (h_2-h_3)/1000 #kW/kg

m_p = 20/w_net #mass flow rate propane, kg/s
Q_in  = q_in * m_p  #kW
Q_out = q_out * m_p #kW

h_5 = PropsSI('H','T',T_5,'P',p_a,'water')
h_6 = PropsSI('H','T',T_6,'P',p_a,'water')
h_7 = PropsSI('H','T',T_7,'P',p_a,'water')
h_8 = PropsSI('H','T',T_8,'P',p_a,'water')

m_h = Q_in / (h_5 - h_6)
m_c = Q_out / (h_8 - h_7)

N_carnot = (T_5-T_7)/T_5

print(f"state 1:  T: {T_1-273.15:.2f}C, P: {p_1/1000:.2f}kPa, h: {h_1/1000:.3f}kJ/kg, s: {s_1/1000:.3f}kJ/kg-K")
print(f"state 2:  T: {T_2-273.15:.2f}C, P: {p_2/1000:.2f}kPa, h: {h_2/1000:.3f}kJ/kg, s: {s_2/1000:.3f}kJ/kg-K")
print(f"state 2s: T: {T_2s-273.15:.2f}C, P: {p_2/1000:.2f}kPa, h: {h_2s/1000:.3f}kJ/kg, s: {s_1/1000:.3f}kJ/kg-K")
print(f"state 3:  T: {T_3-273.15:.2f}C, P: {p_3/1000:.2f}kPa, h: {h_3/1000:.3f}kJ/kg, s: {s_3/1000:.3f}kJ/kg-K")
print(f"state 4:  T: {T_4-273.15:.2f}C, P: {p_4/1000:.2f}kPa, h: {h_4/1000:.3f}kJ/kg, s: {s_4/1000:.3f}kJ/kg-K")
print(f"state 4s: T: {T_4s-273.15:.2f}C, P: {p_4/1000:.2f}kPa, h: {h_4s/1000:.3f}kJ/kg, s: {s_3/1000:.3f}kJ/kg-K")

print(f"")

print(f"turbine work: {w_out:.2f} kW/kg")
print(f"turbine entropy generation: {(s_2-s_1)/1000:.4f}kJ/kg-K")
print(f"pump work: {w_in:.2f} kW/kg")
print(f"heat input: {q_in:.2f} kW/kg")
print(f"heat rejection: {q_out:.2f} kW/kg")

print(f"")

print(f"thermal efficiency: {((w_out+w_in)/q_in)*100:.2f}%")

print(f"")

print(f"total heat in: {Q_in:.2f}kW")
print(f"total heat out: {Q_out:.2f}kW")
print(f"cycle mass flow rate: {m_p:.3f}kg/s")
print(f"hot side mass flow rate: {m_h:.4f}kg/s")
print(f"cold side mass flow rate: {m_c:.4f}kg/s")

print(f"")

print(f"max possible efficiency (carnot efficiency): {N_carnot*100:.2f}%")

print(f"")


###########################################################################################################

T0  = 273 #K
Tcrit = PropsSI('Tcrit','propane')
Pcrit = PropsSI('Pcrit','propane')

Tr = range(int(T0),int(Tcrit+1),1)

Sarr = [0]*len(Tr)*2
Tarr = [0]*len(Tr)*2

j = 0
for i in Tr:
    if i == int(Tcrit+1):
        Scrit = PropsSI('S','T',Tcrit,'P',Pcrit,'propane')
        Tarr[j] = i
        Sarr[j] = Scrit
    else:
        Sif = PropsSI('S','T',i,'Q',0,'propane')
        Sig = PropsSI('S','T',i,'Q',1,'propane')
        Tarr[j] = i
        Sarr[j] = Sif
        Tarr[len(Tr)*2-(j+1)] = i
        Sarr[len(Tr)*2-(j+1)] = Sig
    j = j+1

plot.subplot(1,3,1)
plot.plot(Sarr,Tarr,label="vapor dome")

s_inter = PropsSI('S','T',T_1,'Q',0,'propane') #add a point at T1,x=0 to approximate the isobar
cycleT = [T_1,T_2,T_3,T_4,T_1,T_1]
cycleS = [s_1,s_2,s_3,s_4,s_inter,s_1]
markerS = [s_1,s_2,s_3,s_4]
markerT = [T_1,T_2,T_3,T_4]

plot.plot(cycleS,cycleT,label="cycle T-s",color="r")
plot.plot(markerS,markerT,'o',color="r")
plot.xlabel("s (J/kg-K)")
plot.ylabel("T (K)")
plot.legend()
plot.title("T-s diagram")
plot.grid()

###########################################################################################################

T_1max = 90+273.15 #K
T_1r = range(int(T_1),int(T_1max),1)

Tvariance = []
N_arr = []

for i in T_1r:
    h = PropsSI('H','T',i,'Q',1,'propane')
    wnet = h-h_2 + w_in
    qin = h-h_4
    N = wnet / qin
    N_arr.append(N*100)
    Tvariance.append(i)

plot.subplot(1,3,2)
plot.plot(Tvariance,N_arr,label="thermal efficiency")
plot.xlabel("Temperature (K)")
plot.ylabel("Thermal Efficiency (%)")
plot.legend()
plot.title("Variance in Thermal Efficiency with Temperature")
plot.grid()

###########################################################################################################

P_1new = 2000 #kPa
hnew = PropsSI('H','P',(P_1new*(10**3)),'T',T_1max,'propane')
snew = PropsSI('S','P',(P_1new*(10**3)),'T',T_1max,'propane')
h2snew = PropsSI('H','S',snew,'P',p_2,'propane')
h2new = hnew - N_t*(hnew-h2snew)
Nnew = (hnew - h2new + w_in) / (hnew - h_4)
print(f"new h1, 90C, 2MPa: {hnew/1000:.2f}kJ/kg")
print(f"new h2 from new h1: {h2new/1000:.2f}kJ/kg")
print(f"thermal efficiency @ 90C, 2MPa: {Nnew*100:.2f}%")

Pvariance = []
N_arr2 = []
P_1max = 3500 #kPa
Pr = range(P_1new,P_1max+100,100)

for i in Pr:
    h = PropsSI('H','P',(i*(10**3)),'T',T_1max,'propane')
    s = PropsSI('S','P',(i*(10**3)),'T',T_1max,'propane')
    h2s = PropsSI('H','S',s,'P',p_2,'propane')
    h2 = h - N_t*(h-h2s)
    N = (h - h2 + w_in) / (h - h_4)
    Pvariance.append(i)
    N_arr2.append(N*100)

plot.subplot(1,3,3)
plot.plot(Pvariance,N_arr2,label="thermal efficiency")
plot.xlabel("Pressure (kPa)")
plot.ylabel("Thermal Efficiency (%)")
plot.legend()
plot.title("Variance in Thermal Efficiency with Pressure")
plot.grid()

plot.show()
