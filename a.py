import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

e = 1.602e-19
kB = 1.38e-23

st.title("Real-Time Single Electron Transistor (SET) Simulator")
st.subheader("Coulomb Blockade • Stability Diagram • Monte Carlo Tunneling")

# User Inputs
st.sidebar.header("SET Parameters")
Cs = st.sidebar.slider("Source Capacitance Cs (aF)", 1.0, 50.0, 10.0) * 1e-18
Cd = st.sidebar.slider("Drain Capacitance Cd (aF)", 1.0, 50.0, 10.0) * 1e-18
Cg = st.sidebar.slider("Gate Capacitance Cg (aF)", 1.0, 50.0, 10.0) * 1e-18
Rt = st.sidebar.slider("Tunnel Resistance Rt (MΩ)", 1.0, 200.0, 100.0) * 1e6
Temp = st.sidebar.slider("Temperature (K)", 0.1, 10.0, 1.5)

Ctot = Cs + Cd + Cg
Ec = e**2 / (2 * Ctot)

st.write(f"### Charging Energy: {Ec/e:.3f} eV")

# --------------------------------------------
# 1. Coulomb Oscillations
# --------------------------------------------
st.header("1. Coulomb Blockade Oscillations")

Vg = np.linspace(-0.02, 0.02, 500)
n_values = (Cg*Vg/e)
conductance = np.cos(2*np.pi*n_values)**2

fig1, ax1 = plt.subplots()
ax1.plot(Vg, conductance)
ax1.set_xlabel("Gate Voltage (V)")
ax1.set_ylabel("Conductance (arb.)")
ax1.set_title("Conductance Oscillations (Coulomb Blockade)")
st.pyplot(fig1)

# --------------------------------------------
# 2. I–V Characteristics
# --------------------------------------------
st.header("2. I–V Characteristics under Coulomb Blockade")

Vd = np.linspace(-0.01, 0.01, 400)
I = []

for V in Vd:
    dF = e*V - Ec
    rate = (dF / (e**2 * Rt)) / (1 - np.exp(-dF / (kB * Temp)))
    I.append(rate * e)

I = np.nan_to_num(I)

fig2, ax2 = plt.subplots()
ax2.plot(Vd, I)
ax2.set_xlabel("Drain Voltage (V)")
ax2.set_ylabel("Current (A)")
ax2.set_title("I–V Characteristics")
st.pyplot(fig2)

# --------------------------------------------
# 3. Coulomb Diamond Stability Diagram
# --------------------------------------------
st.header("3. Coulomb Diamond (Stability Diagram)")

Vg_grid = np.linspace(-0.02, 0.02, 200)
Vd_grid = np.linspace(-0.02, 0.02, 200)

diamond = np.zeros((200, 200))

for i, vg in enumerate(Vg_grid):
    for j, vd in enumerate(Vd_grid):
        n = int((Cg*vg)/e)
        F = Ec * (n - Cg*vg/e)**2 + e*vd
        diamond[j, i] = F

fig3, ax3 = plt.subplots()
ax3.imshow(diamond, extent=[Vg_grid.min(), Vg_grid.max(),
                            Vd_grid.min(), Vd_grid.max()],
           origin="lower", aspect="auto")
ax3.set_xlabel("Gate Voltage (V)")
ax3.set_ylabel("Drain Voltage (V)")
ax3.set_title("Coulomb Diamond Stability Diagram")
st.pyplot(fig3)

# --------------------------------------------
# 4. Monte Carlo Tunneling Simulation
# --------------------------------------------
st.header("4. Monte Carlo Single Electron Tunneling")

steps = 500
prob = 0.5
electrons = []
count = 0

for _ in range(steps):
    dF = Ec * (np.random.rand() - 0.5)
    rate = abs(dF)/(Rt * e**2)
    if np.random.rand() < min(rate*1e-12, 1):
        count += 1
    electrons.append(count)

fig4, ax4 = plt.subplots()
ax4.plot(electrons)
ax4.set_xlabel("Time Steps")
ax4.set_ylabel("Tunneled Electrons")
ax4.set_title("Monte Carlo Electron Tunneling")
st.pyplot(fig4)

st.success("Simulation Complete — Adjust parameters from sidebar to see effects in real time!")
