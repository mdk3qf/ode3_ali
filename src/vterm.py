from scipy.integrate import solve_ivp
import numpy as np
from math import sqrt
from matplotlib import pyplot as plt

g = 9.81

# air resistance coefficient
air_k = 0.12

def hit_terminal(t, y):
    vx = y[1]
    vy = y[3]
    v = sqrt(vx*vx + vy*vy)
    drag = air_k*v*abs(v)
    # Trigger when drag is within threshold of weight
    return abs(drag - y[4]*g) - 1e-4
hit_terminal.terminal = True
hit_terminal.direction = -1  # trigger when crossing from above threshold to below

def func(t,y):
    m = y[4]  # append mass into state vector
    vx = y[1]
    vy = y[3]
    v  = sqrt(vx*vx + vy*vy)

    f0 = vx
    f1 = -air_k*v*vx/m
    f2 = vy
    f3 = -air_k*v*vy/m - g
    f4 = 0.0                     # mass
    return [f0,f1,f2,f3,f4]

def compute_vt(mass):
    # initial conditions
    y0 = [0, 10, 2000, 10, mass]

    # integrate
    sol = solve_ivp(func, [0,1000], y0, events=[hit_terminal], max_step=0.05)
    # returns (x, vx, y, vy, mass) at each time step

    if len(sol.t_events[0]) == 0:
        # Check final state
        vx_final = sol.y[1,-1] #(x, vx, y, vy, mass)
        vy_final = sol.y[3,-1]
        v_final = sqrt(vx_final**2 + vy_final**2)
        drag_final = air_k * v_final * abs(v_final)
        print(f"m={mass:.3f} kg, NO EVENT: drag={drag_final:.3f}, mg={mass*g:.3f}, v={v_final:.3f}")

    # extract vx,vy at event
    if len(sol.t_events[0]) == 0:
        return None

    tfinal = sol.t_events[0][0]
    vx = sol.y[1,sol.t == tfinal][0] #(x, vx, y, vy, mass)
    vy = sol.y[3,sol.t == tfinal][0]
    vt = sqrt(vx*vx + vy*vy)
    return vt

# sweep masses from 0.001 to 10 kg
masses = np.logspace(-3, 1, 30)
vts = []

for m in masses:
    vt = compute_vt(m)
    vts.append(vt)
    print(f"m={m:.3f} kg, vt={vt:.3f} m/s")

# plot
plt.plot(masses, vts, '-o', color='black', markersize=3)
plt.xscale('log')
plt.xlabel("Mass [kg]")
plt.ylabel("Terminal velocity [m/s]")
plt.title("Terminal velocity vs Mass")
plt.xscale('log')
plt.grid(alpha=0.5)
plt.tight_layout()
plt.savefig("vt_vs_mass.pdf")