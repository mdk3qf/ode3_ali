from scipy.integrate import solve_ivp
import numpy as np
from math import sqrt

# note: params can also be passed to the function describing the system
#       see solve_ivp docs
# m, air_k specified as parameters
params=[1.0, 0.12]
#params=[1.0, 0.0]  # no air resistance

# "event function" used to terminate the integration
# this function terminates the integral on a zero crossing
# with a negative slope
def hit_ground(t,y):
    return y[2]
hit_ground.terminal = True
hit_ground.direction = -1


def func(t,y):
    g = 9.81
    m = params[0]
    air_k = params[1]
    v = sqrt(y[1]*y[1] + y[3]*y[3])
    f0 = y[1]                         # f_ri
    f1 = -air_k*v*y[1]/m              # f_vi
    f2 = y[3]                         # f_rj
    f3 = -air_k*v*y[3]/m - g          # f_vj
    return [f0,f1,f2,f3]

#starting coordinates
y0=[0,10,100,10]   # x=0, vx=10 m/s, y=0, vy=10 m/s

#t = np.linspace(0,1.6,200)
tlim=15
t = np.linspace(0,tlim,300)
#sol = solve_ivp(func, [0,tlim], y0, t_eval=t, events=[hit_terminal_velocity])
sol = solve_ivp(func, [0,tlim], y0, t_eval=t, events=[hit_ground])
yf=sol.y  # array of coordiantes at each time step
#print(f"Hitting ground at t = {sol.t_events[0][0]:.3f} seconds")
print(f"Hitting terminal velocity at t = {sol.t_events[0][0]:.3f} seconds")
#print(yf)
#for y in sol:
#    print(y)

from matplotlib import pyplot as plt

x=np.array(yf[0])
vx=np.array(yf[1])
y=np.array(yf[2])
vy=np.array(yf[ 3])
ke = 0.5*params[0]*(vx*vx+vy*vy)
pe = params[0]*9.81*y
E_tot = ke + pe

drag_acc_terminal = params[1]*np.sqrt(vx*vx + vy*vy)*vy
grav_acc = params[0]*9.81

for i in range(len(drag_acc_terminal)):
    if abs(drag_acc_terminal[i] + grav_acc) < 1e-4:
        terminal_velocity = sqrt(vx[i]**2 + vy[i]**2)
        print(f"Terminal velocity = {terminal_velocity:.3f} m/s at t = {sol.t[i]:.3f} seconds")
        break

plt.subplot(2, 2, 1)
plt.plot(x,y,'-')
plt.xlabel('x values [m]')
plt.ylabel('y values [m]')
plt.title('y vs x')

plt.subplot(2, 2, 2)
plt.plot(sol.t,ke,'-')
plt.xlabel('time [s]')
plt.ylabel('KE [J]')
plt.title('KE vs t')
plt.axhline(0.5*params[0]*(terminal_velocity**2), color='r', linestyle='--', label='$KE_{terminal}$')
plt.legend()

plt.subplot(2,2,3)
plt.plot(sol.t, pe, '-')
plt.xlabel('time [s]')
plt.ylabel('Potential Energy [J]')
plt.title(r'PE vs t')

plt.subplot(2,2,4)
plt.plot(sol.t, E_tot, '-')
plt.xlabel('time [s]')
plt.ylabel('Total Energy [J]')
#plt.yscale('linear')
plt.ylim(20, max(E_tot)*1.1)
plt.title(r'E$_{\text{tot}}$ vs t')

plt.tight_layout()
plt.savefig('vterm_pt1.png', dpi=300)
