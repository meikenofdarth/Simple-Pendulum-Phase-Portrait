# pendulum_phase.py
import numpy as np
import matplotlib.pyplot as plt


# Parameters
g = 9.81    # gravity (m/s^2)
L = 1.0     # length (m)
h = 0.01    # time step
t_max = 20  # total time to integrate


def pendulum_rhs(t, y, g=g, L=L, b=0.0, A=0.0, Omega=0.0):

    # Returns dy/dt for y = [theta, omega].
    theta, omega = y
    dtheta = omega
    domega = - (g / L) * np.sin(theta) - b * omega + A * np.sin(Omega * t)
    return np.array([dtheta, domega])

def rk4_step(f, t, y, h):
    # One RK4 step for dy/dt = f(t, y).
    k1 = f(t, y)
    k2 = f(t + 0.5*h, y + 0.5*h*k1)
    k3 = f(t + 0.5*h, y + 0.5*h*k2)
    k4 = f(t + h, y + h*k3)
    return y + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)

def integrate_rk4(y0, t0, tmax, h, f):
    # Integrate using RK4 from t0 to tmax starting at y0.
    steps = int(np.ceil((tmax - t0) / h))
    t_vals = np.linspace(t0, t0 + steps*h, steps+1)
    ys = np.zeros((steps+1, 2))
    ys[0] = y0
    for i in range(steps):
        ys[i+1] = rk4_step(lambda tt, yy: f(tt, yy), t_vals[i], ys[i], h)
    return t_vals, ys

def wrap_theta(theta):
    return (theta + np.pi) % (2*np.pi) - np.pi

# Prepare grid for vector field and energy contours 
theta_min, theta_max = -2*np.pi, 2*np.pi  
omega_min, omega_max = -6.0, 6.0
N_theta, N_omega = 41, 41

theta_vals = np.linspace(theta_min, theta_max, N_theta)
omega_vals = np.linspace(omega_min, omega_max, N_omega)
TH, OM = np.meshgrid(theta_vals, omega_vals)

# Vector field: dtheta/dt = omega, domega/dt = - (g/L) * sin(theta)
U = OM
V = - (g / L) * np.sin(TH)

# Normalize vectors 
speed = np.sqrt(U**2 + V**2)
speed[speed == 0] = 1.0
U_norm = U / speed
V_norm = V / speed

# Energy (for contours)
E = 0.5 * OM**2 + (g / L) * (1 - np.cos(TH))
E_sep = 2 * g / L  # separatrix energy when L=1 -> E_sep=2*g

# ---- Plotting ----
fig, ax = plt.subplots(figsize=(9,6))

# Streamplot (smooth flow lines)
ax.streamplot(TH, OM, U, V, density=1.0, linewidth=0.7, arrowsize=1)


# ax.quiver(TH, OM, U_norm, V_norm, pivot='mid', alpha=0.6, scale=20)
ax.quiver(TH, OM, U_norm, V_norm, pivot='mid', alpha=0.5, scale=40, width=0.003)


# Overlay trajectories for multiple initial conditions
initial_conditions = [
    (0.2, 0.0),
    (0.5, 0.0),
    (1.5, 0.0),
    (2.8, 0.0),
    (np.pi - 0.3, 0.0),  # near the unstable inverted point
    (0.0, 2.5),         
    (-1.0, 1.0)
]

colors = plt.cm.viridis(np.linspace(0,1,len(initial_conditions)))
for (theta0, omega0), c in zip(initial_conditions, colors):
    t_vals, ys = integrate_rk4(np.array([theta0, omega0]), 0.0, t_max, h, pendulum_rhs)
    thetas = wrap_theta(ys[:,0])   # compact plotting
    omegas = ys[:,1]
    ax.plot(thetas, omegas, '-', linewidth=1.2, color=c, label=f"$\\theta_0$={theta0:.2f}, $\\omega_0$={omega0:.2f}")

# Energy contours (including separatrix)
levels = np.linspace(0, 8, 25)
cs = ax.contour(TH, OM, E, levels=levels, linewidths=0.6, alpha=0.6)
# highlight separatrix
ax.contour(TH, OM, E, levels=[E_sep], colors='red', linewidths=2.0, linestyles='--', label='separatrix')

# Mark equilibria
ax.plot(0, 0, 'ko', label='stable equilibrium')
ax.plot(np.pi, 0, 'rx', label='unstable equilibrium')

ax.set_xlim(theta_min, theta_max)
ax.set_ylim(omega_min, omega_max)
ax.set_xlabel(r'$\theta$ (rad)')
ax.set_ylabel(r'$\omega$ (rad/s)')
ax.set_title('Phase Portrait of the Simple Pendulum (RK4 trajectories + field + energy contours)')
ax.grid(True)
ax.legend(loc='upper right', fontsize='small', ncol=1)
plt.tight_layout()
plt.show()

# # ---- Optional: energy vs time for one trajectory to check conservation ----
# t_vals, ys = integrate_rk4(np.array([0.8, 0.0]), 0.0, t_max, h, pendulum_rhs)
# thetas = ys[:,0]; omegas = ys[:,1]
# energy_t = 0.5 * omegas**2 + (g / L) * (1 - np.cos(thetas))
# plt.figure(figsize=(6,3))
# plt.plot(t_vals, energy_t)
# plt.title("Energy vs time (RK4) for one trajectory")
# plt.xlabel("t")
# plt.ylabel("E(t)")
# plt.grid(True)
# plt.tight_layout()
# plt.show()
