import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import PchipInterpolator

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.size": 18,
    "axes.titlesize": 20,
    "axes.labelsize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 14,
    "figure.titlesize": 22
})

#Rayleigh
Ra = float(input("Insert Rayleigh number: "))
print(f"Loading files for Ra = {int(Ra)}...")

#Domain parameters
domainWidth = 1.0
aspectRatio = 4.0
fine_N = 500
x_fine = np.linspace(0, domainWidth, fine_N)
y_fine = np.linspace(0, domainWidth * aspectRatio, fine_N)
X_fine, Y_fine = np.meshgrid(x_fine, y_fine)


##################################################################

################## COLORMAPS - SEPARATE FIGURES ##################

##################################################################

quantities = [
    (f"TemperatureDistribution_Ra_{int(Ra)}.txt", "Temperature"),
    (f"VelocityUDistribution_Ra_{int(Ra)}.txt", "u-Velocity"),
    (f"VelocityVDistribution_Ra_{int(Ra)}.txt", "v-Velocity")
]

for filename, quantity_name in quantities:
    # Load data
    data = pd.read_csv(filename, sep=r"\s+", header=None)
    x, y, value = data[0].values, data[1].values, data[2].values

    # Interpolation
    value_fine = griddata((x, y), value, (X_fine, Y_fine), method='cubic')

    fig, ax = plt.subplots(figsize=(6, 10))
    img = ax.contourf(X_fine, Y_fine, value_fine, levels=100, cmap="jet")
    contours = ax.contour(X_fine, Y_fine, value_fine, colors='black', levels=30, linewidths=0.8)
    # ax.clabel(contours, inline=True, fontsize=8)
    ax.set_xlim([0, domainWidth])
    ax.set_ylim([0, domainWidth * aspectRatio])
    ax.set_aspect('equal')
    #ax.axis("off")
    #fig.suptitle(f"{quantity_name} Field", fontsize=22)
    cbar = fig.colorbar(img, ax=ax, shrink=0.9, pad=0.02)
    output_name = f"{quantity_name.replace('-', '')}_Colormap_Ra_{int(Ra)}.pdf"
    plt.tight_layout()
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


##############################################################

################## TURBULENT KINETIC ENERGY ##################

##############################################################

# === Load data
df = pd.read_csv(f"Turbulent_KE_{int(Ra)}.txt", sep=r"\s+", header=None,
                 names=["x", "y", "u_var", "v_var", "uv_fluct", "tke", "ke"])
points = df[["x", "y"]].values

# === Interpolation
fields = {
    "TKE": df["tke"].values,
    "KE": df["ke"].values,
    "u'u'": df["u_var"].values,
    "v'v'": df["v_var"].values,
    "u'v'": df["uv_fluct"].values
}

interpolated_fields = {
    name: griddata(points, values, (X_fine, Y_fine), method='cubic')
    for name, values in fields.items()
}

def plot_field(Z, filename, num_levels=14, cmap="Blues", force_positive=False):

    levels = np.linspace(np.nanmin(Z), np.nanmax(Z), num_levels)

    fig, ax = plt.subplots(figsize=(6, 10))
    cf = ax.contourf(X_fine, Y_fine, Z, levels=levels, cmap=cmap)
    ax.contour(X_fine, Y_fine, Z, levels=levels, colors='black', linewidths=0.8)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim([0, domainWidth])
    ax.set_ylim([0, domainWidth * aspectRatio])
    ax.set_aspect('equal')

    fig.colorbar(cf, ax=ax, shrink=0.9, pad=0.02)

    plt.tight_layout()
    plt.savefig(f"{filename}.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


plot_field(interpolated_fields["TKE"], f"TKE_Ra_{int(Ra)}")
plot_field(interpolated_fields["KE"],  f"KE_Ra_{int(Ra)}")
plot_field(interpolated_fields["u'u'"],  f"u_var_Ra_{int(Ra)}")
plot_field(interpolated_fields["v'v'"],  f"v_var_Ra_{int(Ra)}")
plot_field(interpolated_fields["u'v'"],  f"uv_fluct_Ra_{int(Ra)}")




########################################################################

################## STREAM FUNCTION - SEPARATE FIGURES ##################

########################################################################

df_psi = pd.read_csv(f"StreamFunction_Ra_{int(Ra)}.txt", sep=r"\s+", header=None, names=["x", "y", "psi"])
df_u = pd.read_csv(f"VelocityUDistribution_Ra_{int(Ra)}.txt", sep=r"\s+", header=None, names=["x", "y", "u"])
df_v = pd.read_csv(f"VelocityVDistribution_Ra_{int(Ra)}.txt", sep=r"\s+", header=None, names=["x", "y", "v"])
df_psi_sorted = df_psi.sort_values(by=["y", "x"], ascending=[True, True])
Z = df_psi_sorted.pivot(index="y", columns="x", values="psi").values

x_vals = np.sort(df_psi["x"].unique())
y_vals = np.sort(df_psi["y"].unique())
X, Y = np.meshgrid(x_vals, y_vals)

# === Interpolation
scale = 3
x_all = np.linspace(x_vals.min(), x_vals.max(), len(x_vals) * scale)
y_all = np.linspace(y_vals.min(), y_vals.max(), len(y_vals) * scale)
Xi, Yi = np.meshgrid(x_all, y_all)
Ui = griddata((df_u["x"], df_u["y"]), df_u["u"], (Xi, Yi), method="cubic")
Vi = griddata((df_v["x"], df_v["y"]), df_v["v"], (Xi, Yi), method="cubic")
speed = np.sqrt(Ui**2 + Vi**2)

# === Figure 1: Isolines ===
fig1, ax1 = plt.subplots(figsize=(6, 10))
levels = np.linspace(Z.min(), Z.max(), 20)
cs = ax1.contour(X, Y, Z, levels=levels, colors='black', linewidths=0.8, linestyles='solid')
ax1.set_aspect('equal')
ax1.axis('off')
plt.tight_layout()
plt.savefig(f"StreamFunction_Contours_Ra_{int(Ra)}.pdf", dpi=300, bbox_inches='tight')
plt.show()
plt.close()

# === Figure 2: Streamlines vectors ===
fig2, ax2 = plt.subplots(figsize=(6, 10))
stream = ax2.streamplot(Xi, Yi, Ui, Vi, color=speed, cmap='viridis', density=1.5, linewidth=1, arrowsize=1)
ax2.set_aspect('equal')
ax2.axis('off')
plt.tight_layout()
plt.savefig(f"StreamFunction_Streamlines_Ra_{int(Ra)}.pdf", dpi=300, bbox_inches='tight')
plt.show()
plt.close()

# === Figure 3: Colormap ===
fig3, ax3 = plt.subplots(figsize=(6, 10))
levels = np.linspace(Z.min(), Z.max(), 100)
cf = ax3.contourf(X, Y, Z, levels=levels, cmap='Blues')
#cbar3 = fig3.colorbar(cf, ax=ax3, label=r"$\psi$")
ax3.set_aspect('equal')
ax3.axis('off')
plt.tight_layout()
plt.savefig(f"StreamFunction_Colormap_Ra_{int(Ra)}.pdf", dpi=300, bbox_inches='tight')
plt.show()
plt.close()



#######################################################################

################## VELOCITY PROFILES AT MID SECTIONS ##################

#######################################################################

Ra_list = [1e8, 6.4e8, 1.07e9]
Ra_labels = [r"$Ra = 8\times10^7$", r"$Ra = 6.4\times10^8$", r"$Ra = 1.07\times10^9$"]
colors = ['tab:blue', 'tab:orange', 'tab:green']
styles = ['-', '--', '-.']

# === PROFILE U(y) at x = L/2 ===
plt.figure(figsize=(6, 8))
for Ra, label, color, style in zip(Ra_list, Ra_labels, colors, styles):
    fname = f"UyL2_Ra_{int(Ra)}.txt"
    data = np.loadtxt(fname)
    data = data[(data[:, 0] != 0.05) & (data[:, 0] != 3.95)]
    data = data[np.argsort(data[:, 0])]

    y_vals = data[:, 0]
    u_vals = data[:, 1]

    interp_u = PchipInterpolator(y_vals, u_vals)
    y_smooth = np.linspace(y_vals.min(), y_vals.max(), 300)
    u_smooth = interp_u(y_smooth)

    plt.plot(u_smooth, y_smooth, linestyle=style, color=color, linewidth=2, label=label)

plt.grid(True, linestyle='--', alpha=0.5)
plt.xlabel('u(y) at x = L/2')
plt.ylabel('y')
plt.xlim(-1200, 1200)
plt.legend(loc='center right', fontsize=10, frameon=True)
plt.tight_layout()
plt.savefig("U_Profile_Comparison.pdf")
plt.show()


# === PROFILE V(x) at y = H/2 ===
plt.figure(figsize=(8, 6))
for Ra, label, color, style in zip(Ra_list, Ra_labels, colors, styles):
    fname = f"VxL2_Ra_{int(Ra)}.txt"
    data = np.loadtxt(fname)
    data = data[(data[:, 0] != 0.05) & (data[:, 0] != 0.95)]
    data = data[np.argsort(data[:, 0])]

    x_vals = data[:, 0]
    v_vals = data[:, 1]

    interp_v = PchipInterpolator(x_vals, v_vals)
    x_smooth = np.linspace(x_vals.min(), x_vals.max(), 300)
    v_smooth = interp_v(x_smooth)

    plt.plot(x_smooth, v_smooth, linestyle=style, color=color, linewidth=2, label=label)

plt.grid(True, linestyle='--', alpha=0.5)
plt.xlabel('x')
plt.ylabel('v(x) at y = H/2')
plt.legend()
plt.tight_layout()
plt.savefig("V_Profile_Comparison.pdf")
plt.show()



###########################################################

################## KINETIC ENERGY BUDGET ##################

###########################################################

# Load data: columns = t, Ek, R_diff, R_conv, R_press, R_bous, R_total, dEk_num, error
data = np.loadtxt(f"EnergyTerms_Ra_{int(Ra)}.txt")

# Extract individual arrays
t = data[:, 0]             # Time
Ek = data[:, 1]            # Total kinetic energy
R_diff = data[:, 2]        # Viscous diffusion term
R_conv = data[:, 3]        # Convective transport term
R_press = data[:, 4]       # Pressure work term
R_bous = data[:, 5]        # Buoyancy forcing
R_total = data[:, 6]       # Sum of all budget terms
dEk_num = data[:, 7]       # Numerical derivative of Ek
error = data[:, 8]         # Absolute error

# === Plot 1: Kinetic Energy vs Time ===
plt.figure(figsize=(10, 6))
plt.plot(t, Ek, label='Kinetic Energy $E_k$', color='navy')
plt.xlabel('Time')
plt.ylabel('Kinetic Energy')
#plt.title('Time Evolution of Kinetic Energy')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'kinetic_energy_Ra_{int(Ra)}.pdf')
plt.show()

# === Plot 2: Buoyancy and Viscous Diffusion ===
plt.figure(figsize=(10, 6))
plt.plot(t, R_bous, label='Buoyancy Forcing', linestyle='-', color='purple')
plt.plot(t, R_diff, label='Viscous Diffusion', linestyle='--', color='red')
plt.xlabel('Time')
plt.ylabel('Energy Contribution')
#plt.title('Kinetic Energy Budget: Buoyancy and Viscous Diffusion')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'budget_buoyancy_viscous_Ra_{int(Ra)}.pdf')
plt.show()

# === Plot 3: Convective and Pressure Terms ===
plt.figure(figsize=(10, 6))
plt.plot(t, R_conv, label='Convective Transport', linestyle='-.', color='green')
plt.plot(t, R_press, label='Pressure Work', linestyle=':', color='orange')
plt.xlabel('Time')
plt.ylabel('Energy Contribution')
#plt.title('Kinetic Energy Budget: Convective and Pressure Terms')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'budget_convective_pressure_Ra_{int(Ra)}.pdf')
plt.show()

# === Plot 4: Compare Kinetic Energy vs Total Budget Term ===
plt.figure(figsize=(10, 6))
plt.plot(t, Ek, label='Kinetic Energy $E_k$', color='blue')
plt.plot(t, R_total, label=r'Total Budget Term $\mathcal{R}_{\text{tot}}$', color='black', linestyle='--')
plt.xlabel('Time')
plt.ylabel('Value')
#plt.title('Kinetic Energy vs Total Budget Term')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'Ek_vs_Rtotal_Ra_{int(Ra)}.pdf')
plt.show()

# === Plot 4: dEk_num vs R_total ===
plt.figure(figsize=(10, 6))
plt.plot(t, dEk_num, label='Numerical $\\frac{dE_k}{dt}$', color='darkcyan')
plt.plot(t, R_total, label='Analytical $\\mathcal{R}_{\\text{tot}}$', linestyle='--', color='black')
plt.xlabel('Time')
plt.ylabel('Rate of Change')
#plt.title('Comparison: Numerical vs Analytical Derivative')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'dEk_vs_Rtotal_Ra_{int(Ra)}.pdf')
plt.show()

# === Plot 5: Error vs Time ===
plt.figure(figsize=(10, 6))
plt.plot(t, error, label='|dEk_num - R_total|', color='crimson')
plt.xlabel('Time')
plt.ylabel('Absolute Error')
#plt.title('Error Between Numerical and Analytical dEk/dt')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'error_vs_time_Ra_{int(Ra)}.pdf')
plt.show()



#####################################################

################## PROBES ANALYSIS ##################

#####################################################

# Load the probes history
data = np.loadtxt(f"ProbesHistory_Ra_{int(Ra)}.txt")

time = data[:, 0]

# Bottom left probe
T_bl = data[:, 1]
u_bl = data[:, 2]
v_bl = data[:, 3]

# Bottom right probe
T_br = data[:, 4]
u_br = data[:, 5]
v_br = data[:, 6]

# Top left probe
T_tl = data[:, 7]
u_tl = data[:, 8]
v_tl = data[:, 9]

# Top right probe
T_tr = data[:, 10]
u_tr = data[:, 11]
v_tr = data[:, 12]

def plot_probe(time, T, u, v, label, filename):
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    fig.suptitle(f"{label} – Temporal Evolution", fontsize=20)

    axes[0].plot(time, T, color='tab:red', linestyle='-',
                 label='Temperature')
    axes[0].set_ylabel('T')
    axes[0].grid(True)
    axes[0].legend(fontsize=10, loc='upper right')

    axes[1].plot(time, u, color='tab:blue', linestyle='--',
                 label='u-Velocity')
    axes[1].set_ylabel('u')
    axes[1].grid(True)
    axes[1].legend(fontsize=10, loc='upper right')

    axes[2].plot(time, v, color='tab:green', linestyle='-.',
                 label='v-Velocity')
    axes[2].set_ylabel('v')
    axes[2].set_xlabel('Time [s]')
    axes[2].grid(True)
    axes[2].legend(fontsize=10, loc='upper right')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

plot_probe(time, T_bl, u_bl, v_bl, "Bottom Left Probe", f"Probe_BL_Ra_{int(Ra)}.pdf")
plot_probe(time, T_br, u_br, v_br, "Bottom Right Probe", f"Probe_BR_Ra_{int(Ra)}.pdf")
plot_probe(time, T_tl, u_tl, v_tl, "Top Left Probe", f"Probe_TL_Ra_{int(Ra)}.pdf")
plot_probe(time, T_tr, u_tr, v_tr, "Top Right Probe", f"Probe_TR_Ra_{int(Ra)}.pdf")
