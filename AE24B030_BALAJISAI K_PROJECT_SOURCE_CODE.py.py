import numpy as np
import matplotlib.pyplot as plt

n = int(input('''Enter 1 for benchmark simulation of circular firefront
Enter 2 for merger simulation of 1-circular and 1-line firefront
Enter 3 for merger simulation of 2-circular and 1-line firefront
Enter 4 for merger simulation of 3-circular '''))
#Domain Boundaries
L_x = 100
L_y = 100

#Resolution
N_x = 40
N_y = 40

#grid spacing based on resolution
del_x = L_x/(N_x-1)
del_y = L_y/(N_y-1)

#Grid space
x = np.linspace(-L_x/2, L_x/2, N_x)
y = np.linspace(0, L_y, N_y)

if n ==1:
    phi=np.zeros((N_x,N_y))
    x0=0
    y0=50
    for i in range(N_x):
        for j in range(N_y):
            dist_squared = (x0-x[i])**2+(y0-y[j])**2 + 1e-10
            dist = np.sqrt(dist_squared)-10
            phi[i][j] = dist
elif n==2:
    phi=np.zeros((N_x,N_y))
    x0=0
    y0=45
    for i in range(N_x):
        for j in range(N_y):
            dist_squared = (x0-x[i])**2+(y0-y[j])**2 + 1e-10
            dist = np.sqrt(dist_squared)-10
            
            phi_line_dist = y[j]-25
            phi[i][j] = np.minimum(phi_line_dist, dist)
elif n==4:
    phi=np.zeros((N_x,N_y))
    x0=-25
    y0=0
    x1 = 0
    y1=0
    x2 = 25
    y2=0
    for i in range(N_x):
        for j in range(N_y):
            dist_squared_1 = (x0-x[i])**2+(y0-y[j])**2 + 1e-10
            dist_1 = np.sqrt(dist_squared_1)-10

            dist_squared_2 = (x1-x[i])**2+(y1-y[j])**2 + 1e-10
            dist_2 = np.sqrt(dist_squared_2)-10

            dist_squared_3 = (x2-x[i])**2+(y2-y[j])**2 + 1e-10
            dist_3 = np.sqrt(dist_squared_3)-10

            phi[i][j] = np.minimum(np.minimum(dist_1,dist_2),dist_3)

elif n==3:
    phi=np.zeros((N_x,N_y))
    x0=-20
    y0=50
    x2 = 20
    y2=50
    for i in range(N_x):
        for j in range(N_y):
            dist_squared_1 = (x0-x[i])**2+(y0-y[j])**2 + 1e-10
            dist_1 = np.sqrt(dist_squared_1)-20

            phi_line_dist = y[j]-25

            dist_squared_3 = (x2-x[i])**2+(y2-y[j])**2 + 1e-10
            dist_3 = np.sqrt(dist_squared_3)-20

            phi[i][j] = np.minimum(np.minimum(dist_1,phi_line_dist),dist_3)

r0 = 0.165 
cf = 3.24
V_x = 0
V_y = 3
phi_at_times = []
def F(phi,n):
    U_x = np.zeros((N_x, N_y))
    U_y = np.zeros((N_x, N_y))
    F = np.zeros((N_x,N_y))
    for i in range(1,N_x-1):
        for j in range(1,N_y-1):
            x_der = (phi[i+1][j]-phi[i-1][j])/(2*del_x)
            y_der = (phi[i][j+1]-phi[i][j-1])/(2*del_y)
            mag = np.sqrt(x_der**2+y_der**2)
            n_x = x_der/mag
            n_y = y_der/mag
            if n==2 or n==3 or n==4:
                V_dot_n = V_x*n_x+V_y*n_y
                
                U_x[i,j] = r0*(1+cf*V_dot_n)*n_x
                U_y[i,j] = r0*(1+cf*V_dot_n)*n_y
            elif n==1:
                Vmag = np.sqrt(V_x**2 + V_y**2)

                cos_theta = (V_x*n_x + V_y*n_y)/(Vmag + 1e-10)

                theta = np.arccos(np.clip(cos_theta, -1, 1))

                if abs(theta) < np.pi/2:
                    spread = r0*(1 + cf*np.sqrt(Vmag*(cos_theta**1.5)))
                else:
                    a = 0.5
                    spread = r0*(a + (1-a)*np.sin(theta))

                U_x[i,j] = spread*n_x
                U_y[i,j] = spread*n_y

            if U_x[i,j]>0:
                lim_x_der = (phi[i][j]-phi[i-1][j])/del_x
            elif U_x[i,j]<=0:
                lim_x_der = (phi[i+1][j]-phi[i][j])/del_x  

            if U_y[i,j]>0:
                lim_y_der = (phi[i][j]-phi[i][j-1])/del_y
            elif U_y[i,j]<=0:
                lim_y_der = (phi[i][j+1]-phi[i][j])/del_y 

            F_phi = U_x[i,j]*lim_x_der+U_y[i,j]*lim_y_der
            F[i][j] = F_phi
    return F
def bcs(phi):
    for i in [0]:
        for j in range(N_y):
            phi[i][j] = phi[i+1][j]
    for i in [N_x-2]:
        for j in range(N_y):
            phi[i+1][j] = phi[i][j]
    for j in [0]:
        for i in range(N_x):
            phi[i][j] = phi[i][j+1]
    for j in [N_y-2]:
        for i in range(N_x):
            phi[i][j+1] = phi[i][j]
    return phi

time_step =0.1
total_time = 30
t = 0

phi = bcs(phi)
phi_at_times.append(np.copy(phi))

while t<total_time:
    F_ini = F(phi,n)

    phi_inter = phi-time_step*F_ini
    phi_inter = bcs(phi_inter)
    F_inter = F(phi_inter,n)

    phi_next = 0.5*phi+0.5*(phi_inter-time_step*F_inter)
    phi_next = bcs(phi_next)
    phi = phi_next
    t +=time_step
    if round(t,1)%6==0:
        phi_at_times.append(phi)

print("Simulation complete! Generating plots...")

if n == 1:
    # --- BENCHMARK PLOT 1: All time steps overlaid ---
    plt.figure(figsize=(8, 8))
    for idx, phi_saved in enumerate(phi_at_times):
        plt.contour(x, y, phi_saved.T, levels=[0], colors='red', linewidths=1.5)
    
    plt.title("Benchmark Simulation: Overlaid Progression (0s to 30s)")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.xlim(-50, 50)
    plt.ylim(0, 100)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.gca().set_aspect('equal')
    plt.show()

    # --- BENCHMARK PLOT 2: 6 separate subplots ---
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("Benchmark Simulation: Individual Time Steps", fontsize=16)
    
    for idx, ax in enumerate(axes.flat):
        if idx < len(phi_at_times):
            current_time = idx * 6
            ax.contour(x, y, phi_at_times[idx].T, levels=[0], colors='red', linewidths=1.5)
            ax.set_title(f"Time = {current_time} s")
            ax.set_xlabel("X (m)")
            ax.set_ylabel("Y (m)")
            ax.set_xlim(-50, 50)
            ax.set_ylim(0, 100)
            ax.grid(True, linestyle='--', alpha=0.5)
            ax.set_aspect('equal')
            
    # The fix: Apply tight_layout, then manually adjust the spacing
    plt.tight_layout()
    fig.subplots_adjust(top=0.9, hspace=0.4, wspace=0.3)
    plt.show()

elif n == 2:
    # --- MERGER PLOT: 6 separate subplots ---
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("Merger Simulation: Line and Spot Fire", fontsize=16)
    
    for idx, ax in enumerate(axes.flat):
        if idx < len(phi_at_times):
            current_time = idx * 6
            ax.contour(x, y, phi_at_times[idx].T, levels=[0], colors='red', linewidths=1.5)
            ax.set_title(f"Time = {current_time} s")
            ax.set_xlabel("X (m)")
            ax.set_ylabel("Y (m)")
            ax.set_xlim(-50, 50)
            ax.set_ylim(0, 100)
            ax.grid(True, linestyle='--', alpha=0.5)
            ax.set_aspect('equal')
            
    # The fix: Apply tight_layout, then manually adjust the spacing
    plt.tight_layout()
    fig.subplots_adjust(top=0.9, hspace=0.4, wspace=0.3)
    plt.show()

elif n == 4:
    # --- THREE-FRONT MERGER PLOT: All time steps overlaid ---
    plt.figure(figsize=(8, 8))
    
    # Loop through every saved time step and plot its zero-contour
    for idx, phi_saved in enumerate(phi_at_times):
        plt.contour(x, y, phi_saved.T, levels=[0], colors='red', linewidths=1.5)
    
    plt.title("Three-Front Merger Simulation (0s to 30s)")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.xlim(-50, 50)
    plt.ylim(0, 100)
    plt.grid(True, linestyle='--', alpha=0.5)
    
    # Forces the X and Y axes to have the same scale
    plt.gca().set_aspect('equal') 
    
    plt.show()

elif n == 3:
    # --- UNBURNED POCKET PLOT: 4-panel grid (2x2) ---
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    fig.suptitle("Unburned Fuel Pocket Simulation", fontsize=16)
    
    # We select 4 specific time steps to mimic the 4 frames in the paper.
    # Assuming your list saved t = 0, 6, 12, 18, 24, 30 (indices 0 through 5)
    indices_to_plot = [0, 1, 2,3] # Ploting t=0s, 12s, 24s, and 30s
    
    for ax_idx, list_idx in enumerate(indices_to_plot):
        # Make sure the index actually exists in your saved list before plotting
        if list_idx < len(phi_at_times):
            ax = axes.flat[ax_idx]
            current_time = list_idx * 6
            
            ax.contour(x, y, phi_at_times[list_idx].T, levels=[0], colors='red', linewidths=1.5)
            ax.set_title(f"Time = {current_time} s")
            ax.set_xlabel("X (m)")
            ax.set_ylabel("Y (m)")
            ax.set_xlim(-50, 50)
            ax.set_ylim(0, 100)
            ax.grid(True, linestyle='--', alpha=0.5)
            ax.set_aspect('equal')
            
    plt.tight_layout()
    fig.subplots_adjust(top=0.92, hspace=0.3, wspace=0.3)
    plt.show()
