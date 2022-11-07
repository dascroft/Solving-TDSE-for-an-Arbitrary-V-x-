import numpy as np
import matplotlib.pyplot as plt
from findiff import FinDiff
from scipy.sparse.linalg import inv
from scipy.sparse import eye, diags
import matplotlib.animation as animation

plt.rcParams["axes.labelsize"] = 16

# Input parameters
Nx = 500
xmin = -5
xmax = 5

Nt = 250
tmin = 0
tmax = 20
k = 1

fig, ax = plt.subplots()
line, = ax.plot([], [], color="C0", lw=2)
x_array = np.linspace(xmin, xmax, Nx)
t_array = np.linspace(tmin, tmax, Nt)
v_x = k * x_array ** 2
psi = np.exp(-(x_array+2)**2)

class TDSE(object):
    def __init__(self):
        
        self.Nx = 500
        self.xmin = -5
        self.xmax = 5
        self.Nt = 250
        self.tmin = 0
        self.tmax = 20
        self.v_x = k * x_array ** 2
        self.psi = np.exp(-(x_array+2)**2)
        self.k = 1
        self.x_array = x_array
    def solve(self, x_array, t_array):
  
        # Calculate finite difference elements
        dt = t_array[1] - t_array[0]
        dx = x_array[1] - x_array[0]
        
        # Convert to a diagonal matrix
        v_x_matrix = diags(self.v_x)

        # Calculate the Hamiltonian matrix
        H = -0.5 * FinDiff(0, dx, 2).matrix(x_array.shape) + v_x_matrix

        # Apply boundary conditions to the Hamiltonian
        H[0, :] = H[-1, :] = 0
        H[0, 0] = H[-1, -1] = 1

        # Calculate U
        I_plus = eye(self.Nx) + 1j * dt / 2. * H
        I_minus = eye(self.Nx) - 1j * dt / 2. * H
        U = inv(I_minus).dot(I_plus)

        # Iterate over each time, appending each calculation of psi to a list
        self.psi_list = []
        for t in t_array:
            self.psi = U.dot(self.psi)
            self.psi[0] = self.psi[-1] = 0
            self.psi_list.append(np.abs(self.psi))

    def run(self, psi):
        line.set_data(self.x_array, np.abs(psi)**2)
        return line,
        
    def animate(self, x_array):
        
        ax.set_xlabel("x [arb units]")
        ax.set_ylabel("$|\Psi(x, t)|$", color="C0")

        ax_twin = ax.twinx()
        ax_twin.plot(x_array, self.v_x, color="C1")
        ax_twin.set_ylabel("V(x) [arb units]", color="C1")
   
        ax.grid()
        xdata, ydata = [], []

        ax.set_xlim(x_array[0], x_array[-1])
        ax.set_ylim(0, 1)
        ani = animation.FuncAnimation(fig, self.run, self.psi_list, interval=10)
        ani.save("particle_in_a_well.gif", fps=120, dpi=300)


TDSE = TDSE()
TDSE.solve(x_array, t_array)
TDSE.animate(x_array)



    
