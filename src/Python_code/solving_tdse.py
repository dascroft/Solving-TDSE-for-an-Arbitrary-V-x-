import numpy as np
import matplotlib.pyplot as plt
from findiff import FinDiff
from scipy.sparse.linalg import inv
from scipy.sparse import eye, diags
import matplotlib.animation as animation


class TDSE(object):
    def __init__(self, **kwargs):
        
        self.Nx = 500
        self.xmin = float(kwargs.get("xmin",'-5'))
        self.xmax = float(kwargs.get("xmax",'5'))
        
        self.Nt = int(kwargs.get("Nt", '250'))
        self.tmin = float(kwargs.get("tmin",'0')) 
        self.tmax = float(kwargs.get("tmax",'20'))
        
        self.k = float(kwargs.get("k",'1'))
        
        self.x_array = np.linspace(self.xmin, self.xmax, self.Nx)
        self.t_array = np.linspace(self.tmin, self.tmax, self.Nt)
        
        self.vx = kwargs.get("Form",'self.k * self.x_array ** 2')
        self.psi = np.exp(-(self.x_array+2)**2)
        
        self.left_wall_pstn = float(kwargs.get("Left_wall_position", '-4'))
        self.right_wall_pstn = float(kwargs.get("Right_wall_position", '4'))
        self.wall_height = kwargs.get("Wall_height", '1')
        
        self.barrier_width = float(kwargs.get("Barrier_width", '1'))
        self.barrier_pstn = kwargs.get("Barrier_position", 'none')   
        self.barrier_height = kwargs.get("Barrier_Height", '1')

        self.tracker = 0
              
        if 'square' in self.vx:
            TDSE.square(self)
   
        
    def square(self):
        x = np.linspace(self.xmin, self.xmax, 500)
        self.vx = np.zeros(len(x))

        self.vx[x<self.left_wall_pstn] = self.wall_height
        self.vx[x>self.right_wall_pstn] = self.wall_height
        
        if self.barrier_pstn != "none":
            self.barrier_pstn = float(self.barrier_pstn)  
            barrier_left = self.barrier_pstn - 0.5*self.barrier_width
            barrier_right = self.barrier_pstn + 0.5*self.barrier_width           
            self.v_x[(barrier_left<x) & (x<barrier_right)] = self.barrier_height
            
        self.tracker = 1


    def solve(self, x_array, t_array):
  
        # Calculate finite difference elements
        dt = self.t_array[1] - self.t_array[0]
        dx = self.x_array[1] - self.x_array[0]
        
        # Convert to a diagonal matrix
        if self.tracker == 1:
            vx_matrix = diags(self.vx)
        else:
            vx_matrix = diags(eval(self.vx))

        # Calculate the Hamiltonian matrix
        H = -0.5 * FinDiff(0, dx, 2).matrix(x_array.shape) + vx_matrix

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
        return self.psi_list

    def run(self, psi):
        self.line.set_data(self.x_array, np.abs(psi)**2)
        return self.line,
    
    def plot(self):
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set(xlabel="x [arb units]", ylabel="time [arb units]")
        
        psi_mag_squared = np.abs(TDSE.solve(self.x_array, self.t_array))**2  
        c = ax.pcolor(self.x_array, self.t_array, psi_mag_squared, shading="auto")
        
        cbar = fig.colorbar(c, ax=ax)
        cbar.set_label("$|\Psi(x, t)|^2$")
        plt.show()
        
    def animate(self):
        fig, ax = plt.subplots()       
        plt.rcParams["axes.labelsize"] = 16
        
        ax.set_xlabel("x [arb units]")
        ax.set_ylabel("$|\Psi(x, t)|$", color="C0")
              
        ax_twin = ax.twinx()
        if self.tracker == 1:
            ax_twin.plot(self.x_array, self.vx, color="C1")
        else:
            ax_twin.plot(self.x_array, eval(self.vx), color="C1")
            
        ax_twin.set_ylabel("V(x) [arb units]", color="C1")
   
        self.line, = ax.plot([], [], color="C0", lw=2)
        ax.grid()

        ax.set_xlim(self.x_array[0], self.x_array[-1])
        ax.set_ylim(0, 1)
        
        ani = animation.FuncAnimation(fig, self.run, TDSE.solve(self.x_array, self.t_array), interval=10)
        ani.save("particle_in_a_well.gif", fps=120, dpi=300)


TDSE = TDSE()
#TDSE.animate()
TDSE.plot()