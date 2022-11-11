
import numpy as np
import matplotlib.pyplot as plt
from findiff import FinDiff
from scipy.sparse.linalg import inv
from scipy.sparse import eye, diags
import matplotlib.animation as animation
import string
import argparse
from scipy import signal

# Variant functions for vx
def form1(x,k,p):
    return k*(x**p)

def form2(x,k,p):
    return k*np.exp(p*x)

def form3(x,k,p):
    return k*np.sin(p*x)

def form4(x,k,p):
    return k*np.cos(p*x)

def form5(x, k ,p):
    return signal.square(k * np.pi/p * x)

# Function selector (to be called)
def fsel(x,k,p,n):
    if n == 1:
        return form1(x,k,p)
    elif n == 2: 
        return form2(x,k,p)
    elif n == 3:
        return form3(x,k,p)
    elif n == 4:
        return form4(x,k,p)
    elif n == 5:
        return form5(x,k,p)
    else:
        raise ValueError(f"The input {number} does not correspond to a valid function, please input an integer from 1-5.")
        return

# Set X and T arrays
Nx = 500

Nt = 250
tmin = 0
tmax = 20

t_array = np.linspace(tmin, tmax, Nt)

def gifify(xmin,xmax,k,p,n):
    x_array = np.linspace(xmin, xmax, Nx)
    #Define VX and Psi
    v_x = fsel(x_array,k,p,n)
    psi = np.exp(-(x_array+2)**p)

    # Calculate finite difference elements
    dt = t_array[1] - t_array[0]
    dx = x_array[1] - x_array[0]

    # Convert to a diagonal matrix
    v_x_matrix = diags(v_x)

    # Calculate the Hamiltonian matrix
    H = -0.5 * FinDiff(0, dx, 2).matrix(x_array.shape) + v_x_matrix

    # Apply boundary conditions to the Hamiltonian
    H[0, :] = H[-1, :] = 0
    H[0, 0] = H[-1, -1] = 1

    # Calculate U
    I_plus = eye(Nx) + 1j * dt / 2. * H
    I_minus = eye(Nx) - 1j * dt / 2. * H
    U = inv(I_minus).dot(I_plus)

    # Iterate over each time, appending each calculation of psi to a list
    psi_list = []
    for t in t_array:
        psi = U.dot(psi)
        psi[0] = psi[-1] = 0
        psi_list.append(np.abs(psi))

    fig, ax = plt.subplots()

    ax.set_xlabel("x [arb units]")
    ax.set_ylabel("$|\Psi(x, t)|$", color="C0")

    ax_twin = ax.twinx()
    ax_twin.plot(x_array, v_x, color="C1")
    ax_twin.set_ylabel("V(x) [arb units]", color="C1")

    line, = ax.plot([], [], color="C0", lw=2)
    ax.grid()
    xdata, ydata = [], []
    
    def run(psi):
        line.set_data(x_array, np.abs(psi)**2)
        return line
    
    ax.set_xlim(x_array[0], x_array[-1])
    ax.set_ylim(0, 1)

    ani = animation.FuncAnimation(fig, run, psi_list, interval=10)
    ani.save("particle_in_a_well.gif", fps=120, dpi=300)

def argue():
    parser = argparse.ArgumentParser()
    parser.add_argument("--form", type=int, help="The form of the equation. Takes value 1-5.", required=True)
    parser.add_argument("--xmin", type=int, help="The lower bound of the potential. Takes integer.", required=True)
    parser.add_argument("--xmax", type=int, help="The upper bound of the potential. Takes integer greater than xmin.", required=True)
    parser.add_argument("--k", type=int, help="Coefficient #1.", required=True)
    parser.add_argument("--p", type=int, help="Coefficient #2.", required=True)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = argue()
    gifify(xmin = args.xmin,xmax = args.xmax,k = args.k, p = args.p, n = args.form)
