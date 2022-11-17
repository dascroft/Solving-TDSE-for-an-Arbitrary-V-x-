import numpy as np
import matplotlib.pyplot as plt
from findiff import FinDiff
from scipy.sparse.linalg import inv
from scipy.sparse import eye, diags
import matplotlib.animation as animation
import argparse
import warnings

class TDSE(object):
    '''Solves the time dependant schrodinger equation, in an arbitrary potential.
        Parameters, potential and output can be customised with key word arguments:
        Output = "" - determines form of the output, options are "anim" for a gif, or "plot" for a 3d plot
            
        xmin = "" - set minimum x value, default -5
        xmax = "" - set maximum x value, default 5
        Nx = ""   - set number of x values between xmin and xmax, default 500
        
        tmin = "" - set starting value for time, default 0
        tmax = "" - set the time the program solves up to, default 20
        Nt = ""   - set number of time values between tmin and tmax, default 250
        
        vx = ""   - define the potential to solve in, enter as python code using x as the variable, and "np." for numpy. Default is x**2
                    Note, to get a square well; enter 'square' for vx.  
        
        Used only for square well:
        Left_wall_position = ""  - sets the position of the left wall of the well, default -4
        Right_wall_position = "" - sets the position of the right wall of the well, default 4
        Wall_height = ""         - sets the hight of the well, default 1
        
        Barrier_position = ""    - sets the position of the center of the barrier, default none  (defaults to no barrier)
        Barrier_height = ""      - sets the height of the barrier (if present), default 1
        '''
    def __init__(self, **kwargs):
        try:
            self.xmin = kwargs.get("xmin",'-5')
            if self.xmin == None:
                self.xmin = -5.0
            self.xmin = float(self.xmin)
        except:
            warnings.warn(f'{kwargs.get("xmin")} is not valid as xmin, using default of -5\n') 
            self.xmin = float(-5)
            pass
        try:
            self.xmax = kwargs.get("xmax",'5')
            if self.xmax == None:
                self.xmax = 5.0
            self.xmax = float(self.xmax)
        except:
            warnings.warn(f'{kwargs.get("xmax")} is not valid as xmax, using default of 5\n') 
            self.xmax = float(5)
            pass
        try:
            self.Nx = int(kwargs.get("Nx",'500')) #number of x values between xmin & xmax
        except:
            warnings.warn(f'{kwargs.get("Nx")} is not valid as Nx, using default of 500\n') 
            self.Nx = int(500)
            pass

       #values determining the time range
        try:
            self.tmin = kwargs.get("tmin",'0')
            if self.tmin == None:
                self.tmin = 0
            self.tmin = float(self.tmin)
        except:
            warnings.warn(f'{kwargs.get("tmin")} is not valid as tmin, using default of 0\n') 
            self.tmin = float(0)

            pass        
        try:
            self.tmax = kwargs.get("tmax",'20')
            if self.tmax == None:
                self.tmax = 20
            self.tmax = float(self.tmax)
        except:
            warnings.warn(f'{kwargs.get("tmax")} is not valid as tmax, using default of 20\n') 
            self.tmax = float(20)
            pass
        try:
            self.Nt = int(kwargs.get("Nt",'250'))#number of t values between tmin & tmax
        except:
            warnings.warn(f'{kwargs.get("Nt")} is not valid as Nt, using default of 25\n') 
            self.Nt = int(250)
            pass

        try:
            self.k = kwargs.get("k",'1')
            if self.k == None:
                self.k = 1
            self.k = float(self.k)
        except:
            warnings.warn(f'{kwargs.get("k")} is not valid as k, using default of k\n') 
            self.k = float(1)
            pass  
        try:
            self.p = kwargs.get("p",'2')
            if self.p == None:
                self.p = 2
            self.p = float(self.p)
        except:
            warnings.warn(f'{kwargs.get("p")} is not valid as p, using default of 2\n') 
            self.p = float(2)
            pass       

        #generates arrays for x & t values
        self.x_array = np.linspace(self.xmin, self.xmax, self.Nx)
        self.t_array = np.linspace(self.tmin, self.tmax, self.Nt)
        

        self.vx = kwargs.get('vx',"self.k * self.x_array ** self.p")
        self.psi = np.exp(-(self.x_array+2)**2)
        
        if self.vx == None:
            self.vx = "self.k * self.x_array ** self.p"
        
        #properties of wall of square well        
        try:
            self.left_wall_pstn = float(kwargs.get("Left_wall_position",'-4'))
        except:
            warnings.warn(f'{kwargs.get("Left_wall_position")} is not valid as Left_wall_position, using default of -4\n') 
            self.left_wall_pstn = float(-4)
            pass
        try:
            self.right_wall_pstn = float(kwargs.get("Right_wall_position",'4'))
        except:
            warnings.warn(f'{kwargs.get("Right_wall_position")} is not valid as Right_wall_position, using default of 4\n') 
            self.right_wall_pstn = float(4)
            pass
        try:
            self.wall_height = float(kwargs.get("Wall_height",'1'))
        except:
            warnings.warn(f'{kwargs.get("Wall_height")} is not valid as Wall_height, using default of 1\n') 
            self.wall_height = float(1)
            pass
        
        #properties of optional central barrier for square well
        try:
            self.barrier_width = float(kwargs.get("Barrier_width",'1'))
        except:
            warnings.warn(f'{kwargs.get("Barrier_width")} is not valid as Barrier_width, using default of 1\n') 
            self.barrier_width = float(1)
            pass
        try:
            self.barrier_height = float(kwargs.get("Barrier_height",'1'))
        except:
            warnings.warn(f'{kwargs.get("Barrier_height")} is not valid as Barrier_height, using default of 1\n') 
            self.barrier_height = float(1)
            pass
        self.barrier_position = kwargs.get("Barrier_position", 'none')  #defaults to 'none', for no barrier

        self.tracker = 0  #used to track use of special cases 
        
        #checks if user has asked for square well, runs TDSE.square if so
        if self.vx != None and 'square' in self.vx:
            try:
                self.barrier_position = float(self.barrier_position)  
            except:
                warnings.warn(f'{self.barrier_position} is not valid as Barrier_position, using default of 0\n') 
                self.barrier_position = float(0)
                pass
        
        #output clause
        self.output = kwargs.get("output", 'anim')
        if self.output == None:
            self.output = "anim"

        
    def square(self):
        x = np.linspace(self.xmin, self.xmax, self.Nx)
        self.vx = np.zeros(len(x))  #initially sets vx to a flat line at 0 
        
        #set values in vx array to wall height for x below left position or above right position 
        self.vx[x<self.left_wall_pstn] = self.wall_height
        self.vx[x>self.right_wall_pstn] = self.wall_height
        
        if self.barrier_position != "none":  #check if user has added barrier by changing from default
            self.barrier_position = float(self.barrier_position)  
            
            #calculated x coords of left and right side of barrier
            barrier_left = self.barrier_position - 0.5*self.barrier_width
            barrier_right = self.barrier_position + 0.5*self.barrier_width  
            
            #set values of vx between the barriers left and right sides to the specified barrier height
            self.v_x[(barrier_left<x) & (x<barrier_right)] = self.barrier_height
            
        self.tracker = 1  #used to note a special case has been used
        
    def checks(self):
        if self.xmin >= self.xmax:
            raise ValueError('xmax is greater than/equal to xmin')
        
        if self.tmin >= self.tmax:
            raise ValueError('tmax is greater than/equal to tmin')
        
        if self.left_wall_pstn >= self.right_wall_pstn:
            raise ValueError('left side of barrier is greater than/equal to right side of barrier')
            
        if self.wall_height < 0:
            raise ValueError('wall height is negative')
        if self.tracker == 1:    
            if float(self.barrier_height) < 0:
                raise ValueError('barrier height is negative')

    def solve(self):
  
        # Calculate finite difference elements
        dt = self.t_array[1] - self.t_array[0]
        dx = self.x_array[1] - self.x_array[0]
        
        # Convert to a diagonal matrix, eval needs to not be used if the square well special case is used
        if self.tracker == 1:
            vx_matrix = diags(self.vx)
        else:
            vx_matrix = diags(eval(self.vx))

        # Calculate the Hamiltonian matrix
        H = -0.5 * FinDiff(0, dx, 2).matrix(self.x_array.shape) + vx_matrix

        # Apply boundary conditions to the Hamiltonian
        H[0, :] = H[-1, :] = 0
        H[0, 0] = H[-1, -1] = 1

        # Calculate U
        I_plus = eye(self.Nx) + 1j * dt / 2. * H
        I_minus = eye(self.Nx) - 1j * dt / 2. * H
        U = inv(I_minus).dot(I_plus)

        # Iterate over each time, appending each calculation of psi to a list
        self.psi_list = []
        for t in self.t_array:
            self.psi = U.dot(self.psi)
            self.psi[0] = self.psi[-1] = 0
            self.psi_list.append(np.abs(self.psi))
        return self.psi_list

    def run(self, psi):
        self.line.set_data(self.x_array, np.abs(psi)**2)
        return self.line,
    
    def plot(self):
        self.checks()
        
        #creates a blank plot and adds axis labels
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set(xlabel="x [arb units]", ylabel="time [arb units]")
        
        #calculate magnitude squared of psi
        psi_mag_squared = np.abs(self.solve())**2  
        c = ax.pcolor(self.x_array, self.t_array, psi_mag_squared, shading="auto") #create a heatmap plot of x,t,mag^2 of psi
        
        #adds scale for the colour
        cbar = fig.colorbar(c, ax=ax)
        cbar.set_label("$|\Psi(x, t)|^2$")
        
        plt.savefig("plot.png")
        
    def animate(self):
        TDSE.checks(self)
        
        #creates a blank plot
        fig, ax = plt.subplots()       
        plt.rcParams["axes.labelsize"] = 16
        
        #set x label, and y label for psi
        ax.set_xlabel("x [arb units]")
        ax.set_ylabel("$|\Psi(x, t)|$", color="C0")
        ax_twin = ax.twinx()
        
        #plot vx
        if self.tracker == 1:
            ax_twin.plot(self.x_array, self.vx, color="C1")
        else:
            ax_twin.plot(self.x_array, eval(self.vx), color="C1")
        
        #set label for vx
        ax_twin.set_ylabel("V(x) [arb units]", color="C1")
   
        self.line, = ax.plot([], [], color="C0", lw=2)
        ax.grid()

        ax.set_xlim(self.x_array[0], self.x_array[-1])
        ax.set_ylim(0, 1)
        
        #animate the time evolution of psi and save as gif
        ani = animation.FuncAnimation(fig, self.run, TDSE.solve(self), interval=10)
        ani.save("particle_in_a_well.gif", fps=120, dpi=300)


def argue():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vx", type=str, help="A custom function, for a more complex wavefunction. Takes the function as a string.")
    parser.add_argument("--xmin", type=int, help="The lower bound of the potential. Takes integer.")
    parser.add_argument("--xmax", type=int, help="The upper bound of the potential. Takes integer greater than xmin.", required=False)
    parser.add_argument("--tmin", type=int, help="The lower bound of the time. Takes integer.", required=False)
    parser.add_argument("--tmax", type=int, help="The upper bound of the time. Takes integer greater than tmin.", required=False)
    parser.add_argument("--k", type=int, help="Coefficient #1.", required=False)
    parser.add_argument("--p", type=int, help="Coefficient #2.", required=False)
    parser.add_argument("--output", type=str, help = "The format of the output: .GIF ('anim') or plot ('plot').")
    args = parser.parse_args()
    
    #replacing x with self.x_array
    if args.vx != None:
        args.vx = args.vx.replace("x","self.x_array")
    else:
        args.vx == args.vx
    return args

if __name__ == "__main__":
    args = argue()
    TDSE = TDSE(vx = args.vx, xmin = args.xmin,xmax = args.xmax,tmin = args.tmin,tmax=args.tmax,k = args.k, p = args.p, output = args.output)
    if TDSE.output == "anim":
        TDSE.animate()
    elif TDSE.output == "plot":
        TDSE.plot()
    else:
        warnings.warn(f'The input {args.output} is not a valid output option. Please enter either "anim" or "plot".\n')
    print("")
