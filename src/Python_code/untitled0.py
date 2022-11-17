import argparse
import solving_tdse as tdse

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

    return args

if __name__ == "__main__":
    args = argue()
    kwargs_list = ['vx = ', 'xmin = ', 'xmax = ', 'tmin = ', 'tmax = ', 'k = ', 'p = ', 'output =']            
    tdse.TDSE(**vars(args))

