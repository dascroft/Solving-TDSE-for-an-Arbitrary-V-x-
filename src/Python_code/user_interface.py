from sympy import var
from sympy import sympify
import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin, tan, pi, exp


x=np.linspace(-2*pi,2*pi,1000)
default_potential = x**2


def user_input():
    input_yn = input("Do you want to input a potential? (yes or no)  ")
    if input_yn.lower()=="y" or input_yn.lower()=="yes":
        file_or_command_line = input("Do you have a file with your potential? (yes/no) ")
        if file_or_command_line.lower() == "y" or file_or_command_line.lower() == "yes":
            user_file = input("What is your file path? ")
            readfile = open(user_file, 'r')
            file = readfile.read()
            from numpy import cos
            func = lambda x: eval(file)
            v = func(x)
            plt.plot (x, v)#place holder
            '''
            Put solving_tdse.py here to solve and plot the solutions
            '''
        elif file_or_command_line.lower() == "n" or file_or_command_line.lower() =="no":
            print("Follow command line potential generator")
        else:
            print("Not a valid input please either input yes or no ")
            user_input()
    elif input_yn.lower()=="n" or input_yn.lower()=="no":
        v = default_potential 
        plt.plot (x, v)#place holder
        print("Default potential V(x) = x^2")
    else:
        print("Not a valid input please either input yes or no ")
        user_input()
            
        #else:
            #print("FOLLOW COMMAND LINE MANUAL INPUTS")
        
user_input()
    

