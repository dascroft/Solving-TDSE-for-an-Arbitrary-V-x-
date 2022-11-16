# Solving the TDSE for a Arbitrary Potential V

Using tdse_solver.py, for a specified potention V(x), it is then solved outputting the solution to the Hamiltonian as a gif, saved to a desired path. 
The user can specify what function of V(x) they desire and specify the associated parameters, i.e $kx^2$, with parameters $k$ being the coefficent and $p$ being the power. Using the Crank–Nicolson method, the wavefunction is then resolved within the specified potential and the default output is given as an animated 'output.gif'.

### Intructions for using the program:
Open a program capable of running `.py` files (to use command line, use Anaconda Powershell Prompt for best results.)
* For manual access: 
  * Open file `src\Python_code\solving_tdse.py` in the program of choice. 
  * To input values, scroll to the bottom and replace the "`args.XXX`" terms in  with the chosen values - use the contained "help" instructions above to determine which are which.
* For command line access:
  * From the program root folder, use `cd src\Python_code` to set the run location to the same folder as `solving_tdse.py`.
  * Use command `python solving_tdse.py` to run the file from command line - the program should create a new `.gif` file in `src\Python_code`.
  * To input new values, run `python solving_tdse.py --help`, to have Powershell list all possible input variables. From there, input the desired variables in standard command line format. For example: 
```python
python solving_tdse.py --vx "2*x**2"
```
The example above generates a $V(x) = 2x^2$ potential and the TDSE is solved within this potential, The ouputted solution can be represented as either a .gif or heat map plot directed by the user.

### Contents:

* `.gitignore` contains files that should be ignored by git
* `LICENSE` the project license telling users who install your package the terms under which they can use your package
* `README.md` A [markdown](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) document telling users about the project
* `requirements.txt` contains the requirements for the project, you can install these with `pip install -r requirements.txt`
* `docs/` contains the documentation - we won't discuss this further here.
* `.github/workflows/python_test.yml` contains a [YAML](https://yaml.org/) file which determined how github Action are run
* `src/` contains the python code required to run the project
* `tests/` contains the tests of the python package(unfinished)



