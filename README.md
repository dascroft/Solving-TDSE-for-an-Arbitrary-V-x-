# Solving the TDSE for a Arbitrary Potential V

Using tdse_solver.py, for a specified potention V(x) is then solved outputting the solution to the Hamiltonian as a gif, saved to a desired path. 
The user can specify what function of V(x) they desire and specify the associated parameters, i.e $kx^2$, with parameters $k$ being the coefficent and $p$ being the power. Using the Crank–Nicolson method, the wavefunction is then resolved within the specified potential and the output is given as animated '<output>.gif'.

### Intructions for using the program:





* `.gitignore` contains files that should be ignored by git
* `LICENSE` the project license telling users who install your package the terms under which they can use your package
* `README.md` A [markdown](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) document telling users about the project
* `pyproject.toml` tells build tools (like pip and build) what is required to build your project.
* `requirements.txt` contains the requirements for the project, you can install these with `pip install -r requirements.txt`
* `docs/` contains the documentation - we won't discuss this further here.
* `.github/workflows/python_test.yml` contains a [YAML](https://yaml.org/) file which determined how github Action are run
* `src/` contains the python package itself
* `tests/` contains the tests of the python package



