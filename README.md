# PathtrackingTechniques
A pathtracking techinique implemented, which I had to research and use for my dissertation. The method used involves creating a homotopy between the images of a polynomial system, f(x1,...,xn)=(f1,f2,...,fn), whose roots (solutions to f=0) we want to find and a similar polynomial system, g, with the same number of roots as f, whose roots we already know. I refer you to the text by Daniel J. Bates et al "Numerically Solving Polynomial Systems with Bertini" for a complete mathematical description of the methods implemented here.

- newton.py: Simple definitions of Newton's method applied to single variable equations to be later used.

- single_variable_pathtracking.py: Code for and example of single variable pathtracking, solutions verified via newton's method.

- multivariable_pathtracking.py: Code for, and example of, pathtracking for a system of 2 variables and 2 equations (prerequisite for pathtracking is no. of variables=system dimesnion).

- pathtracking3.pdf: A pdf publish of multivariable_pathtracking.py in a jupyter notebook.

- multivariable_pathtracking_with_symbolic_tools: Same as multivariable_pathtracking.py, except it uses symbolic tools to allow for a easier alternation of the system of equations.
