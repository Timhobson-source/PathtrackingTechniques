import pathtracking.newton as newt

# EXAMPLE: Path tracking algorithm for a single polynomial equation.

# Try looking to find the solutions to the equation: f(z) = 5z^2 +z-37 = 0.
# From the discriminant 1^2 - 4*5*(-37) > 0, we see that there are two real solutions

# Newton's method solutions

# g(z)=z^2 -1 = 0 is an equation we know the solutions to, and has the same number of solutions
# to f(z)=0:


def g(z):
    return z**2-1  # roots 1,-1.


def dg(z):
    return 2*z

# f(z)=0 equation whose solutions we don't know:


def f(z):
    return 5*z**2 + z - 37


def df(z):
    return 10*z + 1

# calculating the solutions to f=0 to error 10e-5 by newton's method:


fsol_1 = newt.iterate_newton(f, df, x_0=5, error=0.0001)
fsol_2 = newt.iterate_newton(f, df, x_0=13, error=0.0001)

#####################################################################################
# Pathtracking solutions
####################################################################################

# Need to define our Hessian, H, as:

gamma = 0.44+0.77j  # prevents paths between the roots crossing.


def H(z, t):
    return gamma*t*g(z)+(1-t)*f(z)


def JH(z, t):  # jacobian of H with respect to z i.e. dH/dz in single variable case
    return gamma*t*dg(z)+(1-t)*df(z)


def Dt_H(z, t):
    return gamma*g(z)-f(z)  # dH/dt

# Have to solve IVP: p'(t) = -(JH(p(t),t)**(-1))*Dt_H(p(t),t)=F(p(t),t), p(1)= 1,
# where 1 is a root of gamma. We want p(0).


def euler(x, t, step, func):
    ''' Single explicit euler step, with the step in the negative t-direction'''
    # x: initial value
    # t: initial time
    # step>0: step in neg t-direction
    return x-step*func(x, t)


def F(x, t):
    return -(JH(x, t)**(-1))*Dt_H(x, t)


def LAonestep(p, t, step):
    '''
    Returns a first order prediction of p(t-step), it uses explicit euler
    method and two newton predictive corrections.
    '''
    z_0 = euler(p, t, step, F)
    z_1 = z_0-H(z_0, t)/JH(z_0, t+step)
    p_1 = z_1-H(z_1, t+step)/JH(z_1, t+2*step)
    return p_1


# solving the equation:
p = 1  # root1 of g(z)=0
q = -1  # root2 of g(z)=0
t = 1.0  # intial 'time'
N = 500  # number of iterations of LAonestep
for i in range(1, N+1):
    r = LAonestep(q, t, -1/N)
    s = LAonestep(p, t, -1/N)
    p = s
    q = r
    t -= 1/N
print('Pathtracking solutions: ', p, q)
print("Newton's method solutions (for checking)", fsol_1, fsol_2)

# It works!

if __name__ == '__main__':
    pass
