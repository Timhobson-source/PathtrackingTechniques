"""Polynomial system path tracking with symbolic tools usage.

You may easily edit this code to your own 2 variable dim=2 polynomial system of equations
with two solutions.
This can be done by editing 'f1' and 'f2' in the definiton of f(x,y) below.

Generalization of this code to solve higher dimensional systems is possible but not as simple.

"""
import numpy as np
import sympy as sp
import math


def f(x, y):
    # system of equations we want to solve (can be of higher dimension)
    return np.array([x**2 - y, x * y + 1])  # f(x,y)=[f1,f2]


# X = np.linspace(-10, 10, num=100)
# Y = np.linspace(-10, 10, num=100)
# u, v = np.meshgrid(X, Y)
# f1 = f(u, v)[0]
# f2 = f(u, v)[1]

# # plotting surfaces z=f1, z=f2
# fig = plt.figure(figsize=(8, 8))
# ax1 = fig.add_subplot(211, projection='3d')
# ax2 = fig.add_subplot(212, projection='3d')
# surf1 = ax1.plot_surface(u, v, f1, color="r")
# surf2 = ax2.plot_surface(u, v, f2, color="b")
# ax1.set_xlabel("$x$")
# ax2.set_xlabel("$x$")
# ax1.set_ylabel("$y$")
# ax2.set_ylabel("$y$")
# ax1.set_zlabel("$f 1$")
# ax2.set_zlabel("$f 2$")
# plt.show()


# Path tracking section

gamma = 0.44+0.77j


def g(x, y):
    # g(x,y)=0 has known solutions: (1,1),(1,-1),(-1,1) and (-1,-1).
    return np.array([x**2-1, y**2-1])


def H(x, y, t):
    # Homotopy between f and g: 0<=t<=1
    return gamma*t*g(x, y)+(1-t)*f(x, y)


x, y, t = sp.symbols('x y t')

mat_f = sp.Matrix(f(x, y))
Jacobf = mat_f.jacobian([x, y])
Jf = sp.lambdify((x, y), Jacobf)  # Jacobian of f wrt (x,y)

mat_g = sp.Matrix(g(x, y))
Jacobg = mat_g.jacobian([x, y])
Jg = sp.lambdify((x, y), Jacobg)  # jacobian of g wrt (x,y)

mat_H = sp.Matrix(H(x, y, t))
JacobH = mat_H.jacobian([x, y])
JH = sp.lambdify((x, y, t), JacobH)  # Jacobian of H wrt (x,y)

Dt_H = sp.lambdify([x, y, t], sp.diff(H(x, y, t), t))


def euler2D(x0, y0, t0, step, func):
    # 2D explicit euler iteration
    # x,y: 2D initial spacial values
    # t: intial time
    # func: function returning 2D array satisfying d(x,y)/dt=func(x,y)
    return np.array([x0, y0])-step*func(x0, y0, t0)


def F(x, y, t):
    return -np.dot(np.linalg.inv(JH(x, y, t)), Dt_H(x, y, t))


def iterate2D(x0, y0, t0, dt):
    # single path tracking iteration
    # path tracking: t=1 -> t=0 (dt<0)
    z0 = euler2D(x0, y0, t0, dt, F)
    z1 = z0-np.dot(np.linalg.inv(JH(z0[0], z0[1], t0+dt)), H(z0[0], z0[1], t0))
    z2 = z1 - \
        np.dot(np.linalg.inv(
            JH(z1[0], z1[1], t0+2*dt)), H(z1[0], z1[1], t0+dt))
    return z2


if __name__ == '__main__':
    # Using above functions to implement pathtracking algorithm and find the solutions:

    X1 = complex(1, 0)
    Y1 = complex(1, 0)
    X2 = complex(1, 0)
    Y2 = complex(-1, 0)
    X3 = complex(-1, 0)
    Y3 = complex(-1, 0)

    T = 1
    N = 1000
    for i in range(0, N):
        s1 = iterate2D(X1, Y1, T, dt=-1/N)
        s2 = iterate2D(X2, Y2, T, dt=-1/N)
        s3 = iterate2D(X3, Y3, T, dt=-1/N)
        X1 = s1[0]
        Y1 = s1[1]
        X2 = s2[0]
        Y2 = s2[1]
        X3 = s3[0]
        Y3 = s3[1]
        T -= 1/N

    # approximated solutions with pathtracking:

    print(round(X1, 3), round(Y1, 3), "\t\t\tf(sol1): ",
          round(np.linalg.norm(f(X1, Y1))))  # sol1
    print(round(X2, 3), round(Y2, 3), "\t\t\tf(sol2): ",
          round(np.linalg.norm(f(X2, Y2))))  # sol2
    print(round(X3, 3), round(Y3, 3), "\t\t\tf(sol3): ",
          round(np.linalg.norm(f(X3, Y3))))  # sol3

    # exact solutions comparison from analytic methods (for the given example)
    sol1 = [complex(0.5, 0.5*math.sqrt(3)), complex(-0.5, 0.5*math.sqrt(3))]
    sol2 = [-1, 1]
    sol3 = [complex(0.5, -0.5*math.sqrt(3)), complex(-0.5, -0.5*math.sqrt(3))]
    print("sol1: ", sol1)
    print("\nsol2: ", sol2)
    print("\nsol3: ", sol3)
    # As we can see, it works! And has generated all the solutions.
