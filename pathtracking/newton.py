def newton(g, dg, x=1):
    '''Returns simple newton iteration for the roots of g(x)'''
    # g: function
    # dg: derivative of g
    return x - g(x)/dg(x)


def iterate_newton(g, dg, x_0=1, error=10**(-3)):
    '''Iterates the newton method to approximate one of the roots of g(x) to a given error'''
    # g: function
    # dg: derivative of g
    # x_0: intial 'guess' of root
    # error: minimum error to return approximation to
    e = 1
    x = x_0
    while e > error:
        y = newton(g, dg, x)
        e = abs(y-x)
        x = y
    return x
