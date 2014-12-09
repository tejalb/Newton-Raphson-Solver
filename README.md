Newton-Raphson-Solver
=====================

This package implements a Newton-Raphson solver. Here is a description of the included files:

newton.py:
Implements the class newton, which returns a new object to find the roots of
f(x) = 0 using Newton Raphson method. Optional arguments are

tol: tolerance, iterate until |f(x)| < tol
maxiter: maximum number of iterations

dx: step size to compute the approximate Jacobian

Df: The analytic Jacobian. If not provided, the approximate Jacobian is used.

r: If the approximate root does not lie within a radius r of the initial guess x0, an exception is raised. i.e. if || x_k - x0 || > r

A warning is given if the slope of the function is close to zero, because convergence is very slow or not possible in this case, for example if we start our initial guess at a stationary point.
____________________

functions.py:

Contains a class Polynomial which is a callable polynomial object. Can be called as p=Polynomial([1, 2, 3]).
The function ApproximateJacobian returns an approximation of the Jacobian Df(x)as a numpy matrix.

____________________

testNewton.py:

Contains several tests for newton.py

testLinear: To test linear polynomial
testQuadraticbad: Test quadratic with both roots equal, so slope is 0 at the root.
testQuadraticgood: Test quadratic with distinct roots.
test_np: Testing with a numpy array instead of matrix.
test_2d: Testing a 2D function.
test_analytic_2d: Testing a 2d function with analytic Jacobian given.
test_analytic_1d: Testing a 1d function with analytic Jacobian given.
test_radius: Testing that the code throws an exception if the approximated root is not within a radius r of the initial guess.
test_radius1: Testing that the code throws an exception if the approximated root is not within a radius r of the initial guess.
test_analytic_numeric_1d: Testing that the analytic Jacobian is more accurate than the approximate one.
test_nonlinear: Testing a 2d nonlinear function.
test_3d: Testing a 3d function.
testNewtonStep: Testing that a single step of Newton performs correctly.
testNewtonStep1: Testing that a single step of Newton performs correctly.
testConverge: Test that the code throws an exception if it doesn't converge in the given number of iterations.

________________

testFunctions:

testApproxJacobian1: Testing that the approximate jacobian is correct for a 1d function.
testApproxJacobian2: Testing that the approximate jacobian is correct for a 2d function.
testApproxJacobian3: Testing that the approximate jacobian is correct for a 3d function.
testPolynomial: Test the Polynomial class.
testAnalytic: Test that analytic Jacobian gives the correct result by comparing with the
approximate Jacobian.
