
# newton - Newton-Raphson solver
#


import numpy as N
import functions as F

class Newton(object):
    def __init__(self, f, tol=1.e-6, maxiter=20, dx=1.e-6, Df=None, r=None):
        """Return a new object to find roots of f(x) = 0 using Newton's method.
        tol:     tolerance for iteration (iterate until |f(x)| < tol)
        maxiter: maximum number of iterations to perform
        dx:      step size for computing approximate Jacobian"""
        self._f = f
        self._tol = tol
        self._maxiter = maxiter
        self._dx = dx
        if r!=None:
            self._r=r
        else:
            self._r=10
        if Df!=None:
            self._Df=Df
        else: self._Df=0
    def solve(self, x0):
        """Return a root of f(x) = 0, using Newton's method, starting from
        initial guess x0"""
        x = x0
        for i in xrange(self._maxiter):
            fx = self._f(x)
            if N.linalg.norm(fx) < self._tol:
                return x
            x = self.step(x, fx)
            if N.linalg.norm(x-x0) > self._r: #pick a suitable  r, default = 10
                print 'norm(x0-x) is %f' % N.linalg.norm(x-x0)
                raise Exception (" norm(xo-x) is not within the threshold")
        if N.linalg.norm(fx) > self._tol:
            raise Exception( "Did not converge in chosen number of iterations. Change initial guess or increase number of iterations or tolerance")
        return x

    def step(self, x, fx=None):
        """Take a single step of a Newton method, starting from x
        If the argument fx is provided, assumes fx = f(x)"""
        
        if fx is None:
            fx = self._f(x)
        if self._Df==0:
            Df_x = F.ApproximateJacobian(self._f, x, self._dx)
        #    if foo:
        #        print 'Using numerical Jacobian'
        #        foo=0
        else: 
            #Df_x=self._Df(x)
            Df_x=F.AnalyticJacobian(self._Df, x)
         #   if foo:
        #print 'Using analytic Jacobian'
         #       foo=0
        if N.linalg.norm(Df_x)<1.e-6:
            print "Warning: Slope is almost zero"
        h = N.linalg.solve(N.matrix(Df_x), N.matrix(fx))
        return x - h
