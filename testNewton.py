#!/usr/bin/env python

import newton
import unittest
import numpy as N
import math
import functions as F

class TestNewton(unittest.TestCase):
    def testLinear(self):
        #f = lambda x : 3.0 * x + 6.0
        f=F.Polynomial([0,3,6])

        solver = newton.Newton(f, tol=1.e-15, maxiter=10)
        x = solver.solve(2.0)
        self.assertEqual(x, -2.0)
#        print x

    def testQuadraticbad(self):
#        f = lambda x : x*x + 6*x + 9
        f= F.Polynomial([1,6,9])
        _Df=F.Polynomial([0,2,6])
        solver = newton.Newton(f, tol=1.e-15, maxiter=200,Df=_Df)
        x = solver.solve(-1.0)
        self.assertAlmostEqual(x, -3.0)
#        print x

    def testQuadraticgood(self):
#        f = lambda x : x*x - 7*x + 12
        f=F.Polynomial([1,-7,12])
        solver = newton.Newton(f, tol=1.e-14, maxiter=1000)
        x = solver.solve(3.6)
        self.assertAlmostEqual(x, 4.0)
#        print x
              
    def test_np(self):
#        f = N.array([F.Polynomial([1.,0.,-4.]), F.Polynomial([0.,1.,-4.])])
        f= lambda x: N.array([[x*x-7*x+12],])
        solver = newton.Newton(f, tol=1.e-6, maxiter=30)
        x0=N.array([ [3.6],])
        x = solver.solve(x0)
        N.testing.assert_array_almost_equal(x, N.array([[4.],]))
#   
    def test_2d(self):
#        f = N.array([F.Polynomial([1.,0.,-4.]), F.Polynomial([0.,1.,-4.])])
        f= lambda x: N.matrix([[x[0,0]-x[1,0]-3],[ x[0,0]+x[1,0]-5]])

        solver = newton.Newton(f, tol=1.e-6, maxiter=30)
        x0=N.matrix([[3],[2]])
        x = solver.solve(x0)
#        print x
        N.testing.assert_array_almost_equal(x, N.matrix([[4.], [1.]]))

    def test_analytic_2d(self):
#        f = N.array([F.Polynomial([1.,0.,-4.]), F.Polynomial([0.,1.,-4.])])
        f= lambda x: N.matrix([[x[0,0]-2*x[1,0]-3],[ x[0,0]+x[1,0]-6]])
        _Df= lambda x: N.matrix([[1, -2],[1, 1]])
        solver = newton.Newton(f, tol=1.e-6, maxiter=30, Df=_Df)
        x0=N.matrix([[3],[2]])
        x = solver.solve(x0)
#        print x
        N.testing.assert_array_almost_equal(x, N.matrix([[5.], [1.]]))
   
    def test_analytic_1d(self):
        f= F.Polynomial([1,0,-4])
#        print f
        _Df= lambda x: 2*x
        solver = newton.Newton(f, tol=1.e-8, maxiter=30,Df=_Df)
        x = solver.solve(-3.0)
        self.assertAlmostEqual(x, -2.0)
        print x

    def test_radius(self):
        f= F.Polynomial([1,0,-4])
#        print f
        _Df= lambda x: 2*x
        solver = newton.Newton(f, tol=1.e-8, maxiter=30,Df=_Df,r=50)
        x = solver.solve(-50.0)
        self.assertAlmostEqual(x, -2.0)
#        print x

    def test_radius1(self):
        f= lambda x: x/(x*x + 1)
        solver = newton.Newton(f, tol=1.e-8, maxiter=30,r=5)
        x = solver.solve(0.5)
        self.assertAlmostEqual(x, 0.)
#        print x

   
    def test_analytic_numeric_1d(self):
        f= F.Polynomial([1,0,-4])
#        print f
        _Df= lambda x: 2*x
        solver_an = newton.Newton(f, tol=1.e-8, maxiter=1,Df=_Df)
        x_an = solver_an.step(15.0)
        solver_num=  newton.Newton(f, tol=1.e-8, maxiter=1)
        x_num = solver_num.step(15.)
#        self.assertAlmostEqual(x_an, -2.0)
#        self.assertAlmostEqual(x_num, -2.)
        self.assertTrue(abs(x_an-2)< abs(x_num-2))
        print x_an, x_num
   



    def test_nonlinear(self):
#        f = N.array([F.Polynomial([1.,0.,-4.]), F.Polynomial([0.,1.,-4.])])
        f= lambda x: N.matrix([[math.tan(x[1,0]) ],[x[0,0]-4 ]])

        solver = newton.Newton(f, tol=1.e-6, maxiter=100)
        x0=N.matrix([[3.],[3.]])
        x = solver.solve(x0)
#        print x
        N.testing.assert_array_almost_equal(x, N.matrix([[4.], [math.pi]]))
   


    def test_3d(self):
#        f = N.array([F.Polynomial([1.,0.,-4.]), F.Polynomial([0.,1.,-4.])])
        f= lambda x: N.matrix([[x[0,0]+x[1,0]+x[2,0]-5],[ x[0,0]-x[1,0]-2*x[2,0]-4],[3*x[1,0]+x[2,0]-3]])

        solver = newton.Newton(f, tol=1.e-6, maxiter=30)
        x0=N.matrix([[3.],[2.],[0.]])
        x = solver.solve(x0)
#        print x
        N.testing.assert_array_almost_equal(x, N.matrix([[30./7.], [8./7.],[-3./7.]]))



    def testNewtonStep(self):
        f = lambda x : math.sin(x)
        x0=6
        solver = newton.Newton(f,tol=1.e-15, maxiter=1)
        x=solver.step(x0)
#        print x 
        self.assertTrue(abs(x-2*math.pi)<abs(2*math.pi-x0))

    def testNewtonStep1(self):
        f = lambda x : x*x-1
        x0=5
        solver = newton.Newton(f,tol=1.e-15, maxiter=1)
        x=solver.step(x0)
#        print 'In NewtonSTep1 %f' % x
        self.assertTrue(abs(-1-x)<abs(-1-x0))

    def testConverge(self):
        f = lambda x : x*x-1
        solver = newton.Newton(f, tol=1.e-15, maxiter=30)
        x = solver.solve(0.5)    
        self.assertEqual(x, 1.0)
#        print x

if __name__ == "__main__":
    unittest.main()
