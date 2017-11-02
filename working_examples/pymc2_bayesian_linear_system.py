import pymc as pm
import numpy as np
from scipy.optimize import lsq_linear, nnls

def model(a_matrix, b_vector): 
      
    x_coeffs = [pm.Uniform('x_coeffs_%i' % i, 0.0, 5.00) for i in range(len(x_true))]
       
    @pm.deterministic(plot=False)
    def linear_solver(x_coeffs=x_coeffs, a_matrix=a_matrix):
        solver_solution = a_matrix.dot(x_coeffs)
        return solver_solution
      
    @pm.stochastic(observed=True)
    def likelihood(value=b_vector, fit_results = linear_solver, sigma=err):
        chiSq = np.sum(np.square(fit_results - value) / np.square(sigma))
        return - chiSq / 2
      
    return locals()

x_true  = np.array([3,2,1])
a       = np.array([[10, -5, 0], [15, 20, 0],
                    [12, 5, -20]])
b       = a.dot(x_true)
err     = 0.05

#Least squares
non_negative_least_squares = nnls(a, b)[0]

#Bayesian
MDL1 = pm.MCMC(model(a, b))
MDL1.sample(10000, 5000, 1)

#Invert matrix
x_inversion = np.dot(np.linalg.inv(np.dot(a.T, a)), np.dot(a.T, b))

#Compare
for i in range(len(b)):
    bayesian_fit_coeff = MDL1.stats()['x_coeffs_'+str(i)]['mean']
    print  x_true[i], bayesian_fit_coeff, non_negative_least_squares[i], x_inversion[i]


