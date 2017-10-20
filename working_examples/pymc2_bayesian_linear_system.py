import pymc as pm
import numpy as np

a = np.array([[1,2,3,4],
             [5,6,7,8],
             [9,10,11,12],
             [1,2,3,4],
             [0.1,0.2,0.3,0.4]])

print a.shape

print np.delete(a, [0,3], axis=0)
print np.delete(a, [0,3], axis=1)


# def generate_pymc_distr_array(variable_name, size):
#     
#     x_coeffs = np.empty(size)
#     for i in range(size):
#         x_coeffs[i] = pm.Uniform(variable_name % i,  0.0, 5.00)
#         
#     return x_coeffs
# 
# def model(a_matrix, b_vector): 
#      
#     x_coeffs = [pm.Uniform('x_coeffs_%i' % i, 0.0, 5.00) for i in range(len(x_true))]
#       
#     @pm.deterministic(plot=False)
#     def linear_solver(x_coeffs=x_coeffs, a_matrix=a_matrix):
#         solver_solution = a_matrix.dot(x_coeffs)
#         return solver_solution
#      
#     @pm.stochastic(observed=True)
#     def likelihood(value=b_vector, fit_results = linear_solver, sigma=err):
#         chiSq = np.sum(np.square(fit_results - value) / np.square(sigma))
#         return - chiSq / 2
#      
#     return locals()
# 
# def model2(a_matrix, b_vector): 
#     
#     x_coeffs = generate_pymc_distr_array('x_coeffs_%i', 3)
#     
#     @pm.deterministic(plot=False)
#     def linear_solver(x_coeffs=x_coeffs, a_matrix=a_matrix):
#         solver_solution = a_matrix.dot(np.array(x_coeffs))
#         return solver_solution
#     
#     @pm.stochastic(observed=True)
#     def likelihood(value=b_vector, fit_results = linear_solver, sigma=err):
#         chiSq = np.sum(np.square(fit_results - value) / np.square(sigma))
#         return - chiSq / 2
#     
#     return locals()
# 
# def model3(a_matrix, b_vector): 
#            
#     @pm.deterministic(plot=False)
#     def linear_solver(x_coeffs=x_coeffs, a_matrix=a_matrix):
#         solver_solution = a_matrix.dot(x_coeffs)
#         return solver_solution
#      
#     @pm.stochastic(observed=True)
#     def likelihood(value=b_vector, fit_results = linear_solver, sigma=err):
#         chiSq = np.sum(np.square(fit_results - value) / np.square(sigma))
#         return - chiSq / 2
#      
#     return locals()
# 
# def lo_sacamos_de_aqui(x_array):
#     print 'hola'
#     x_coeffs = [pm.Uniform('x_coeffs_%i' % i, 0.0, 5.00) for i in range(len(x_array))]
#     
#     return
#     
# def model4(a_matrix, b_vector): 
#     
#     @pm.deterministic
#     def x_coeffs(b_vector=b_vector):
#         x_coeffs = lo_sacamos_de_aqui(b_vector)
#         return x_coeffs
#         
#     @pm.deterministic
#     def linear_solver(x_coeffs=x_coeffs, a_matrix=a_matrix):
#         solver_solution = a_matrix.dot(x_coeffs)
#         return solver_solution
#      
#     @pm.stochastic(observed=True)
#     def likelihood(value=b_vector, fit_results = linear_solver, sigma=err):
#         chiSq = np.sum(np.square(fit_results - value) / np.square(sigma))
#         return - chiSq / 2
#      
#     return locals()
# 
# 
# x_true  = np.array([3,2,1])
# a       = np.array([[10, -5, 0], [15, 20, 0],
#                     [12, 5, -20]])
# b       = a.dot(x_true)
# err     = 0.05
# 
# MDL1 = pm.MCMC(model(a, b))
# MDL1.sample(10000, 5000, 1)
# 
# # MDL2 = pm.MCMC(model2(a, b))
# # MDL2.sample(10000, 5000, 1)
# # 
# # x_coeffs = [pm.Uniform('x_coeffs_%i' % i, 0.0, 5.00) for i in range(len(x_true))]
# # MDL3 = pm.MCMC([model3(a, b)] + x_coeffs)
# # MDL3.sample(10000, 5000, 1)
# 
# # MDL4 = pm.MCMC([model4(a, b)])
# # MDL4.sample(10000, 5000, 1)
# 
# 
# print 'True results', b
# print 'Fit results a', MDL1.stats()['linear_solver']['mean']
# # print 'Fit results b', MDL2.stats()['linear_solver']['mean']
# # print 'Fit results c', MDL3.stats()['linear_solver']['mean']
# # print 'Fit results d', MDL4.stats()['linear_solver']['mean']

