import numpy as np
import lmfit
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

def gaussian_curve(variables_dict, x, zerolev):
    
    y_model = variables_dict['A'] * np.exp(-(x-variables_dict['mu'])*(x-variables_dict['mu'])/(2*variables_dict['sigma']*variables_dict['sigma'])) + zerolev
    
    return y_model

def lmfit_gaussian_Residual(params, x, y, zerolev, err):
 
    return (gaussian_curve(params.valuesdict(), x, zerolev) - y) / err

m, n = 0.2, 0.1

x_real = np.linspace(-0.5, 0.5, 100)
zerolevel = m * x_real + n

real_values = {}
real_values['A']        = 1
real_values['mu']       = 0
real_values['sigma']    = 0.1

y_real = gaussian_curve(real_values, x_real, zerolev=zerolevel) 

x_obs = x_real + np.random.normal(0.0, 0.005, len(x_real))
y_obs = y_real + np.random.normal(0.0, 0.006, len(y_real))


index_sample = np.searchsorted(x_obs, [-0.15, -0.08, -0.05, 0.0, 0.03, 0.11, 0.22]) 

params = lmfit.Parameters()
params.add('A', value = 1)
params.add('mu', value = 0)
params.add('sigma', value = 1, min=0.0)

left_continuum_indx  = (x_obs >= -0.5) & (x_obs <= -0.3)
right_continuum_indx = (x_obs >= 0.3) & (x_obs <= 0.5)

x_cont   = np.concatenate((x_obs[left_continuum_indx], x_obs[right_continuum_indx]))
y_cont   = np.concatenate((y_obs[left_continuum_indx], y_obs[right_continuum_indx]))

# left_index = np.searchsorted(x_obs, x_obs[])
linear_continuum    = m * x_cont + n
err_continuum       = np.std(linear_continuum - y_cont)

x_sample        = x_obs[index_sample]
y_sample        = y_obs[index_sample]
zerolevel       = zerolevel[index_sample]
sample_length   = len(x_sample) 

Single_fit_Output = lmfit.minimize(lmfit_gaussian_Residual, params, args=(x_sample, y_sample, zerolevel, err_continuum))
print 'single', Single_fit_Output.params.valuesdict()


iterations                  = 10
output_vector_dict          = {}
output_vector_dict['A']     = np.zeros(iterations)
output_vector_dict['mu']    = np.zeros(iterations)
output_vector_dict['sigma'] = np.zeros(iterations)

for i in range(iterations):
    if i == 0:
        y_new = y_sample
    else:
        y_new = y_sample + np.random.normal(0.0, err_continuum, sample_length)
        
    fit_Output      = lmfit.minimize(lmfit_gaussian_Residual, params, args=(x_sample, y_new, zerolevel, err_continuum))
    iteration_ouput = fit_Output.params.valuesdict()
    print iteration_ouput
    
    output_vector_dict['A'][i]     = iteration_ouput['A']
    output_vector_dict['mu'][i]    = iteration_ouput['mu']
    output_vector_dict['sigma'][i] = iteration_ouput['sigma']


# print output_vector_dict['A']
# print output_vector_dict['mu']
# print output_vector_dict['sigma']
mean_values = {}
mean_values['A']        = np.mean(output_vector_dict['A'])
mean_values['mu']       = np.mean(output_vector_dict['mu'])
mean_values['sigma']    = np.mean(output_vector_dict['sigma'])

print 'Mean values'
print mean_values


fit_gaussian = gaussian_curve(fit_Output.params.valuesdict(), x_obs, zerolev=m * x_obs + n) 
fit_gaussian2 = gaussian_curve(mean_values, x_obs, zerolev=m * x_obs + n) 


plt.plot(x_real, y_real)
plt.plot(x_obs, y_obs)
plt.plot(x_obs, fit_gaussian)
plt.plot(x_obs, fit_gaussian2)

plt.plot(x_sample, y_sample,'ro')
# plt.plot(x_cont, y_cont,'bo')

plt.show()