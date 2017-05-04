
import itertools

import matplotlib.pyplot as plt, seaborn as sb
import numpy as np
import pymc as pm


#set random seed for reproducibility
np.random.seed(12345)

#x axis range
x_true      = np.arange(5,400,10)*1e3

# Parameters for gaussian
Amplitude   = 0.2
size_true   = 1.8
zero_level  = 0.1

# Gaussian function
gauss       = lambda x_true,amp,size,ps: amp*np.exp(-1*(np.pi**2/(3600.*180.)*size*x_true)**2/(4.*np.log(2.)))+ps
f_true      = gauss(x_true=x_true,amp=Amplitude, size=size_true, ps=zero_level )

# add noise to the data points
noise       = np.random.normal(size=len(x_true)) * .02 
f           = f_true + noise 
f_error     = np.ones_like(f_true)*0.05*f.max()

# add noise to observed x_true values
x_obs       = pm.rnormal(mu=x_true, tau=(1e4)**-2)

# define the model/function to be fitted.
def model(x_obs, f): 
    amp     = pm.Uniform('amp',  0.05,   0.4,    value = 0.15)
    size    = pm.Uniform('size', 0.5,    2.5,    value = 1.0)
    ps      = pm.Normal('ps',   0.13,   40,     value = 0.15)
    x_pred  = pm.Normal('x_true', mu=x_obs, tau=(1e4)**-2) # this allows error in x_obs
    

    @pm.deterministic(plot=False)
    def gauss(x_true=x_pred, amp=amp, size=size, ps=ps):
        e = -1*(np.pi**2*size*x_true/(3600.*180.))**2/(4.*np.log(2.))
        return amp*np.exp(e)+ps
    
    y   = pm.Normal('y', mu=gauss, tau=1.0/f_error**2, value=f, observed=True)
    
    return locals()

palette = itertools.cycle(sb.color_palette('BuGn'))

print sb.color_palette('colorblind', 1)
print sb.color_palette('colorblind', 4)

#Running the MCMC hammer
MDL = pm.MCMC(model(x_obs, f))
MDL.use_step_method(pm.AdaptiveMetropolis, MDL.x_pred) # use AdaptiveMetropolis to "learn" how to step
MDL.sample(20000, 10000, 10)  # run chain longer since there are more dimensions

#Plotting graphs
y_min = MDL.stats()['gauss']['quantiles'][2.5]
y_max = MDL.stats()['gauss']['quantiles'][97.5]
y_fit = MDL.stats()['gauss']['mean']
x_fit = MDL.x_pred.trace().mean(0)
plt.plot(x_true,f_true,'b', marker='.',color=next(palette), ls='-', lw=1, label='True')
plt.errorbar(x_obs, f, yerr=f_error, marker='.',color=next(palette), ls='None', label='Observed')
plt.plot(x_fit,y_fit,'k', marker='+', ls='None', color=next(palette), ms=5, mew=1, label='Fit')
plt.plot(x_fit*1.1,y_fit*1.1,'k', marker='+', ls='None', color=next(palette), ms=5, mew=1, label='Fit2')

plt.fill_between(x_fit, y_min, y_max, alpha=0.5, color=next(palette))
plt.legend()

plt.show()


# # define the model/function to be fitted.
# def model(x_true, f): 
#     amp = pm.Uniform('amp', 0.05, 0.4, value= 0.15)
#     size = pm.Uniform('size', 0.5, 2.5, value= 1.0)
#     ps = pm.Normal('ps', 0.13, 40, value=0.15)
# 
#     @pm.deterministic(plot=False)
#     def gauss(x_true=x_true, amp=amp, size=size, ps=ps):
#         e = -1*(np.pi**2*size*x_true/(3600.*180.))**2/(4.*np.log(2.))
#         return amp*np.exp(e)+ps
#     y = pm.Normal('y', mu=gauss, tau=1.0/f_error**2, value=f, observed=True)
#     return locals()
# 
# MDL = pm.MCMC(model(x_true,f))
# MDL.sample(20000, 10000, 1)
# 
# # extract and plot results
# y_min = MDL.stats()['gauss']['quantiles'][2.5]
# y_max = MDL.stats()['gauss']['quantiles'][97.5]
# y_fit = MDL.stats()['gauss']['mean']
# plt.plot(x_true,f_true,'b', marker='None', ls='-', lw=1, label='True')
# plt.errorbar(x_true,f,yerr=f_error, color='r', marker='.', ls='None', label='Observed')
# plt.plot(x_true,y_fit,'k', marker='+', ls='None', ms=5, mew=1, label='Fit')
# plt.fill_between(x_true, y_min, y_max, color='0.5', alpha=0.5)
# plt.legend()
# 
# pm.Matplot.plot(MDL)
# 
# plt.show()