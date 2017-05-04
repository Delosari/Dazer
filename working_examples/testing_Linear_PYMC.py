import pymc     as pm
import pymc3    as pm3
import numpy    as np
import matplotlib.pyplot as plt

def model_gauss(x, y, y_error): 
    amp     = pm.Uniform('amp', 0.05, 0.4, value= 0.15)
    size    = pm.Uniform('size', 0.5, 2.5, value= 1.0)
    ps      = pm.Normal('ps', 0.13, 40, value=0.15)

    @pm.deterministic(plot=False)
    def gauss_likelihood(x=x, amp=amp, size=size, ps=ps):
        e = -1*(np.pi**2*size*x/(3600.*180.))**2/(4.*np.log(2.))
        return amp*np.exp(e)+ps
    
    y = pm.Normal('y', mu=gauss_likelihood, tau=1.0/y_error**2, value=y, observed=True)
    
    return locals()

def model_gauss_err(x, y, x_errObs, y_errObs): 
    amp     = pm.Uniform('amp', 0.05, 0.4, value= 0.15)
    size    = pm.Uniform('size', 0.5, 2.5, value= 1.0)
    ps      = pm.Normal('ps', 0.13, 40, value=0.15)
    
    @pm.deterministic(plot=False)
    def yerr_from_xerr(x_errObs=x_errObs, amp=amp, size=size, ps=ps):
        return amp * np.exp(-1*(np.pi**2/(3600.*180.)*size*x_errObs)**2/(4.*np.log(2.))) + ps
    
    @pm.deterministic(plot=False)
    def gauss_likelihood(x=x, amp=amp, size=size, ps=ps):
        e = -1*(np.pi**2*size*x/(3600.*180.))**2/(4.*np.log(2.))
        return amp*np.exp(e)+ps
        
    y = pm.Normal('y', mu=gauss_likelihood, tau=1.0/(y_errObs+yerr_from_xerr)**2, value=y, observed=True)
    
    return locals()

def model_lineal_jake(x, y):  
    
    # Define the variables needed for the routine, with their prior distributions
    alpha = pm.Uniform('alpha', -100, 100)
    
    @pm.stochastic(observed=False)
    def beta(value=0):
        return -1.5 * np.log(1 + value ** 2)

    @pm.stochastic(observed=False)
    def sigma(value=1):
        return -np.log(abs(value))
    
    # Define the form of the model and likelihood
    @pm.deterministic
    def y_likelihood(x=x, alpha=alpha, beta=beta):
        return alpha + beta * x

    y = pm.Normal('y', mu=y_likelihood, tau=1. / sigma ** 2, observed=True, value=y)

    return locals()

np.random.seed(12345)

#-------------Linear fit Jake----------------

alpha_true, beta_true = 25, 0.5
x_true      = 100 * np.random.random(20)
y_true      = alpha_true + beta_true * x_true

x           = np.random.normal(x_true, 10)
y           = np.random.normal(y_true, 10)

MDL = pm.MCMC(model_lineal_jake(x, y))
MDL.sample(20000, 10000, 1)

y_min       = MDL.stats()['y_likelihood']['quantiles'][2.5]
y_max       = MDL.stats()['y_likelihood']['quantiles'][97.5]
y_fit       = MDL.stats()['y_likelihood']['mean']
alpha_min   = MDL.stats()['alpha']['quantiles'][2.5]
alpha_max   = MDL.stats()['alpha']['quantiles'][97.5]
alpha_mean  = MDL.stats()['alpha']['mean']
beta_min    = MDL.stats()['beta']['quantiles'][2.5]
beta_max    = MDL.stats()['beta']['quantiles'][97.5]
beta_mean   = MDL.stats()['beta']['mean']

# print MDL.trace('y_likelihood')[:]
# print type( MDL.trace('y_likelihood')[:]),  MDL.trace('y_likelihood')[:].shape
# xfit = np.linspace(-20, 120, 10)
# y_otrofit = MDL.stats()['alpha'] + MDL.stats()['beta'] * xfit



plt.plot(x, y, 'ok', label='Obs data')
# plt.errorbar(x,y,xerr=x_error,yerr=y_error, color='r', marker='.', ls='None',label='Observed')
plt.plot(x, y_fit,'-', ms=5, mew=1, label='Fit')
plt.plot(x, alpha_mean + x * beta_mean, '-', label='Fit mean params')

plt.plot(x, alpha_max + x * beta_max, '-', label='Fit max')
plt.plot(x, alpha_min + x * beta_min, '-', label='Fit min')

x_regression = np.linspace(x.min(), x.max(), 100)
y_regression = MDL.stats()['alpha']['mean'] + MDL.stats()['beta']['mean'] * x_regression
# plt.fill_between(x, y_min, y_max, color='0.5', alpha=0.5)
plt.legend()
plt.show()

# #-------------Gaussian fit----------------
# 
# x = np.arange(5,400,10) * 1e3
# 
# amp_true    = 0.2
# size_true   = 1.8
# ps_true     = 0.1
# 
# gauss       = lambda x,amp,size,ps: amp * np.exp(-1*(np.pi**2/(3600.*180.)*size*x)**2/(4.*np.log(2.))) + ps
# y_true      = gauss(x=x,amp=amp_true, size=size_true, ps=ps_true )
# 
# noise       = np.random.normal(size=len(x)) * .02 
# y           = y_true + noise 
# y_error     = np.ones_like(y_true)*0.05*y.max()
# x_error     = np.ones_like(x)*0.02*x.max()
# 
# #Define the model_gauss/function to be fitted.
# MDL = pm.MCMC(model_gauss(x, y, y_error))
# # MDL = pm.MCMC(model_gauss_err(x, y, x_error, y_error))
# MDL.sample(20000, 10000, 1)
# 
# #Plot the function
# y_min = MDL.stats()['gauss_likelihood']['quantiles'][2.5]
# y_max = MDL.stats()['gauss_likelihood']['quantiles'][97.5]
# y_fit = MDL.stats()['gauss_likelihood']['mean']
# plt.plot(x,y_true,'b', marker='None', ls='-', lw=1, label='True')
# plt.errorbar(x,y,xerr=x_error,yerr=y_error, color='r', marker='.', ls='None',label='Observed')
# plt.plot(x,y_fit,'k', marker='+', ls='None', ms=5, mew=1, label='Fit')
# plt.fill_between(x, y_min, y_max, color='0.5', alpha=0.5)
# plt.legend()
# plt.show()
# 
# # #Plot the data
# # pm.Matplot.plot(MDL)
# # plt.show()


