'''
Created on Oct 26, 2015

@author: vital
'''

    import numpy as np
    import pymc as pm
    import pylab as plt
    import seaborn as sb
    
    #True m, n parameters
    n_true = 10
    m_slope = 0.5
    
    #Generate points x and y coordinates
    x = np.linspace(0, 10, 100)
    y = m_slope * x + n_true
    
    #Add noise to data
    x += 0.5 * (np.random.rand(len(x)) * 2.0 - 1.0)
    y += 1.0 * (np.random.rand(len(x)) * 3 - 1.0)
    
    #Priors
    m   = pm.Normal('m', 0.0, tau=1e-5, value=0)
    n   = pm.Normal('n', 0.0, tau=1e-5, value=0)
    tau = pm.Gamma('tau', 0.01, 0.01, value=0.01)
    
    #Theoretical values
    y_model = m*x + n
    
    #Defining likelihood
    likelihood = pm.Normal('likelihood', y_model, tau=tau, value=y, observed=True) 
    
    #Launch MCMC
    mcmc = pm.MCMC(locals())
    mcmc.sample(20000, burn=10000)
    
    #Display results
    print m.summary()
    print n.summary()
    
    #Display results
    m_Inference =  mcmc.trace('m').stats()['mean']
    n_Inference =  mcmc.trace('n').stats()['mean'] 
    
    mError_Inference = mcmc.trace('m').stats()['standard deviation']
    nError_Inference = mcmc.trace('n').stats()['standard deviation'] 
    
    sb.color_palette('colorblind')
    plt.ylabel('y')
    plt.xlabel('x')
    plt.title('Inference regression: ' + r'$m_{inf}=$' + str(round(m_Inference,3)) + r'$\pm$'+str(round(mError_Inference, 3)) + ' ' + r'$n_{inf}=$' + str(round(n_Inference,3)) + r'$\pm$'+str(round(nError_Inference,3)))
    plt.scatter(x, y)
    plt.plot(x, m_Inference * x + n_Inference)
    plt.show()