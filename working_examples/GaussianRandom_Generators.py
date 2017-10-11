'''
Created on Jan 11, 2016

@author: vital
'''

import  numpy               as np
import  matplotlib.pyplot   as plt
from    scipy.integrate     import simps

def GaussianCurve(x, A, mu, sigma):
    
    y = A * np.exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma))
    
    return y


x = np.linspace(-10, 10, 300)

y1 = GaussianCurve(x, 5, 0, 1)
y2 = GaussianCurve(x, 10, 0, 1)

plt.plot(x, y1)
plt.plot(x, y2)


A1 = simps(y1, x)
A2 = simps(y2, x)

print 'A1', A1
print 'A2', A2
print 'Difference', A1*2 - A2

plt.show()


# mu, sigma = 0, 1 # mean and standard deviation
# s = np.random.normal(mu, sigma, 10000)
# s2 = np.random.randn(10000)
# 
# 
# s = np.random.normal(10, 2, 10000)
# s2 = 2*np.random.randn(10000) + 10
# 
# count, bins, ignored = plt.hist(s, 30)
# count, bins, ignored = plt.hist(s2, 30, color='red')
# 
# 
# 
# # plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2) ), linewidth=2, color='r')
# plt.show()