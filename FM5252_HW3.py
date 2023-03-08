# S = underlying price
# K = strike price
# r = risk-free rate
# T = tenor
# sigma = initial volatility, or ideal volatility
# OP = option price
# call,put = option_type

import numpy as np
from math import pi
import scipy.stats as spstat

def black_scholes(S, K, T, r, sigma, option_type):
    d1 = (np.log(S/K) + (r + 0.5*sigma**2)*T) / (sigma*np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)
    if option_type == 'call':
        return S*spstat.norm.cdf(d1) - K*np.exp(-r*T)*spstat.norm.cdf(d2)
    elif option_type == 'put':
        return K*np.exp(-r*T)*spstat.norm.cdf(-d2) - S*spstat.norm.cdf(-d1)

def bisection_implied_volatility(S, K, T, r, OP, option_type):
    tol=1e-6
    max_iter=100
    vol_min = 0.01
    vol_max = 50.0

    for i in range(max_iter):
        vol_mid = (vol_min + vol_max) / 2

        price_mid = black_scholes(S, K, T, r, vol_mid, option_type)
        price_low = black_scholes(S, K, T, r, vol_min, option_type)

        if abs(price_mid - OP) < tol:
            return vol_mid

        if (price_low - OP) * (price_mid - OP) > 0:
            vol_min = vol_mid
            
        else:
            vol_max = vol_mid

    return vol_mid

def newton_implied_volatility(S, K, T, r, OP, option_type):
    tol=1e-6
    max_iter=100
    x_0 = .5
    sigma = x_0

    for i in range(max_iter):
        price = black_scholes(S, K, T, r, sigma, option_type)
        d1 = (np.log(S/K) + (r + 0.5*sigma**2)*T) / (sigma*np.sqrt(T))
        vega = S*(np.exp(-.5*d1**2)/np.sqrt(2*pi))*np.sqrt(T)

        diff = OP - price

        if abs(diff) < tol:
            return sigma
        else:
            sigma = sigma + diff/vega
            return sigma

    return sigma

def SVIsmile(k,parameter):
    a = parameter[0]
    b = parameter[1]
    rho = parameter[2]
    m = parameter[3]
    sigma = parameter[4]
    return a + b*(rho*(k - m) + np.sqrt((k - m)**2 + sigma**2))

parameter = [100,100,0.05,0.5,0.05]
k = np.linspace(-3,3,5)

import matplotlib.pyplot as plt

plt.plot(k, SVIsmile(k,parameter), label = parameter, marker = '.')
plt.xlabel('Moneyness')
plt.ylabel('Variance')
plt.show()