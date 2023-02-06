#    s (float): current price
#    k (float): strike price
#    r (float): risk free rate
#    vol (float): underlying volatility
#    t (float): remaining time to expiry

import numpy as np
import scipy.stats as spstats
from scipy.stats import norm


def call_price(
    s: float, k: float, r: float, vol: float, t: float
) -> float:
    d1 = (np.log(s / k) + (r + (vol ** 2) / 2) * t) / (vol * (np.sqrt(t)))
    d2 = d1 - (vol * (np.sqrt(t)))
    return s * norm.cdf(d1) - k * np.exp(-r * t) * norm.cdf(d2)


def put_price(
    s: float, k: float, r: float, vol: float, t: float
) -> float:
    d1 = (np.log(s / k) + (r + (vol ** 2) / 2) * t) / (vol * (np.sqrt(t)))
    d2 = d1 - (vol * (np.sqrt(t)))
    return k * np.exp(-r * t) * norm.cdf(-d2) - s * norm.cdf(-d1)

def call_delta(
    s: float, k: float, r: float, vol: float, t: float
) -> float:
    d1 = (np.log(s / k) + (r + (vol ** 2) / 2) * t) / (vol * (np.sqrt(t)))
    return norm.cdf(d1)

def put_delta(
    s: float, k: float, r: float, vol: float, t: float
) -> float:
    d1 = (np.log(s / k) + (r + (vol ** 2) / 2) * t) / (vol * (np.sqrt(t)))
    return -norm.cdf(-d1)

def gamma(
    s: float, k: float, r: float, vol: float, t: float
) -> float:
    d1 = (np.log(s / k) + (r + (vol ** 2) / 2) * t) / (vol * (np.sqrt(t)))
    return norm.pdf(d1) / (s * vol * np.sqrt(t))

def vega(
    s: float, k: float, r: float, vol: float, t: float
) -> float:
    d1 = (np.log(s / k) + (r + (vol ** 2) / 2) * t) / (vol * (np.sqrt(t)))
    return s * norm.pdf(d1) * np.sqrt(t)

def call_theta(
    s: float, k: float, r: float, vol: float, t: float
) -> float:
    d1 = (np.log(s / k) + (r + (vol ** 2) / 2) * t) / (vol * (np.sqrt(t)))
    d2 = d1 - (vol * (np.sqrt(t)))
    return (-s * norm.pdf(d1) * vol) / (2 * np.sqrt(t)) - r * k * np.exp(-r * t) * norm.cdf(d2)

def put_theta(
    s: float, k: float, r: float, vol: float, t: float
) -> float:
    d1 = (np.log(s / k) + (r + (vol ** 2) / 2) * t) / (vol * (np.sqrt(t)))
    d2 = d1 - (vol * (np.sqrt(t)))
    return (-s * norm.pdf(d1) * vol) / (2 * np.sqrt(t)) + r * k * np.exp(-r * t) * norm.cdf(-d2)

def call_rho(
    s: float, k: float, r: float, vol: float, t: float
) -> float:
    d1 = (np.log(s / k) + (r + (vol ** 2) / 2) * t) / (vol * (np.sqrt(t)))
    d2 = d1 - (vol * (np.sqrt(t)))
    return k * t * np.exp(-r * t) * norm.cdf(d2)

def put_rho(
    s: float, k: float, r: float, vol: float, t: float
) -> float:
    d1 = (np.log(s / k) + (r + (vol ** 2) / 2) * t) / (vol * (np.sqrt(t)))
    d2 = d1 - (vol * (np.sqrt(t)))
    return -k * t * np.exp(-r * t) * norm.cdf(-d2)