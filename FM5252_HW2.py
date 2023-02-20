#    s (float): current price
#    k (float): strike price
#    r (float): risk free rate
#    sigma (float): underlying volatility
#    t (float): remaining time to expiry
#    n (float): number of steps that should be generated

import numpy as np

# European call option 
def crr_eu_call(
    s: float, k: float, r: float, t: float, sigma: float, n: float
) -> float:

    dt = t/n
    u = np.exp(sigma*np.sqrt(dt))
    d = 1/u
    p = (np.exp(r*dt)-d)/(u-d)

    def call_value(t,j,s): 
        # 't' is current time step, 'j' is index of the node in 't'
        # both of them are default to zero
        if t == n:
            return max(s-k,0)
        else:
            up = s*u
            down = s*d

            up_option = call_value(t+1,j,up)
            down_option = call_value(t+1,j+1,down)

            return np.exp(-r*dt)*(p*up_option+(1-p)*down_option)

    call_price = call_value(0,0,s)

    return call_price

# European put option 
def crr_eu_put(
    s: float, k: float, r: float, t: float, sigma: float, n: float
) -> float:

    dt = t/n
    u = np.exp(sigma*np.sqrt(dt))
    d = 1/u
    p = (np.exp(r*dt)-d)/(u-d)

    def put_value(t,j,s): 
        if t == n:
            return max(k-s,0)
        else:
            up = s*u
            down = s*d

            up_option = put_value(t+1,j,up)
            down_option = put_value(t+1,j+1,down)

            return np.exp(-r*dt)*(p*up_option+(1-p)*down_option)

    put_price = put_value(0,0,s)

    return put_price

# American call option 
def crr_us_call(
    s: float, k: float, r: float, t: float, sigma: float, n: float
) -> float:

    dt = t/n
    u = np.exp(sigma*np.sqrt(dt))
    d = 1/u
    p = (np.exp(r*dt)-d)/(u-d)

    def call_value(t,j,s): 
        if t == n:
            return max(s-k,0)
        else:
            up = s*u
            down = s*d

            up_option = call_value(t+1,j,up)
            down_option = call_value(t+1,j+1,down)

            previous_exercise = max(s*(u**j)*(d**(n-j))-k,0)
            current_exercise = np.exp(-r*dt)*(p*up_option+(1-p)*down_option)

            return max(previous_exercise,current_exercise)

    call_price = call_value(0,0,s)

    return call_price

# American put option
def crr_us_put(
    s: float, k: float, r: float, t: float, sigma: float, n: float
) -> float:

    dt = t/n
    u = np.exp(sigma*np.sqrt(dt))
    d = 1/u
    p = (np.exp(r*dt)-d)/(u-d)

    def put_value(t,j,s): 
        if t == n:
            return max(k-s,0)
        else:
            up = s*u
            down = s*d

            up_option = put_value(t+1,j,up)
            down_option = put_value(t+1,j+1,down)

            previous_exercise = max(k-s*(u**j)*(d**(n-j)),0)
            current_exercise = np.exp(-r*dt)*(p*up_option+(1-p)*down_option)

            return max(previous_exercise,current_exercise)

    put_price = put_value(0,0,s)

    return put_price

# Greek values for Europran and American, call and put options

def eu_call_greeks(
        s: float, k: float, r: float, t: float, sigma: float, n: float
) -> float:

    dt = t/n
    u = np.exp(sigma*np.sqrt(dt))
    d = 1/u
    p = (np.exp(r*dt)-d)/(u-d)

    def call_value(t,j,s): 
        if t == n:
            return max(s-k,0)
        else:
            up = s*u
            down = s*d

            up_option = call_value(t+1,j,up)
            down_option = call_value(t+1,j+1,down)

            return np.exp(-r*dt)*(p*up_option+(1-p)*down_option)

    def delta(t,j,s):
        if t == n:
            return max(u*s*d**(n-1-j)-k,0) - max(d*s*u**(n-1-j)-k,0)
        else:
            return np.exp(-r*dt)*(p*delta(t+1,j+1,s)+(1-p)*delta(t+1,j,s))/(u-d)

    def gamma(t,j,s):
        if t == n:
            return 0
        else:
            return np.exp(-r*dt)*((p*gamma(t+1,j+1,s)+(1-p)*gamma(t+1,j,s))/(u-d)-delta(t+1,j+1,s)/(s*u**(n-t-1)*d**j-s*u**(n-t-1)*d**(j+1)))

    def theta(t,j,s):
        if t == n:
            return 0
        else:
            return (call_value(t+1,j,s)-call_value(t,j,s))/dt

    def vega(t,j,s):
        if t == n:
            return 0
        else:
            return np.exp(-r*dt)*((p*vega(t+1,j+1,s)+(1-p)*vega(t+1,j,s))/(u-d)-(call_value(t+1,j+1,s)-call_value(t+1,j-1,s))/(2*sigma*np.sqrt(dt)*s*u**(n-t-1)*d**(j-1)))

    def rho(t,j,s):
        if t == n:
            return 0
        else:
            return (call_value(t+1,j,s)-call_value(t,j,s))/dt

    delta_0 = delta(0,0,s)
    gamma_0 = delta(0,0,s)
    theta_0 = theta(0,0,s)
    vega_0 = vega(0,0,s)
    rho_0 = rho(0,0,s)
    return 'delta =',delta_0,'gamma =',gamma_0,'theta =',theta_0,'vega =',vega_0,'rho =',rho_0

def eu_put_greeks(
        s: float, k: float, r: float, t: float, sigma: float, n: float
) -> float:

    dt = t/n
    u = np.exp(sigma*np.sqrt(dt))
    d = 1/u
    p = (np.exp(r*dt)-d)/(u-d)

    def put_value(t,j,s): 
        if t == n:
            return max(k-s,0)
        else:
            up = s*u
            down = s*d

            up_option = put_value(t+1,j,up)
            down_option = put_value(t+1,j+1,down)

            return np.exp(-r*dt)*(p*up_option+(1-p)*down_option)

    def delta(t,j,s):
        if t == n:
            return max(k-u*s*d**(n-1-j),0) - max(k-d*s*u**(n-1-j),0)
        else:
            return np.exp(-r*dt)*(p*delta(t+1,j+1,s)+(1-p)*delta(t+1,j,s))/(u-d)

    def gamma(t,j,s):
        if t == n:
            return 0
        else:
            return np.exp(-r*dt)*((p*gamma(t+1,j+1,s)+(1-p)*gamma(t+1,j,s))/(u-d)-delta(t+1,j+1,s)/(s*u**(n-t-1)*d**j-s*u**(n-t-1)*d**(j+1)))

    def theta(t,j,s):
        if t == n:
            return 0
        else:
            return (put_value(t+1,j,s)-put_value(t,j,s))/dt

    def vega(t,j,s):
        if t == n:
            return 0
        else:
            return np.exp(-r*dt)*((p*vega(t+1,j+1,s)+(1-p)*vega(t+1,j,s))/(u-d)-(put_value(t+1,j+1,s)-put_value(t+1,j-1,s))/(2*sigma*np.sqrt(dt)*s*u**(n-t-1)*d**(j-1)))

    def rho(t,j,s):
        if t == n:
            return 0
        else:
            return (put_value(t+1,j,s)-put_value(t,j,s))/dt

    delta_0 = delta(0,0,s)
    gamma_0 = delta(0,0,s)
    theta_0 = theta(0,0,s)
    vega_0 = vega(0,0,s)
    rho_0 = rho(0,0,s)
    return 'delta =',delta_0,'gamma =',gamma_0,'theta =',theta_0,'vega =',vega_0,'rho =',rho_0

def us_call_greeks(
        s: float, k: float, r: float, t: float, sigma: float, n: float
) -> float:

    dt = t/n
    u = np.exp(sigma*np.sqrt(dt))
    d = 1/u
    p = (np.exp(r*dt)-d)/(u-d)

    def call_value(t,j,s): 
        if t == n:
            return max(s-k,0)
        else:
            previous_exercise = max(s*u**j*d**(n-j)-k,0)
            current_exercise = np.exp(-r*dt)*(p*call_value(t+1,j+1,s)+(1-p)*call_value(t+1,j,s))
            return max(previous_exercise,current_exercise)
    
    def delta(t,j,s):
        if t == n:
            return max(u*s*d**(n-1-j)-k,0) - max(d*s*u**(n-1-j)-k,0)
        else:
            return np.exp(-r*dt)*(p*delta(t+1,j+1,s)+(1-p)*delta(t+1,j,s))/(u-d)

    def gamma(t,j,s):
        if t == n:
            return 0
        else:
            return np.exp(-r*dt)*((p*gamma(t+1,j+1,s)+(1-p)*gamma(t+1,j,s))/(u-d)-delta(t+1,j+1,s)/(s*u**(n-t-1)*d**j-s*u**(n-t-1)*d**(j+1)))

    def theta(t,j,s):
        if t == n:
            return 0
        else:
            return (call_value(t+1,j,s)-call_value(t,j,s))/dt

    def vega(t,j,s):
        if t == n:
            return 0
        else:
            return np.exp(-r*dt)*((p*vega(t+1,j+1,s)+(1-p)*vega(t+1,j,s))/(u-d)-(call_value(t+1,j+1,s)-call_value(t+1,j-1,s))/(2*sigma*np.sqrt(dt)*s*u**(n-t-1)*d**(j-1)))

    def rho(t,j,s):
        if t == n:
            return 0
        else:
            return (call_value(t+1,j,s)-call_value(t,j,s))/dt

    delta_0 = delta(0,0,s)
    gamma_0 = delta(0,0,s)
    theta_0 = theta(0,0,s)
    vega_0 = vega(0,0,s)
    rho_0 = rho(0,0,s)
    return 'delta =',delta_0,'gamma =',gamma_0,'theta =',theta_0,'vega =',vega_0,'rho =',rho_0

def us_put_greeks(
        s: float, k: float, r: float, t: float, sigma: float, n: float
) -> float:

    dt = t/n
    u = np.exp(sigma*np.sqrt(dt))
    d = 1/u
    p = (np.exp(r*dt)-d)/(u-d)

    def put_value(t,j,s): 
        if t == n:
            return max(s-k,0)
        else:
            previous_exercise = max(k-s*u**j*d**(n-j),0)
            current_exercise = np.exp(-r*dt)*(p*put_value(t+1,j+1,s)+(1-p)*put_value(t+1,j,s))
            return max(previous_exercise,current_exercise)

    def delta(t,j,s):
        if t == n:
            return max(k-u*s*d**(n-1-j),0) - max(k-d*s*u**(n-1-j),0)
        else:
            return np.exp(-r*dt)*(p*delta(t+1,j+1,s)+(1-p)*delta(t+1,j,s))/(u-d)

    def gamma(t,j,s):
        if t == n:
            return 0
        else:
            return np.exp(-r*dt)*((p*gamma(t+1,j+1,s)+(1-p)*gamma(t+1,j,s))/(u-d)-delta(t+1,j+1,s)/(s*u**(n-t-1)*d**j-s*u**(n-t-1)*d**(j+1)))

    def theta(t,j,s):
        if t == n:
            return 0
        else:
            return (put_value(t+1,j,s)-put_value(t,j,s))/dt

    def vega(t,j,s):
        if t == n:
            return 0
        else:
            return np.exp(-r*dt)*((p*vega(t+1,j+1,s)+(1-p)*vega(t+1,j,s))/(u-d)-(put_value(t+1,j+1,s)-put_value(t+1,j-1,s))/(2*sigma*np.sqrt(dt)*s*u**(n-t-1)*d**(j-1)))

    def rho(t,j,s):
        if t == n:
            return 0
        else:
            return (put_value(t+1,j,s)-put_value(t,j,s))/dt

    delta_0 = delta(0,0,s)
    gamma_0 = delta(0,0,s)
    theta_0 = theta(0,0,s)
    vega_0 = vega(0,0,s)
    rho_0 = rho(0,0,s)
    return 'delta =',delta_0,'gamma =',gamma_0,'theta =',theta_0,'vega =',vega_0,'rho =',rho_0