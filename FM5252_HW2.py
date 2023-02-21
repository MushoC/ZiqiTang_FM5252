#    S (float): current price
#    K (float): strike price
#    r (float): risk free rate
#    sigma (float): underlying volatility
#    T (float): remaining time to expiry
#    N (float): number of steps that should be generated
#    delta_r (float): random little number that customer consider as floating volatility

import numpy as np

# European call option, option price and Greek values
def crr_eu_call(
    S: float, K: float, r: float, T: float, sigma: float, N: float, delta_r: float
) -> float:

    dt = T/N
    u = np.exp(sigma*np.sqrt(dt))
    d = 1/u
    p = (np.exp(r*dt) - d) / (u - d)

    # Recursive function to calculate option value at each node
    def option_value(n, j):
        # 'n' is current time step, 'j' is the index of the node, both of them default to 0
        if n == N:
            return max(S * u**j * d**(N-j) - K, 0)
        else:
            return np.exp(-r*dt) * (p * option_value(n+1, j+1) + (1-p) * option_value(n+1, j))

    # Recursive function to calculate delta
    def delta(n, j):
        if n == N:
            return max(S*u*d**(N-1-j) - K, 0) - max(S*d*u**(N-1-j) - K, 0)
        else:
            return np.exp(-r*dt) * (p * delta(n+1, j+1) + (1-p) * delta(n+1, j)) /(S*u**(N-n-1)*d**j - S*u**(N-n-1)*d**(j+1))

    # Recursive function to calculate gamma
    def gamma(n, j):
        if n == N:
            return 0
        else:
            return 2*(delta(n+1,j+1)-delta(n+1,j))/((S*u**(N-n-1)*d**j - S*u**(N-n-1)*d**(j+2)))
        
    # Recursive function to calculate theta
    def theta(n, j):
        if n == N:
            return 0
        else:
            return (option_value(n+2, j) - option_value(n, j)) / (2*T)
        
    # Recursive function to calculate vega
    def vega(n, j):
        if n == N:
            return 0
        else:
            u1 = np.exp((sigma+delta_r)*np.sqrt(dt))
            d1 = 1/u1

            u2 = np.exp((sigma-delta_r)*np.sqrt(dt))
            d2 = 1/u2

            def option_value_vega_plus(n, j):
              if n == N:
                return max(S * u1**j * d1**(N-j) - K, 0)
              else:
                return np.exp(-r*dt) * (p * option_value_vega_plus(n+1, j+1) + (1-p) * option_value_vega_plus(n+1, j))

            def option_value_vega_minus(n, j):
              if n == N:
                return max(S * u2**j * d2**(N-j) - K, 0)
              else:
                return np.exp(-r*dt) * (p * option_value_vega_minus(n+1, j+1) + (1-p) * option_value_vega_minus(n+1, j))

            return (option_value_vega_plus(n,j)-option_value_vega_minus(n,j)) / (2*delta_r)
        
    # Recursive function to calculate rho
    def rho(n, j):
        if n == N:
            return 0
        else:
            p1 = (np.exp((r+delta_r)*dt) - d) / (u - d)
            p2 = (np.exp((r-delta_r)*dt) - d) / (u - d)

            def option_value_rho_plus(n, j):
              if n == N:
                return max(S * u**j * d**(N-j) - K, 0)
              else:
                return np.exp(-(r+delta_r)*dt) * (p1 * option_value_rho_plus(n+1, j+1) + (1-p1) * option_value_rho_plus(n+1, j))

            def option_value_rho_minus(n, j):
              if n == N:
                return max(S * u**j * d**(N-j) - K, 0)
              else:
                return np.exp(-(r-delta_r)*dt) * (p2 * option_value_rho_minus(n+1, j+1) + (1-p2) * option_value_rho_minus(n+1, j))

            return (option_value_rho_plus(n, j) - option_value_rho_minus(n, j)) / (2*delta_r)

    option_value_0 = option_value(0, 0)
    delta_0 = delta(0, 0)
    gamma_0 = gamma(0, 0)
    theta_0 = theta(0, 0)
    vega_0 = vega(0, 0)
    rho_0 = rho(0, 0)

    return 'option value =',option_value_0,'delta =',delta_0,'gamma =',gamma_0,'theta =',theta_0,'vega =',vega_0,'rho =',rho_0



# European put option, option price and Greek values
def crr_eu_put(
    S: float, K: float, r: float, T: float, sigma: float, N: float, delta_r: float
) -> float:

    dt = T/N
    u = np.exp(sigma*np.sqrt(dt))
    d = 1/u
    p = (np.exp(r*dt) - d) / (u - d)

    def option_value(n, j):
        if n == N:
            return max(K - S * u**j * d**(N-j), 0)
        else:
            return np.exp(-r*dt) * (p * option_value(n+1, j+1) + (1-p) * option_value(n+1, j))

    # Recursive function to calculate delta
    def delta(n, j):
        if n == N:
            return max(K - S*u*d**(N-1-j), 0) - max(K - S*d*u**(N-1-j), 0)
        else:
            return np.exp(-r*dt) * (p * delta(n+1, j+1) + (1-p) * delta(n+1, j)) /(S*u**(N-n-1)*d**j - S*u**(N-n-1)*d**(j+1))

    # Recursive function to calculate gamma
    def gamma(n, j):
        if n == N:
            return 0
        else:
            return 2*(delta(n+1,j+1)-delta(n+1,j))/((S*u**(N-n-1)*d**j - S*u**(N-n-1)*d**(j+2)))
        
    # Recursive function to calculate theta
    def theta(n, j):
        if n == N:
            return 0
        else:
            return (option_value(n+2, j) - option_value(n, j)) / (2*T)
        
    # Recursive function to calculate vega
    def vega(n, j):
        if n == N:
            return 0
        else:
            u1 = np.exp((sigma+delta_r)*np.sqrt(dt))
            d1 = 1/u1

            u2 = np.exp((sigma-delta_r)*np.sqrt(dt))
            d2 = 1/u2

            def option_value_vega_plus(n, j):
              if n == N:
                return max(K - S * u1**j * d1**(N-j), 0)
              else:
                return np.exp(-r*dt) * (p * option_value_vega_plus(n+1, j+1) + (1-p) * option_value_vega_plus(n+1, j))

            def option_value_vega_minus(n, j):
              if n == N:
                return max(K - S * u2**j * d2**(N-j), 0)
              else:
                return np.exp(-r*dt) * (p * option_value_vega_minus(n+1, j+1) + (1-p) * option_value_vega_minus(n+1, j))

            return (option_value_vega_plus(n,j)-option_value_vega_minus(n,j)) / (2*delta_r)
        
    # Recursive function to calculate rho
    def rho(n, j):
        if n == N:
            return 0
        else:
            p1 = (np.exp((r+delta_r)*dt) - d) / (u - d)
            p2 = (np.exp((r-delta_r)*dt) - d) / (u - d)

            def option_value_rho_plus(n, j):
              if n == N:
                return max(K - S * u**j * d**(N-j), 0)
              else:
                return np.exp(-(r+delta_r)*dt) * (p1 * option_value_rho_plus(n+1, j+1) + (1-p1) * option_value_rho_plus(n+1, j))

            def option_value_rho_minus(n, j):
              if n == N:
                return max(K - S * u**j * d**(N-j), 0)
              else:
                return np.exp(-(r-delta_r)*dt) * (p2 * option_value_rho_minus(n+1, j+1) + (1-p2) * option_value_rho_minus(n+1, j))

            return (option_value_rho_plus(n, j) - option_value_rho_minus(n, j)) / (2*delta_r)

    option_value_0 = option_value(0, 0)
    delta_0 = delta(0, 0)
    gamma_0 = gamma(0, 0)
    theta_0 = theta(0, 0)
    vega_0 = vega(0, 0)
    rho_0 = rho(0, 0)

    return 'option value =',option_value_0,'delta =',delta_0,'gamma =',gamma_0,'theta =',theta_0,'vega =',vega_0,'rho =',rho_0

# American call option, option price and Greek values
def crr_us_call(
    S: float, K: float, r: float, T: float, sigma: float, N: float, delta_r: float
) -> float:

    dt = T/N
    u = np.exp(sigma*np.sqrt(dt))
    d = 1/u
    p = (np.exp(r*dt) - d) / (u - d)

    def option_value(n, j):
        if n == N:
            return max(S * u**j * d**(N-j) - K, 0)

        else:
            exercise_value = S * u**j * d**(n-j) - K
            continuation_value = np.exp(-r*dt) * (p * option_value(n+1, j+1) + (1-p) * option_value(n+1, j))

        return max(exercise_value, continuation_value)

    # Recursive function to calculate delta
    def delta(n, j):
        if n == N:
            return max(S*u*d**(N-1-j) - K, 0) - max(S*d*u**(N-1-j) - K, 0)
        else:
            return np.exp(-r*dt) * (p * delta(n+1, j+1) + (1-p) * delta(n+1, j)) /(S*u**(N-n-1)*d**j - S*u**(N-n-1)*d**(j+1))

    # Recursive function to calculate gamma
    def gamma(n, j):
        if n == N:
            return 0
        else:
            return 2*(delta(n+1,j+1)-delta(n+1,j))/((S*u**(N-n-1)*d**j - S*u**(N-n-1)*d**(j+2)))
        
    # Recursive function to calculate theta
    def theta(n, j):
        if n == N:
            return 0
        else:
            return (option_value(n+2, j) - option_value(n, j)) / (2*T)
        
    # Recursive function to calculate vega
    def vega(n, j):
        if n == N:
            return 0
        else:
            u1 = np.exp((sigma+delta_r)*np.sqrt(dt))
            d1 = 1/u1

            u2 = np.exp((sigma-delta_r)*np.sqrt(dt))
            d2 = 1/u2

            def option_value_vega_plus(n, j):
              if n == N:
                return max(S * u1**j * d1**(N-j) - K, 0)
              else:
                exercise_value = S * u1**j * d1**(n-j) - K
                continuation_value = np.exp(-r*dt) * (p * option_value_vega_plus(n+1, j+1) + (1-p) * option_value_vega_plus(n+1, j))
                return max(exercise_value, continuation_value)
            

            def option_value_vega_minus(n, j):
              if n == N:
                return max(S * u2**j * d2**(N-j) - K, 0)
              else:
                exercise_value = S * u2**j * d2**(n-j) - K
                continuation_value = np.exp(-r*dt) * (p * option_value_vega_minus(n+1, j+1) + (1-p) * option_value_vega_minus(n+1, j))

                return max(exercise_value, continuation_value)

            return (option_value_vega_plus(n,j)-option_value_vega_minus(n,j)) / (2*delta_r)
        
    # Recursive function to calculate rho
    def rho(n, j):
        if n == N:
            return 0
        else:
            p1 = (np.exp((r+delta_r)*dt) - d) / (u - d)
            p2 = (np.exp((r-delta_r)*dt) - d) / (u - d)

            def option_value_rho_plus(n, j):
              if n == N:
                return max(S * u**j * d**(N-j), 0) - K
              else:
                exercise_value = S * u**j * d**(n-j) - K
                continuation_value = np.exp(-r*dt) * (p1 * option_value(n+1, j+1) + (1-p1) * option_value(n+1, j))

                return max(exercise_value, continuation_value)

            def option_value_rho_minus(n, j):
              if n == N:
                return max(S * u**j * d**(N-j), 0) - K
              else:
                exercise_value = S * u**j * d**(n-j) - K
                continuation_value = np.exp(-r*dt) * (p2 * option_value(n+1, j+1) + (1-p2) * option_value(n+1, j))

                return max(exercise_value, continuation_value)
                
            return (option_value_rho_plus(n, j) - option_value_rho_minus(n, j)) / (2*delta_r)

    option_value_0 = option_value(0, 0)
    delta_0 = delta(0, 0)
    gamma_0 = gamma(0, 0)
    theta_0 = theta(0, 0)
    vega_0 = vega(0, 0)
    rho_0 = rho(0, 0)

    return 'option value =',option_value_0,'delta =',delta_0,'gamma =',gamma_0,'theta =',theta_0,'vega =',vega_0,'rho =',rho_0

# American put option, option price and Greek values
def crr_us_put(
    S: float, K: float, r: float, T: float, sigma: float, N: float, delta_r: float
) -> float:

    dt = T/N
    u = np.exp(sigma*np.sqrt(dt))
    d = 1/u
    p = (np.exp(r*dt) - d) / (u - d)

    def option_value(n, j):
        if n == N:
            return max(K - S * u**j * d**(N-j), 0)

        else:
            exercise_value = K - S * u**j * d**(n-j)
            continuation_value = np.exp(-r*dt) * (p * option_value(n+1, j+1) + (1-p) * option_value(n+1, j))

        return max(exercise_value, continuation_value)

    # Recursive function to calculate delta
    def delta(n, j):
        if n == N:
            return max(K - S*u*d**(N-1-j), 0) - max(K - S*d*u**(N-1-j), 0)
        else:
            return np.exp(-r*dt) * (p * delta(n+1, j+1) + (1-p) * delta(n+1, j)) /(S*u**(N-n-1)*d**j - S*u**(N-n-1)*d**(j+1))

    # Recursive function to calculate gamma
    def gamma(n, j):
        if n == N:
            return 0
        else:
            return 2*(delta(n+1,j+1)-delta(n+1,j))/((S*u**(N-n-1)*d**j - S*u**(N-n-1)*d**(j+2)))
        
    # Recursive function to calculate theta
    def theta(n, j):
        if n == N:
            return 0
        else:
            return (option_value(n+2, j) - option_value(n, j)) / (2*T)
        
    # Recursive function to calculate vega
    def vega(n, j):
        if n == N:
            return 0
        else:
            u1 = np.exp((sigma+delta_r)*np.sqrt(dt))
            d1 = 1/u1

            u2 = np.exp((sigma-delta_r)*np.sqrt(dt))
            d2 = 1/u2

            def option_value_vega_plus(n, j):
              if n == N:
                return max(K - S * u1**j * d1**(N-j), 0)
              else:
                exercise_value = K - S * u1**j * d1**(n-j)
                continuation_value = np.exp(-r*dt) * (p * option_value_vega_plus(n+1, j+1) + (1-p) * option_value_vega_plus(n+1, j))
                return max(exercise_value, continuation_value)
            

            def option_value_vega_minus(n, j):
              if n == N:
                return max(K - S * u2**j * d2**(N-j), 0)
              else:
                exercise_value = K - S * u2**j * d2**(n-j)
                continuation_value = np.exp(-r*dt) * (p * option_value_vega_minus(n+1, j+1) + (1-p) * option_value_vega_minus(n+1, j))

                return max(exercise_value, continuation_value)

            return (option_value_vega_plus(n,j)-option_value_vega_minus(n,j)) / (2*delta_r)
        
    # Recursive function to calculate rho
    def rho(n, j):
        if n == N:
            return 0
        else:
            p1 = (np.exp((r+delta_r)*dt) - d) / (u - d)
            p2 = (np.exp((r-delta_r)*dt) - d) / (u - d)

            def option_value_rho_plus(n, j):
              if n == N:
                return max(K - S * u**j * d**(N-j), 0)
              else:
                exercise_value = K - S * u**j * d**(n-j)
                continuation_value = np.exp(-r*dt) * (p1 * option_value(n+1, j+1) + (1-p1) * option_value(n+1, j))

                return max(exercise_value, continuation_value)

            def option_value_rho_minus(n, j):
              if n == N:
                return max(K - S * u**j * d**(N-j), 0)
              else:
                exercise_value = K - S * u**j * d**(n-j)
                continuation_value = np.exp(-r*dt) * (p2 * option_value(n+1, j+1) + (1-p2) * option_value(n+1, j))

                return max(exercise_value, continuation_value)

            return (option_value_rho_plus(n, j) - option_value_rho_minus(n, j)) / (2*delta_r)

    option_value_0 = option_value(0, 0)
    delta_0 = delta(0, 0)
    gamma_0 = gamma(0, 0)
    theta_0 = theta(0, 0)
    vega_0 = vega(0, 0)
    rho_0 = rho(0, 0)

    return 'option value =',option_value_0,'delta =',delta_0,'gamma =',gamma_0,'theta =',theta_0,'vega =',vega_0,'rho =',rho_0
