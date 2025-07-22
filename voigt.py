import numpy as np

def gaussian(x, sigma, mu):
    pwr = -np.power((x-mu)/sigma, 2)/2
    y = np.exp(pwr)/(sigma*np.sqrt(2*np.pi))
    return y

def lorentzian(x, gamma, mu):
    y = gamma/(np.power(x-mu, 2)+np.power(gamma, 2))*1/np.pi
    return y

def voigt(x, fL, fG, mu):
    f = np.power(np.power(fG, 5) + 2.69269*np.power(fG, 4)*fL + 2.42843*np.power(fG, 3)*np.power(fL, 2)
                 + 4.47163*np.power(fG, 2)*np.power(fL, 3) + 0.07842*fG*np.power(fL, 4) + np.power(fL, 5), 1/5)
    eta = 1.36603*(fL/f) - 0.47719*np.power(fL/f, 2) + 0.11116*np.power(fL/f, 3)
    y = eta*lorentzian(x, fL, mu) + (1-eta)*gaussian(x, fG, mu)
    return y