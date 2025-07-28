import time
import yaml
import numpy as np
import matplotlib.pyplot as plt
from importhitran import importhitran
from importhitran import importpartfun
from voigt import voigt

print('Executing...')
start_time = time.time()

# Parameter setting and data loading
try:
    with open('config.yaml', 'r') as config_file:
        data = yaml.safe_load(config_file)
        molecules = data["molecules"]
        parameters = data["parameters"]
except FileNotFoundError:
    print('Config file not found')

T = parameters[0]["temperature"]
L = parameters[1]["length"]
PTorr = parameters[2]["PTorr"]
nu0 = parameters[3]["nu0"]
nuf = parameters[4]["nuf"]

vv = np.linspace(nu0, nuf, 100000)
Abs = np.zeros(vv.shape)

for molecule in molecules:

    #molecule = molecules[1]
    molecule_file = molecule[1]["molecule_file"]
    partfun_file = molecule[2]["partfun_file"]
    M = molecule[3]["mass"]
    PTorr_M = molecule[4]["PTorr_M"]
    nd = (PTorr_M/760) * 0.101325 / (1.38 * (10**-23) * T)

    lines = importhitran(molecule_file)
    Ind = (lines.transitionWavenumber > nu0) & (lines.transitionWavenumber < nuf)
    Inten = lines.lineIntensity[Ind]
    E_j_L = lines.lowerStateEnergy[Ind]
    WaveN = lines.transitionWavenumber[Ind]
    gamma_P =np.power(lines.airBroadenedWidth[Ind] * 2 * (296/T), 0.6)

    #Line Intensity temperature Correction
    Tref = 296
    T_part, Q_part = importpartfun(partfun_file)
    QTref = Q_part[np.argmin(abs(T_part - Tref))]
    QT = Q_part[np.argmin(abs(T_part - T))]

    rQ = QTref/QT
    h = 6.626*10**-34
    c0 = 2.9979245800*10**10
    kb = 1.38*10**-23
    c2 = h*c0/kb
    rB = np.divide(np.exp(-c2*E_j_L/T), np.exp(-c2*E_j_L/Tref))
    rE = np.divide(1 - np.exp(-c2*WaveN/T), 1 - np.exp(-c2*WaveN/T))
    Inten = rQ*rB*rE*Inten

    #Generate Synthetic Spectrum
    P = PTorr/760
    FWHM_L = gamma_P*P
    Ind = np.arange(0, np.size(WaveN), 1) #Figure out what the fuck this is
    spec = np.zeros(vv.shape)
    for i in Ind:
        v0 = WaveN[i]
        FWHM_D = (7.17**-7)*v0*np.sqrt(T/M)
        V = voigt(vv, FWHM_L[i], FWHM_D, v0)
        spec = spec + V*Inten[i]
    Abs = Abs + spec*nd*L
        
end_time = time.time()
elapsed_time = end_time - start_time
print(f'Elapsed time: {elapsed_time:.2f} seconds')

plt.plot(vv, Abs)
plt.xlim([nu0, nuf])
plt.ylim([0, max(Abs) + 0.1*max(Abs)])
plt.xlabel('Wavenumber (cm^-1)')
plt.ylabel(r'$\alpha$L')
plt.show()
