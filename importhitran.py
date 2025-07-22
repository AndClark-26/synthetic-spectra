import sys
import numpy as np

class molecule:
    def __init__(self, moleculeNumber, isotopologueNumber, transitionWavenumber, lineIntensity, einsteinACoefficient, airBroadenedWidth, selfBroadenedWidth, 
                 lowerStateEnergy, temperatureDependance, pressureShift, upperVibrationalQuanta, lowerVibrationalQuanta, upperLocalQuanta, lowerLocalQuanta):
        self.moleculeNumber = moleculeNumber
        self.isotopologueNumber = isotopologueNumber
        self.transitionWavenumber = transitionWavenumber
        self.lineIntensity = lineIntensity
        self.einsteinACoefficient = einsteinACoefficient
        self.airBroadenedWidth = airBroadenedWidth
        self.selfBroadenedWidth = selfBroadenedWidth
        self.lowerStateEnergy = lowerStateEnergy
        self.temperatureDependance = temperatureDependance
        self.pressureShift = pressureShift
        self.upperVibrationalQuanta = upperVibrationalQuanta
        self.lowerVibrationalQuanta = lowerVibrationalQuanta
        self.upperLocalQuanta = upperLocalQuanta
        self.lowerLocalQuanta = lowerLocalQuanta

def importhitran(file):
    try:
          with open(file) as f:
                    lines = f.readlines()
                    line_count = len(lines)

                    moleculeNumber = np.empty((line_count, 1), dtype=int)
                    isotopologueNumber = np.empty((line_count, 1), dtype=int)
                    transitionWavenumber = np.empty((line_count, 1), dtype=float)
                    lineIntensity = np.empty((line_count, 1), dtype=float)
                    einsteinACoefficient = np.empty((line_count, 1), dtype=float)
                    airBroadenedWidth = np.empty((line_count, 1), dtype=float)
                    selfBroadenedWidth = np.empty((line_count, 1), dtype=float)
                    lowerStateEnergy = np.empty((line_count, 1), dtype=float)
                    temperatureDependance = np.empty((line_count, 1), dtype=float)
                    pressureShift = np.empty((line_count, 1), dtype=float)
                    upperVibrationalQuanta = np.empty((line_count, 1), dtype=str)
                    lowerVibrationalQuanta = np.empty((line_count, 1), dtype=str)
                    upperLocalQuanta = np.empty((line_count, 1), dtype=str)
                    lowerLocalQuanta = np.empty((line_count, 1), dtype=str)

                    i = 0
                    for line in lines:
                         moleculeNumber[i] = int(line[0:2])
                         isotopologueNumber[i] = int(line[2])
                         transitionWavenumber[i] = float(line[3:15])
                         lineIntensity[i] = float(line[15:25])
                         einsteinACoefficient[i] = float(line[25:35])
                         airBroadenedWidth[i] = float(line[35:40])
                         selfBroadenedWidth[i] = float(line[40:45])
                         lowerStateEnergy[i] = float(line[45:55])
                         temperatureDependance[i] = float(line[55:59])
                         pressureShift[i] = float(line[59:67])
                         upperVibrationalQuanta[i] = line[67:82]
                         lowerVibrationalQuanta[i] = line[82:97]
                         upperLocalQuanta[i] = line[97:112]
                         lowerLocalQuanta[i] = line[112:127]
                         i = i + 1

                    moleculeData = molecule(moleculeNumber, isotopologueNumber, transitionWavenumber, lineIntensity, einsteinACoefficient, airBroadenedWidth, selfBroadenedWidth, 
                                        lowerStateEnergy, temperatureDependance, pressureShift, upperVibrationalQuanta, lowerVibrationalQuanta, upperLocalQuanta, lowerLocalQuanta)
     
                    return moleculeData
    except FileNotFoundError:
         print('Hitran file not found')
         sys.exit(1)

def importpartfun(file):
     try:
          with open(file) as f:
               lines = f.readlines()
               line_count = len(lines)

               T_part = np.empty((line_count - 4, 1), dtype=int)
               Q_part = np.empty((line_count - 4, 1), dtype=float)
               
               i = 0
               for line in lines:
                    if i < 4:
                         i = i + 1
                         continue
                    else:
                         vals = line.split()
                         T_part[i - 4] = int(vals[0])
                         Q_part[i - 4] = float(vals[1])
                         i = i + 1
          return T_part, Q_part
     except FileNotFoundError:
          print('Partiton file function not found')
          sys.exit(1)