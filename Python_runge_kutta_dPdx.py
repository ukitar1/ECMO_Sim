import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


# Q_b = flow cm3/s
# d is diameter in cm
# p is porosity, unitless
# SA is surface area in m2
# _Hgb is hemoglobin, g/dL
# _Hct is hematocrit in fraction
# _T is temperature in Celsius

def calc_VO2(P0, Q_b, d, p, A_f, L, SA, _Hgb, _Hct, _T, alpha, beta):


    def solubility(T = _T, Hct = _Hct):
        exp = math.pow(1.01, 37-T)
        
        return Hct*4.658E-5*exp+(1-Hct)*2.855E-5*exp
    
    _k = solubility()
    
    def diffusivity(T = _T, Hct = _Hct):
        M = 0.283
        exp_k = math.pow(1.01, 37-T)
        exp_D = math.pow(1.025, T-25)
        
        k_pl = 2.855E-5*exp_k
        k_cell = 4.658E-5*exp_k
        
        D_pl = 1.62E-5*exp_D
        D_cell = .76E-5*exp_D
        
        N = k_cell*D_cell/(k_pl*D_pl)
    
        _B = (N-1)/3*(2/(1+M/2*(N-1))+1/(1+(N-1)*(M-1)))
        _X = (1-N*(1-_B))/(N-1-_B)
        _G = Hct*(N-1)/(N+_X)
        
        return D_pl*k_pl*(1+_X*_G)/(_k*(1-_G))
        
    _D = diffusivity()
    
    PO2_g = 760*.95
    
    mu = 0.03 # dyn/cm2*s
    rho = 1 #g/cm3
    nu =  mu/rho #cm2/s
    pH = 7.4
    
    #bovine
    _P50 = 29*math.pow(10, 0.41*(7.4-pH))
    
    #human
    #_P50 = 26.6*math.pow(10, 0.48*(7.4-pH))
    
    
    #cTAL, Keith
    
    # Dave's thesis

    
    def Hill(P, P50 = _P50, n = 2.85): #2.85 for bovine
        X = math.pow(P/P50, n)
        return X/(1 + X)
    
    
    #lambda(P)
    def sink(P, Hgb = _Hgb):
        dSdP = (Hill(P+.05)-Hill(P-.05))/.1
        return 1.34*Hgb/_k*dSdP * 0.01 #dl blood --> ml blood
    
    
    # Python program to implement Runge Kutta method
    def dPdx(x, P, D = _D, P_g = PO2_g):
        term1 = 4*alpha/p*math.pow((1-p)/d, 2-beta)
        term2 = math.pow(A_f*nu/Q_b, 1-beta)
        term3 = math.pow(D/nu, 2/3.)
        term4 = (P_g - P)/math.pow((1+sink(P)), 2/3.)
        
        
        return term1*term2*term3*term4
     
    # Finds value of y for a given x using step size h
    # and initial value y0 at x0.
    def rungeKutta(x0, P0, x, h):
        # Count number of iterations using step size or
        # step height h
        n = (int)((x - x0)/h)
        # Iterate for number of iterations
        P = P0
        for i in range(1, n + 1):
            "Apply Runge Kutta Formulas to find next value of y"
            k1 = h * dPdx(x0, P0)
            k2 = h * dPdx(x0 + 0.5 * h, P + 0.5 * k1)
            k3 = h * dPdx(x0 + 0.5 * h, P + 0.5 * k2)
            k4 = h * dPdx(x0 + h, P + k3)
     
            # Update next value of y
            P = P + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
     
            # Update next value of x
            x0 = x0 + h
        return P
     
    
    # initial conditions: 
    x0 = 0 
    
    h= 0.005

    P_out = rungeKutta(x0 = x0, P0 = P0, x = L, h = h) 
    
    S_out = Hill(P_out); 
    
    VO2= Q_b*60/100*(1.34*_Hgb*(Hill(P_out)-Hill(P0)) + _k*(P_out - P0))
    
    return [S_out, VO2]


Q_list = [1, 2, 3, 4, 5, 6] #LPM

#CONSTANTS, SMO     
dd = 290 #microns
dd = dd *1E-4 #convert to cm
pp = 0.576 #porosity
Area = 36.88 #cm2  
PL = 9.5 #cm
S_A = 1.87 #m^2
alpha_SMO = 0.136
beta_SMO = 0.832
  

#CONSTANTS, Univox
# dd = 365 #microns
# dd = dd *1E-4 #convert to cm
# pp = 0.31 #porosity
# Area = 207.7 #cm2  
# PL = 1.2 #cm
# S_A = 1.8 #m^2

# alpha_vox = 0.256 #PHI
# beta_vox = 0.779 #m

Hb = 13 #g/dL
crit = 0.37 #
tt = 37 #C

#Implantable Artificial Lung 1994
dd = 380 #microns
dd = dd *1E-4 #convert to cm
pp = 0.53 #porosity
Area = 128  #cm2, frontal
PL = 3.5 #cm
S_A = 2.2 #m^2
alpha_TAL = 0.363
beta_TAL = 0.725


mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 30
plt.rcParams['axes.linewidth'] = 3

fig, ax = plt.subplots(figsize = (12, 15), linewidth = 2, edgecolor = 'k')

P_in = 25

#SMO 
# for Q in Q_list:
    
#     Q = Q * 1000/60.
#     [S1,V1] = calc_VO2(25, Q, dd, pp, Area, PL, S_A, Hb, crit, tt, alpha_SMO, beta_SMO) #38%
#     [S2,V2] = calc_VO2(33, Q, dd, pp, Area, PL, S_A, Hb, crit, tt, alpha_SMO, beta_SMO) #59%
#     [S3,V3] = calc_VO2(50, Q, dd, pp, Area, PL, S_A, Hb, crit, tt, alpha_SMO, beta_SMO) #83%

#     Q_LPM = Q *60/1000
#     ax.scatter(Q_LPM, V1, s = 250, color = 'b')
#     ax.scatter(Q_LPM, V2, s = 250, color = 'k')
#     ax.scatter(Q_LPM, V3, s = 250, color = 'r')


# ax.set_xlim([0, 6.2])
# ax.set_ylim([0, 450])
    
# ax.set_title('Univox')
# ax.set_xlabel('Flow (L/min)')
# ax.set_ylabel('Oxygen Transfer (mL/min)')

Q_TAL = [5]
Hb = 11
crit = 0.33

for Q in Q_TAL:
    
    Q = Q * 1000/60
    #[S1,V1] = calc_VO2(41.1, Q, dd, pp, Area, PL, S_A, 12.3, 0.37, tt, alpha_TAL, beta_TAL) #Experiment 1
    [S2,V2] = calc_VO2(35.4, Q, dd, pp, Area, PL, S_A, 10.7, 0.37, tt, alpha_TAL, beta_TAL) #Experiment 2
    [S3,V3] = calc_VO2(35.4, Q, dd, pp, Area, PL, S_A, 11, 0.37, tt, alpha_TAL, beta_TAL) #Experiment 3
    [S4,V4] = calc_VO2(32.5, Q, dd, pp, Area, PL, S_A, 8, 0.30, tt, alpha_TAL, beta_TAL) #Experiment 4
    [S5,V5] = calc_VO2(29.8, Q, dd, pp, Area, PL, S_A, 8.2, 0.31, tt, alpha_TAL, beta_TAL) #Experiment5

    Q_LPM = Q *60/1000
   # ax.scatter(Q_LPM, V1, s = 250, color = 'b')
    ax.scatter(Q_LPM, V2, s = 250, color = 'k')
    ax.scatter(Q_LPM, V3, s = 250, color = 'r')
    ax.scatter(Q_LPM, V4, s = 250, color = 'g')
    ax.scatter(Q_LPM, V5, s = 250, color = 'y')
    
ax.set_xlim([0, 5.2])
ax.set_ylim([0, 250])
    
ax.set_title('TAL')
ax.set_xlabel('Flow (L/min)')
ax.set_ylabel('Oxygen Transfer (mL/min)')
