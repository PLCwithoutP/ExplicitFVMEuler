import math
import numpy as np
import postProcesser as pp

gamma = 1.4
beta = 1/8
K_u = 0.75
K_p = 0.25
sigma = 1 

def defineInputValues(cellValues):
    rho = cellValues[0]
    u_x = cellValues[1] / cellValues[0]
    u_y = cellValues[2] / cellValues[0]
    u_z = cellValues[3] / cellValues[0]
    E = cellValues[4]
    return rho, u_x, u_y, u_z, E

def calculateSpeedOfSound(rho, p):
    p = max(0,p)
    a = math.sqrt(gamma*p/rho)
    return a

def calculateFaceSpeedOfSound(a_L, a_R):
    a_f = (a_L + a_R)/2
    return a_f

def calculateFaceDensity(rho_L, rho_R):
    rho_12 = (rho_L + rho_R)/2
    return rho_12

def calculateLeftRightMach(n, u_L, u_R, a_12):
    M_L = np.dot(n,u_L)/a_12
    M_R = np.dot(n,u_R)/a_12
    return M_L, M_R

def calculateMeanMach(u_L, u_R, a_12):
    meanMach = math.sqrt((np.linalg.norm(u_L)**2 + np.linalg.norm(u_R)**2)/(2*a_12**2))
    return meanMach

def calculateM0(meanMach, inletMach):
    M_0 = math.sqrt(min(1,max(meanMach**2, inletMach**2)))
    return M_0

def calculateScalingFactor(M_0):
    f_a = M_0 * (2 - M_0)
    return f_a

def mach4PlusFunction(mach):
    if (abs(mach) >= 1):
        result = 0.5*(mach + abs(mach))
    else:
        result = 0.25*((mach + 1)**2)*(1 + 16*beta*0.25*(mach - 1)**2)
    return result

def mach4MinusFunction(mach):
    if (abs(mach) >= 1):
        result = 0.5*(mach - abs(mach))
    else:
        result = -0.25*((mach - 1)**2)*(1 + 16*beta*0.25*(mach + 1)**2)
    return result

def calculateFaceMach(M_L, M_R, M_0, meanMach, f_a, p_R, p_L, rho_12, a_12):
    M_12 = mach4PlusFunction(M_L) + mach4MinusFunction(M_R) - ((K_p/f_a) * max(1 - sigma*meanMach*meanMach , 0)*(p_R-p_L)/(rho_12*a_12**2))
    return M_12

def calculateMassFlux(M_12, a_12, rho_L, rho_R):
    if (M_12 > 0):
        result = a_12*M_12*rho_L
    else:
        result = a_12*M_12*rho_R
    return result

def pressure5PlusFunction(mach, f_a):
    alpha = 3/16*(-4 + 5*f_a**2)
    if (abs(mach) >= 1):
        result = (1/(2*mach))*(mach + abs(mach))
    else:
        result = 0.25*((mach + 1)**2)*((2 - mach) + 16*alpha*mach*0.25*(mach - 1)**2)
    return result

def pressure5MinusFunction(mach, f_a):
    alpha = 3/16*(-4 + 5*f_a**2)
    if (abs(mach) >= 1):
        result = (1/(2*mach))*(mach - abs(mach))
    else:
        result = -0.25*((mach - 1)**2)*((-2 - mach) + 16*alpha*mach*0.25*(mach + 1)**2)
    return result

def calculateFacePressure(M_L, M_R, rho_L, rho_R, p_L, p_R, f_a, a_12):
    result = pressure5PlusFunction(M_L, f_a)*p_L + pressure5MinusFunction(M_R, f_a)*p_R - (K_u*pressure5PlusFunction(M_L, f_a)*pressure5MinusFunction(M_R, f_a)*(rho_L + rho_R)*(f_a*a_12**2)*(M_R - M_L))
    return result

def calculateEnthalpy(E, rho, p):
    H = (E + p) / rho
    return H

def decideFinalFlux(M_12, mdot_12, p_12, u_L, u_R, H_L, H_R, n):
    flux1 = [0,0,0,0,0]
    flux2 = [0,0,0,0,0]
    flux_total = [0,0,0,0,0]
    flux2[0] = 0
    flux2[1] = p_12*n[0]
    flux2[2] = p_12*n[1]
    flux2[3] = p_12*n[2]
    flux2[4] = 0
    if (M_12 > 0):
        flux1[0] = mdot_12
        flux1[1] = mdot_12*u_L[0]
        flux1[2] = mdot_12*u_L[1]
        flux1[3] = mdot_12*u_L[2]
        flux1[4] = mdot_12*H_L
    else:    
        flux1[0] = mdot_12
        flux1[1] = mdot_12*u_R[0]
        flux1[2] = mdot_12*u_R[1]
        flux1[3] = mdot_12*u_R[2]
        flux1[4] = mdot_12*H_R
    for i in range(len(flux1)):
        flux_total[i] = flux1[i] + flux2[i]
    return flux_total

def invokeAUSMPlusUp(Q_L, Q_R, n_surface, inletMach):
    rho_L, u_x_L, u_y_L, u_z_L, E_L = defineInputValues(Q_L)
    rho_R, u_x_R, u_y_R, u_z_R, E_R = defineInputValues(Q_R)

    p_L = pp.calculatePressure(Q_L)
    p_R = pp.calculatePressure(Q_R)
    
    u_L = [u_x_L, u_y_L, u_z_L]
    u_R = [u_x_R, u_y_R, u_z_R]

    a_L = calculateSpeedOfSound(rho_L, p_L)
    a_R = calculateSpeedOfSound(rho_R, p_R)

    H_L = calculateEnthalpy(E_L, rho_L, p_L)
    H_R = calculateEnthalpy(E_R, rho_R, p_R)

    a_12 = calculateFaceSpeedOfSound(a_L, a_R)
    rho_12 = calculateFaceDensity(rho_L, rho_R)

    M_L, M_R = calculateLeftRightMach(n_surface, u_L, u_R, a_12)
    meanMach = calculateMeanMach(u_L, u_R, a_12)
    M_0 = calculateM0(meanMach, inletMach)
    f_a = calculateScalingFactor(M_0)

    M_12 = calculateFaceMach(M_L, M_R, M_0, meanMach, f_a, p_R, p_L, rho_12, a_12)
    mdot_12 = calculateMassFlux(M_12, a_12, rho_L, rho_R)
    p_12 = calculateFacePressure(M_L, M_R, rho_L, rho_R, p_L, p_R, f_a, a_12)
    fluxFinal = decideFinalFlux(M_12, mdot_12, p_12, u_L, u_R, H_L, H_R, n_surface)
    return fluxFinal