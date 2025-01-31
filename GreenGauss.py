import meshImporter as mesher
import numpy as np
import postProcesser as pp
import AUSMPlus as ausm

gamma = 1.4

def gradValue(faceCenterValue, cellIndex):
    cellVolume = mesher.calculateElementVolumePoly(cellIndex)
    testElementFaces = []
    gradQ_x = [0,0,0,0,0]
    gradQ_y = [0,0,0,0,0]
    gradQ_z = [0,0,0,0,0]
    totalGrad_x = [0,0,0,0,0]
    totalGrad_y = [0,0,0,0,0]
    totalGrad_z = [0,0,0,0,0]
    testElementFaces = mesher.findFacesOfCell(cellIndex)
    for i in range(len(testElementFaces)):
        faceIndex = testElementFaces[i][0]
        area = mesher.calculateFaceArea(faceIndex)
        faceNormalVec = mesher.calculateFaceNormalVector(faceIndex, cellIndex)
        n_x = faceNormalVec[0]
        n_y = faceNormalVec[1]
        n_z = faceNormalVec[2]
        for i in range(len(faceCenterValue)):
            gradQ_x[i] = (1/cellVolume)*faceCenterValue[i]*n_x*area
            gradQ_y[i] = (1/cellVolume)*faceCenterValue[i]*n_y*area
            gradQ_z[i] = (1/cellVolume)*faceCenterValue[i]*n_z*area
        totalGrad_x[i] = totalGrad_x[i] + gradQ_x[i]
        totalGrad_y[i] = totalGrad_y[i] + gradQ_y[i]
        totalGrad_z[i] = totalGrad_z[i] + gradQ_z[i]
    totalGrad = [totalGrad_x, totalGrad_y, totalGrad_z]
    totalGrad = np.array(totalGrad).T.tolist()
    return totalGrad

def interpolateValue(leftCellValue, rightCellValue):
    interpolatedValue = [0,0,0,0,0]
    for i in range(len(leftCellValue)):
        interpolatedValue[i] = (leftCellValue[i] + rightCellValue[i])/2
    return interpolatedValue

def calculateLeftRightStateVector(leftCellCenterValue, rightCellCenterValue, leftCellIndex, rightCellIndex, faceIndex):
    Q_l = [0,0,0,0,0]
    Q_r = [0,0,0,0,0]
    Q_old_l = leftCellCenterValue
    Q_old_r = rightCellCenterValue
    Q_face = interpolateValue(leftCellCenterValue, rightCellCenterValue)
    gradQ_left = gradValue(Q_face, leftCellIndex)
    gradQ_right = gradValue(Q_face, rightCellIndex)
    r_left = mesher.calculateInterpolationVector(leftCellIndex, faceIndex)
    r_right = mesher.calculateInterpolationVector(rightCellIndex, faceIndex)
    prod_left = np.dot(gradQ_left, r_left)
    prod_right = np.dot(gradQ_right, r_right)
    for i in range(len(Q_l)):
        Q_l[i] = Q_old_l[i] + prod_left[i]
        Q_r[i] = Q_old_r[i] + prod_right[i]
    return Q_l, Q_r

def calculateInletStateVector(rightCellCenterValue, ghostCellCenterValue, n_boundary):
    Q_ghost = ghostCellCenterValue
    Q_R = rightCellCenterValue
    rho = Q_R[0]
    u = Q_R[1]
    v = Q_R[2]
    w = Q_R[3]
    n_x = n_boundary[0]
    n_y = n_boundary[1]
    n_z = n_boundary[2]
    p = pp.calculatePressure(Q_R)
    a = ausm.calculateSpeedOfSound(rho, p)
    rho_ghost = Q_ghost[0]
    u_ghost = Q_ghost[1]
    v_ghost = Q_ghost[2]
    w_ghost = Q_ghost[3]
    p_ghost = pp.calculatePressure(Q_ghost)
    p_boundary = 0.5*(p_ghost + p - rho*a*(n_x*(u_ghost - u) + n_y*(v_ghost - v) + n_z*(w_ghost - w)))
    rho_boundary = rho_ghost + (p_boundary - p_ghost)/(a*a)
    u_boundary = u_ghost - n_x*(p_ghost - p_boundary)/(rho*a)
    v_boundary = v_ghost - n_y*(p_ghost - p_boundary)/(rho*a)
    w_boundary = w_ghost - n_z*(p_ghost - p_boundary)/(rho*a)
    Q_boundary = [rho_boundary, rho_boundary*u_boundary, rho_boundary*v_boundary, rho_boundary*w_boundary, p_boundary]
    Q_boundary[4] = pp.calculateEnergyFromPressure(Q_boundary, p_boundary)
    return Q_boundary

def calculateOutletStateVector(leftCellCenterValue, ghostCellCenterValue, leftCellIndex, faceIndex):
    Q_l = [0,0,0,0,0]
    Q_old_l = leftCellCenterValue
    Q_face = interpolateValue(leftCellCenterValue, ghostCellCenterValue)
    gradQ_left = gradValue(Q_face, leftCellIndex)
    r_left = mesher.calculateInterpolationVector(leftCellIndex, faceIndex)
    prod_left = np.dot(gradQ_left, r_left)
    for i in range(len(Q_l)):
        Q_l[i] = Q_old_l[i] + prod_left[i]*0.5
    return Q_l

def calculateInviscidWallStateVector(leftCellCenterValue, ghostCellCenterValue, leftCellIndex, faceIndex):
    Q_l = [0,0,0,0,0]
    Q_old_l = leftCellCenterValue
    Q_face = interpolateValue(leftCellCenterValue, ghostCellCenterValue)
    gradQ_left = gradValue(Q_face, leftCellIndex)
    r_left = mesher.calculateInterpolationVector(leftCellIndex, faceIndex)
    prod_left = np.dot(gradQ_left, r_left)
    for i in range(len(Q_l)):
        Q_l[i] = Q_old_l[i] + prod_left[i]*0.5
    return Q_l

def defineOutletGhostStateVector(leftCellCenterValue, n_boundary):
    Q_L = leftCellCenterValue
    Q_ghost = [0,0,0,0,0]
    rho = Q_L[0]
    u = Q_L[1]
    v = Q_L[2]
    w = Q_L[3]
    n_x = n_boundary[0]
    n_y = n_boundary[1]
    n_z = n_boundary[2]
    p = pp.calculatePressure(Q_L)
    a = ausm.calculateSpeedOfSound(rho, p)
    p_boundary = 1/gamma
    rho_boundary = rho + (p_boundary - p)/(a*a)
    u_boundary = u + n_x*(p - p_boundary)/(rho*a)
    v_boundary = v + n_y*(p - p_boundary)/(rho*a)
    w_boundary = w + n_z*(p - p_boundary)/(rho*a)
    Q_boundary = [rho_boundary, rho_boundary*u_boundary, rho_boundary*v_boundary, rho_boundary*w_boundary, p_boundary]
    Q_boundary[4] = pp.calculateEnergyFromPressure(Q_boundary, p_boundary)
    for i in range(len(Q_L)):
        Q_ghost[i] = 2*Q_boundary[i] - Q_L[i]
    return Q_ghost

def defineInviscidWallGhostStateVector(leftCellCenterValue, n_boundary):
    Q_L = leftCellCenterValue
    rho = Q_L[0]
    u = Q_L[1]
    v = Q_L[2]
    w = Q_L[3]
    n_x = n_boundary[0]
    n_y = n_boundary[1]
    n_z = n_boundary[2]
    p = pp.calculatePressure(Q_L)
    p_ghost = p
    rho_ghost = rho
    contraVel = n_x*u + n_y*v + n_z*w
    u_ghost = u - 2*contraVel*n_x
    v_ghost = v - 2*contraVel*n_y
    w_ghost = w - 2*contraVel*n_z
    Q_ghost = [rho_ghost, rho_ghost*u_ghost, rho_ghost*v_ghost, rho_ghost*w_ghost, p_ghost]
    Q_ghost[4] = pp.calculateEnergyFromPressure(Q_ghost, p_ghost)
    return Q_ghost