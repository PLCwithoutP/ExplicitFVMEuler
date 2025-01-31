# ====================================================================
# Author    : Mehmet Furkan Tunc
# ID        : 511261165
# Title     : Unstructured Cell-Centered Hexagonal FVM 
# ====================================================================

import AUSMPlus as ausm
import GreenGauss as gG
import meshImporter as mesher
import postProcesser as pp
import numpy as np

Q_global = []
Q_global_old_iter = [] 
p_init_array = []
M_init = []
a_iter = []
p_iter = []
M_iter = []
u_iter = []
v_iter = []
magVel_iter = []

p_inlet = 4.5/1.4
rho_inlet = 2.666
a_inlet = ausm.calculateSpeedOfSound(rho_inlet, p_inlet)
inletMach = 2
initialMach = 0.57
u_inlet = a_inlet*initialMach
E_inlet = pp.calculateEnergyFromPressure([rho_inlet, rho_inlet*u_inlet, 0, 0, 2.5], p_inlet)
Q_inlet_L = [rho_inlet, rho_inlet*u_inlet, 0, 0, E_inlet]

p_init_value = 1/1.4
rho_init = 1
a_init = ausm.calculateSpeedOfSound(rho_init, p_init_value)
u_init = 0
E_init = pp.calculateEnergyFromPressure([rho_init, 0, 0, 0, 1], p_init_value)
Q_init = [rho_init,0,0,0,E_init]

tCurrent = 0
deltaT = 0.0001
tEnd = 0.3
iteration = int(tEnd / deltaT)

for i in range(mesher.totalCells):
    Q_global_old_iter.append(Q_init)
    Q_global.append(Q_init)

for eachQ in Q_global_old_iter:
    p_init_array.append(pp.calculatePressure(eachQ))
    M_init.append(pp.calculateMach(eachQ))
    
pp.createPressureFieldData(tCurrent, p_init_array)
pp.createMachData(tCurrent, M_init)
for i in range(iteration):
    print("##################################################################################################")
    print("Start of the time step")
    print("Current time is:", tCurrent)
    pp.logging("##################################################################################################")
    pp.logging("Start of the time step")
    pp.logging(f"Current time is: {tCurrent}")
    for cellID in range(mesher.totalCells):
    #for cellID in range(2):
        pp.logging(f"Cell number: {cellID}")
        cell = mesher.findFacesOfCell(cellID)
        cellVolume = mesher.calculateElementVolumePoly(cellID)
        product = -deltaT / cellVolume
        Q_cell = Q_global_old_iter[cellID]
        Q_newCell = [0,0,0,0,0]
        deltaQ = [0,0,0,0,0]
        for localFaceID in range(0,6):
            face = cell[localFaceID]
            globalFaceID = face[0]
            faceArea = mesher.calculateFaceArea(globalFaceID)
            faceNormal = mesher.calculateFaceNormalVector(globalFaceID, cellID)
            Q_faceCont = [0,0,0,0,0]
            if (globalFaceID >= mesher.minInletFaceID and globalFaceID <= mesher.maxInletFaceID):
                rightCellID = cellID
                Q_inlet_R = Q_global_old_iter[rightCellID]
                Q_L = gG.calculateInletStateVector(Q_inlet_R, Q_inlet_L, faceNormal)
                flux = ausm.invokeAUSMPlusUp(Q_L, Q_inlet_R, faceNormal, inletMach)
                for i in range(5):
                    Q_faceCont[i] = flux[i]*product*faceArea
                pp.logging(f"Inlet face {globalFaceID} deltaQ: " + str(Q_faceCont))
            elif (globalFaceID >= mesher.minWallFaceID and globalFaceID <= mesher.maxWallFaceID):
                leftCellID = cellID
                Q_wall_R = gG.defineInviscidWallGhostStateVector(Q_global_old_iter[leftCellID], faceNormal)
                Q_L = gG.calculateInviscidWallStateVector(Q_global_old_iter[leftCellID], Q_wall_R, leftCellID, globalFaceID)
                #Q_L = Q_global_old_iter[leftCellID]
                flux = ausm.invokeAUSMPlusUp(Q_L, Q_wall_R, faceNormal, inletMach)
                for i in range(5):
                    Q_faceCont[i] = flux[i]*product*faceArea
                pp.logging(f"Wall face {globalFaceID} deltaQ: " + str(Q_faceCont))
            elif (globalFaceID >= mesher.minOutletFaceID and globalFaceID <= mesher.maxOutletFaceID):
                leftCellID = cellID
                Q_outlet_R = gG.defineOutletGhostStateVector(Q_global_old_iter[leftCellID], faceNormal)
                #Q_outlet_R = Q_init
                Q_L = gG.calculateOutletStateVector(Q_global_old_iter[leftCellID], Q_outlet_R, leftCellID, globalFaceID)
                #Q_L = Q_global_old_iter[leftCellID]
                flux = ausm.invokeAUSMPlusUp(Q_L, Q_outlet_R, faceNormal, inletMach)
                for i in range(5):
                    Q_faceCont[i] = flux[i]*product*faceArea
                pp.logging(f"Outlet face {globalFaceID} deltaQ: " + str(Q_faceCont))
            elif (globalFaceID >= mesher.minEmptyFaceID and globalFaceID <= mesher.maxEmptyFaceID):
                continue
            else:
                leftCellID = cellID
                for i in range(len(faceNormal)):
                    if faceNormal[i] == -1:
                        rightCellID = mesher.findOwnerCells(globalFaceID)
                        break
                    if faceNormal[i] == 1:
                        rightCellID = mesher.findNbourCells(globalFaceID)
                        break
                Q_L, Q_R = gG.calculateLeftRightStateVector(Q_global_old_iter[leftCellID], Q_global_old_iter[rightCellID], leftCellID, rightCellID, globalFaceID)
                #print("faceID:\t", globalFaceID)
                #print("leftCellID:\t", leftCellID)
                #print("rightCellID:\t", rightCellID)
                #print("faceNormal:\t", faceNormal)
                #Q_L, Q_R = Q_global_old_iter[leftCellID], Q_global_old_iter[rightCellID]
                #print("Q_L ", Q_L)
                #print("Q_R ", Q_R)
                flux = ausm.invokeAUSMPlusUp(Q_L, Q_R, faceNormal, inletMach)
                #print("flux:\t", flux)
                for i in range(5):
                    Q_faceCont[i] = flux[i]*product*faceArea
                pp.logging(f"Internal face {globalFaceID} deltaQ: " + str(Q_faceCont))
            for i in range(5):
                deltaQ[i] = deltaQ[i] + Q_faceCont[i]
        pp.logging(f"Total cell {cellID} flux is: {deltaQ}")
        for i in range(5):
            Q_newCell[i] = deltaQ[i] + Q_cell[i]
        Q_global[cellID] = Q_newCell
    tCurrent = tCurrent + deltaT
    Q_global_old_iter = Q_global
    print("End of the time step")
    print("##################################################################################################")
    pp.logging("End of the time step")
    pp.logging("##################################################################################################")
    for eachQ in Q_global:
        p_iter.append(pp.calculatePressure(eachQ))
        M_iter.append(pp.calculateMach(eachQ))
        u_iter.append(pp.calculateVelocity(eachQ)[0])
        v_iter.append(pp.calculateVelocity(eachQ)[1])
        magVel_iter.append(pp.calculateVelocity(eachQ)[0]**2 + pp.calculateVelocity(eachQ)[1]**2)
        a_iter.append(pp.calculateSpeedOfSound(eachQ))

    pp.createPressureFieldData(tCurrent, p_iter)
    pp.createMachData(tCurrent, M_iter)
    pp.createUData(tCurrent, u_iter)
    pp.createVData(tCurrent, v_iter)
    pp.createScalarFieldData(tCurrent, "magVelocity", "[0 1 -1 0 0 0]", magVel_iter)
    pp.createAData(tCurrent, a_iter)
    p_iter = []
    M_iter = []
    u_iter = []
    v_iter = []
    magVel_iter = []
    a_iter = []
