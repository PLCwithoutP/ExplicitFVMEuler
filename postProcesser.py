import math
import numpy as np
import os
import meshImporter as mesher

gamma = 1.4

def calculatePressure(Q):
    rho = Q[0]
    vel = calculateVelocity(Q)
    E = Q[4]
    e_int = E/rho - 0.5*(vel[0]**2 + vel[1]**2 + vel[2]**2)
    p = e_int * rho * (gamma - 1)
    return p

def calculateEnergyFromPressure(Q,p):
    rho = Q[0]
    vel = calculateVelocity(Q)
    magVel = vel[0]**2 + vel[1]**2 + vel[2]**2
    E = p/(gamma - 1) + (0.5*rho*magVel)    
    return E

def calculateVelocity(Q):
    u = Q[1]/Q[0]
    v = Q[2]/Q[0]
    w = Q[3]/Q[0]
    vel = [u, v, w]
    return vel

def calculateMach(Q):
    p = calculatePressure(Q)
    rho = Q[0]
    vel = np.linalg.norm(calculateVelocity(Q))
    a = math.sqrt((p*gamma)/rho)
    Ma = vel/a
    return Ma

def calculateSpeedOfSound(Q):
    p = calculatePressure(Q)
    rho = Q[0]
    a = math.sqrt((p*gamma)/rho)
    return a

def createScalarFieldData(time, fieldName, dimensions, fieldDataList):
    os.system(f"mkdir {time}")
    fieldFile = open(f"{time}/{fieldName}", "w")
    headerString = r""" /*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2406                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/"""
    foamHeaderString = r"""
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       volScalarField;"""
    locationString = f"\n    location    \"{time}\";\n"
    objectString = f"    object      {fieldName};\n"    
    dimensionsString = f"dimensions      {dimensions};\n\n"
    internalFieldHeader = "internalField      nonuniform List<scalar>\n"
    totalCellsString = f"{mesher.totalCells}\n"
    boundaryFieldHeader = "boundaryField\n" + "{\n"
    
    fieldFile.write(headerString)
    fieldFile.write(foamHeaderString)
    fieldFile.write(locationString)
    fieldFile.write(objectString)
    fieldFile.write("""}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
""")
    fieldFile.write(dimensionsString)
    fieldFile.write(internalFieldHeader)
    fieldFile.write(totalCellsString)
    fieldFile.write("(\n")
    for i in range(len(fieldDataList)):
        fieldFile.write(f"{fieldDataList[i]}" + "\n")
    fieldFile.write(")\n" + ";\n\n")
    fieldFile.write(boundaryFieldHeader)
    for i in range(len(mesher.parsedBoundary)):
        fieldFile.write(f"      {mesher.parsedBoundary[i][0]}\n")
        fieldFile.write("      {\n")
        fieldFile.write("            value      nonuniform List<scalar>\n")
        fieldFile.write(f"{mesher.parsedBoundary[i][2]}\n")
        fieldFile.write("(\n")
        for j in range(int(mesher.parsedBoundary[i][2])):
            fieldFile.write("1\n")
        fieldFile.write(")\n" + ";\n" + "      }\n")
    fieldFile.write("}\n")
    fieldFile.close()
    return

def createPressureFieldData(time, fieldDataList):
    createScalarFieldData(time, "p", "[1 -1 2 0 0 0 0]", fieldDataList)
    return

def createMachData(time, fieldDataList):
    createScalarFieldData(time, "Ma", "[0 0 0 0 0 0 0]", fieldDataList)
    return

def createUData(time, fieldDataList):
    createScalarFieldData(time, "u", "[0 1 -1 0 0 0 0]", fieldDataList)
    return

def createVData(time, fieldDataList):
    createScalarFieldData(time, "v", "[0 1 -1 0 0 0 0]", fieldDataList)
    return

def createAData(time, fieldDataList):
    createScalarFieldData(time, "a", "[0 1 -1 0 0 0 0]", fieldDataList)
    return

def logging(data):
    f = open("logFile.txt", "a")
    f.write(data + "\n")
    f.close()
    return
    