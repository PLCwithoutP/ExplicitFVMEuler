import numpy as np

def meshFileImport(path):
    f = open(path, "r")
    data = f.read()
    content = data.split('\n')
    return content

points = meshFileImport("./constant/polyMesh/points")
faces = meshFileImport("./constant/polyMesh/faces")
owners = meshFileImport("./constant/polyMesh/owner")
nbours = meshFileImport("./constant/polyMesh/neighbour")
boundary = meshFileImport("./constant/polyMesh/boundary")

def checkDataType(content):
    for line in content:
        if ("points;" in line):
            dataType = 0
        if ("faces;" in line):
            dataType = 1
        if ("owner;" in line):
            dataType = 2
        if ("neighbour;" in line):
            dataType = 3
        if ("boundary;" in line):
            dataType = 4
    return dataType

def findTotalCellNumbers(content):
    for line in content:
        if ('note' in line):
            totalCellNumber = int(line.split('nCells:')[1].split(' ')[0])
    return totalCellNumber

def parseData(content, dataType):
    parsedData = []
    boundaryNames = []
    boundaryTypes = []
    boundaryNFaces = []
    boundaryStartFaces = []
    singleDataID = 0
    totalNumber = 0
    for line in content:
        singleData = []
        if ((dataType == 0 or dataType == 1)):
            if (len(line.split('(')) != 1):
                if (line.split('(')[1] != ''):
                    singleData.append(singleDataID)
                    coordinates = line.split('(')[1].split(')')[0].split(' ')
                    for i in range(len(coordinates)):
                        singleData.append(coordinates[i])
                    parsedData.append(singleData)
                    singleDataID = singleDataID + 1               
                else:
                    totalNumber = int(content[content.index(line) - 1]);
        elif ((dataType == 2 or dataType == 3)):
            if (line.isnumeric()):
                if ('(' in content[content.index(line) + 1]):
                    totalNumber = int(line)
                else:
                    singleData.append(singleDataID)
                    singleData.append(int(line))
                    parsedData.append(singleData)
                    singleDataID = singleDataID + 1
        else:
            if ('(' in line):
                if ('inGroups' not in line):
                    totalNumber = int(content[content.index(line) - 1])
            if ('type' in line):
                boundaryTypes.append(line.split(' ')[-1].split(';')[0])
            if ('nFaces' in line):
                boundaryNFaces.append(line.split(' ')[-1].split(';')[0])
            if ('startFace' in line):
                boundaryStartFaces.append(line.split(' ')[-1].split(';')[0])
            if (line.split(' ')[-1].isalpha() and line != 'FoamFile'):
                boundaryNames.append(line.split(' ')[-1])
    for i in range(len(boundaryTypes)):
        singleData.append(boundaryNames[i])
        singleData.append(boundaryTypes[i])
        singleData.append(boundaryNFaces[i])
        singleData.append(boundaryStartFaces[i])
        parsedData.append(singleData)
        singleData = []
    return parsedData, totalNumber

def readVertexCoordinates(index):
    xCoordinate = parsedPoints[index][1]
    yCoordinate = parsedPoints[index][2]
    zCoordinate = parsedPoints[index][3]
    coordinate = [xCoordinate, yCoordinate, zCoordinate]
    return coordinate

def calculateFaceArea(faceIndex):
    face = facesWCoords[faceIndex]
    point1 = face[0]
    point2 = face[1]
    point3 = face[2]
    point4 = face[3]
    vector1 = [0, 0, 0]
    vector2 = [0, 0, 0]
    for i in range(0,3):
        point1[i] = float(point1[i])
        point2[i] = float(point2[i])
        point3[i] = float(point3[i])
        point4[i] = float(point4[i])
        vector1[i] = point1[i] - point2[i]
        vector2[i] = point1[i] - point4[i]
    area = (np.cross(vector2, vector1))
    area = np.linalg.norm(area)
    return area

def calculateFaceNormalVector(faceIndex, cellIndex):
    face = facesWCoords[faceIndex]
    point1 = face[0]
    point2 = face[1]
    point3 = face[2]
    point4 = face[3]
    pointCenter = calculateFaceCenterPoint(faceIndex)
    vector1 = [0, 0, 0]
    vector2 = [0, 0, 0]
    vector3 = [0, 0, 0]
    vector4 = [0, 0, 0]
    for i in range(0,3):
        point1[i] = float(point1[i])
        point2[i] = float(point2[i])
        point3[i] = float(point3[i])
        point4[i] = float(point4[i])
        vector1[i] = point1[i] - pointCenter[i]
        vector2[i] = point2[i] - pointCenter[i]
        vector3[i] = point3[i] - pointCenter[i]
        vector4[i] = point4[i] - pointCenter[i]
    faceOrthoVector = 0.5*(np.cross(vector1, vector2) + np.cross(vector2, vector3) 
                + np.cross(vector3, vector4) + np.cross(vector4, vector1))
    faceNormalVector = faceOrthoVector/np.linalg.norm(faceOrthoVector)
    if not (cellIndex == parsedOwners[faceIndex][1]):
        for i in range(len(faceNormalVector)):
            faceNormalVector[i] = (-1)*faceNormalVector[i]
    return faceNormalVector
    
def calculateEdgeVector(faceIndex):
    face = facesWCoords[faceIndex]
    point1 = face[0]
    point2 = face[1]
    vector1 = [0, 0, 0]
    for i in range(0,3):
        point1[i] = float(point1[i])
        point2[i] = float(point2[i])
        vector1[i] = point2[i] - point1[i]
    return vector1

def calculateFaceCenterPoint(faceIndex):
    face = facesWCoords[faceIndex]
    point1 = face[0]
    point2 = face[1]
    point3 = face[2]
    point4 = face[3]
    for i in range(0,3):
        point1[i] = float(point1[i])
        point2[i] = float(point2[i])
        point3[i] = float(point3[i])
        point4[i] = float(point4[i])
    xMax = max(point1[0], point2[0], point3[0], point4[0])
    yMax = max(point1[1], point2[1], point3[1], point4[1])
    zMax = max(point1[2], point2[2], point3[2], point4[2])
    xMin = min(point1[0], point2[0], point3[0], point4[0])
    yMin = min(point1[1], point2[1], point3[1], point4[1])
    zMin = min(point1[2], point2[2], point3[2], point4[2])
    xCoord = (xMax + xMin)/2
    yCoord = (yMax + yMin)/2
    zCoord = (zMax + zMin)/2
    centerPoint = [xCoord, yCoord, zCoord]
    return centerPoint

def calculateElementCenterPoint(cellIndex):
    facesOfElement = []
    xCenter = 0
    yCenter = 0
    zCenter = 0
    facesOfElement = findFacesOfCell(cellIndex)
    for i in range(len(facesOfElement)):
        faceCenter = calculateFaceCenterPoint(int(facesOfElement[i][0]))
        xCenter = xCenter + faceCenter[0]
        yCenter = yCenter + faceCenter[1]
        zCenter = zCenter + faceCenter[2]
    xCenter = xCenter/6
    yCenter = yCenter/6
    zCenter = zCenter/6
    centerPoint = [xCenter, yCenter, zCenter]
    return centerPoint

# =============================================================================
# def findFacesOfCell(cellIndex):
#      facesOfElement = []
#      for i in range(len(parsedOwners)):
#          if (parsedOwners[i][1] == cellIndex):
#              facesOfElement.append(parsedFaces[i])
#      for i in range(len(parsedNbours)):
#          if (parsedNbours[i][1] == cellIndex):
#              facesOfElement.append(parsedFaces[i])
#      return facesOfElement
# =============================================================================

def findFacesOfCell(cellIndex):
    facesOfElement = []

    for i in range(len(parsedFaces)):  # Iterate over faces in the given order
        faceID = i  # Face index in OpenFOAM

        # Check if the face belongs to the given cell
        if parsedOwners[faceID][1] == cellIndex or (faceID < len(parsedNbours) and parsedNbours[faceID][1] == cellIndex):
            facesOfElement.append(parsedFaces[faceID])

    return facesOfElement  # Preserves order from faces file


# Vector'den kaynaklı volume hatalı
def calculateElementVolumePoly(cellIndex):
    facesOfElement = []
    volume = 0
    facesOfElement = findFacesOfCell(cellIndex)
    for face in facesOfElement:
        faceArea = calculateFaceArea(face[0])
        faceCenter = calculateFaceCenterPoint(face[0])
        faceNormal = calculateFaceNormalVector(face[0], cellIndex)
        volume = volume + (faceArea*(np.dot(faceCenter,faceNormal)))
    return abs(volume)/3
    

def calculateInterpolationVector(cellIndex, faceIndex):
    cellCenter = calculateElementCenterPoint(cellIndex)
    faceCenter = calculateFaceCenterPoint(faceIndex)
    vector = [0,0,0]
    for i in range(len(cellCenter)):
        vector[i] = cellCenter[i] - faceCenter[i]
    return vector

def determineLimitsOfBoundary(boundary):
    minFaceID = int(boundary[3])
    numberOfFaces = int(boundary[2])
    maxFaceID = minFaceID + numberOfFaces - 1
    return minFaceID, maxFaceID

def findOwnerCells(faceID):
    leftCellID = parsedOwners[faceID][1]
    return leftCellID

def findNbourCells(faceID):
    rightCellID = parsedNbours[faceID][1]
    return rightCellID

pointDataType = checkDataType(points)
faceDataType = checkDataType(faces)
ownerDataType = checkDataType(owners)
nboursDataType = checkDataType(nbours)
boundaryDataType = checkDataType(boundary)

parsedPoints, totalPoints = parseData(points, pointDataType)
parsedFaces, totalFaces = parseData(faces, faceDataType)
parsedOwners, totalOwners = parseData(owners, ownerDataType)
parsedNbours, totalNbours = parseData(nbours, nboursDataType)
parsedBoundary, totalBoundary = parseData(boundary, boundaryDataType)
totalCells = findTotalCellNumbers(owners)

facesWCoords = []
faceWCoords = []
for face in parsedFaces:
    for j in range(1,5):
        faceWCoords.append(readVertexCoordinates(int(face[j])))
    facesWCoords.append(faceWCoords)
    faceWCoords = []
    
for boundary in parsedBoundary:
    if (boundary[1] == 'wall'):
        minWallFaceID, maxWallFaceID = determineLimitsOfBoundary(boundary)
    elif (boundary[0] == 'inlet'):
        minInletFaceID, maxInletFaceID = determineLimitsOfBoundary(boundary)
    elif (boundary[0] == 'outlet'):
        minOutletFaceID, maxOutletFaceID = determineLimitsOfBoundary(boundary)
    else:
        minEmptyFaceID, maxEmptyFaceID = determineLimitsOfBoundary(boundary)