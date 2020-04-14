# Project 1 - Transferência de Calor e Mecânica dos Sólidos - Caio Fauza, Pedro Paulo Telho, Luiz Vitor Germanos
import math
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from funcoesTermosol import plota, importa, geraSaida
[nn, N, nm, Inc, nc, F, nr, R] = importa('contest.xlsx')
info = dict()
matrix = []
degrees = []
lenghts = []
barInfo = []
converged = False
equationSolver = input("Choose the iterative method to solve the equations\n(Type 'Gauss' or 'Jacobi', without the quotes): ")

#Filling in info according to the incidence
for i in range(1, nm + 1):
    info[i] = [Inc[i-1][0], Inc[i-1][1]]

for i in range(1, len(info) + 1):
    #Calculation of bar sizes, sin and cos
    x1 = N[0][int(info[i][0] - 1)]
    x2 = N[0][int(info[i][1] - 1)]
    y1 = N[1][int(info[i][0] - 1)]
    y2 = N[1][int(info[i][1] - 1)]
    length = math.sqrt((x2-x1)**2 + (y2-y1)**2)
    lenghts.append(length)
    s = (y2-y1)/length
    c = (x2-x1)/length
    barInfo.append([s, c])

    #Creation of the stiffness matrix * (E * A)/L
    matrix.append(np.dot(np.matrix([[c**2, c*s, -c**2, -c*s],
                                    [c*s, s**2, -c*s, -s**2],
                                    [-c**2, -c*s, c**2, c*s],
                                    [-c*s, -s**2, c*s, s**2]]), (Inc[i-1][2]*Inc[i-1][3])/length))

    #Organization of degrees of freedom
    var1 = (int(info[i][0] - 1))*2
    var2 = (int(info[i][1] - 1))*2
    degrees.append([[var1, var1+1], [var2, var2+1]])



#Global array initialization with zeros
globalMatrix = np.zeros((2*int(nn), 2*int(nn)))


#Superposition
for i in range(len(degrees)):
    globalMatrix[degrees[i][0][0]: degrees[i][0][1] + 1, degrees[i]
                 [0][0]: degrees[i][0][1] + 1] += matrix[i][0:2, 0:2]
    globalMatrix[degrees[i][0][0]: degrees[i][0][1] + 1, degrees[i]
                 [1][0]: degrees[i][1][1] + 1] += matrix[i][0:2, 2:4]
    globalMatrix[degrees[i][1][0]: degrees[i][1][1] + 1, degrees[i]
                 [0][0]: degrees[i][0][1] + 1] += matrix[i][2:4, 0:2]
    globalMatrix[degrees[i][1][0]: degrees[i][1][1] + 1, degrees[i]
                 [1][0]: degrees[i][1][1] + 1] += matrix[i][2:4, 2:4]

#Application of boundary conditions
globalMatrixRaw = globalMatrix

for i in range(nr - 1, -1, -1):
    globalMatrix = np.delete(globalMatrix, int(R[i]), 0)
    F = np.delete(F, int(R[i]), 0)
    try:
        globalMatrix = np.delete(globalMatrix, int(R[i]), 1)
    except:
        pass

#Equation solver - Gauss  
def gaussSolver(a, b, it):
    u = np.zeros((len(F)))
    precision = 10**(-14)
    for h in range(it):
        uPrevious = np.copy(u)
        for i in range(len(a)):
            for j in range(len(a)):
                if(j != i):
                    u[i] += a[i][j] * u[j]
            u[i] = (b[i] - u[i]) / a[i][i]
            
        if(h > 0):
            convergency = np.amax(np.abs(np.divide(u - uPrevious, u)))
            if(convergency < precision):
                global converged
                converged = True
                print('Converged. Ready to get results and check if structure collapses...')
                return u, h + 1

    return u, it
        
#Equation solver - Jacobi
def jacobiSolver(a, b, it):
    u = np.zeros((len(F)))
    precision = 10**(-14)
    for h in range(it):
        uPrevious = np.copy(u)
        for i in range(len(a)):
            for j in range(len(a)):
                if(j != i):
                    u[i] += a[i][j] * uPrevious[j]
            u[i] = (b[i] - u[i]) / a[i][i]
            
        if(h > 0):
            convergency = np.amax(np.abs(np.divide(u - uPrevious, u)))
            if(convergency < precision):
                global converged
                converged = True
                print('Converged. Ready to get results and check if structure collapses...')
                return u, h + 1
       
    return u, it
    

if(equationSolver == "Gauss"):
    solutionG = gaussSolver(globalMatrix, F, 5500)
    u = solutionG[0]
    print('Gauss - ' + str(solutionG[1]) + ' iterations were needed.') 
else:
    solutionJ = jacobiSolver(globalMatrix, F, 5500)
    u = solutionJ[0]
    print('Jacobi - ' + str(solutionJ[1]) + ' iterations were needed.')


#Get u in expanded form
uExpanded = np.zeros((2*(nn)))
uIndex = []
for i in range((nn)*2):
    if(i not in R):
        uIndex.append(i)

for i in range(len(u)):
    uExpanded[uIndex[i]] = u[i]     


#Reactions
reactionMatrix = globalMatrixRaw.dot(uExpanded)
reactions = []
for i in range(len(reactionMatrix)):
    if(i in R):
        reactions.append(reactionMatrix[i])


#Deformations
deformations = []
for i in range(1, len(info) + 1):
    s = barInfo[i-1][0]
    c = barInfo[i-1][1]
    u1 = uExpanded[degrees[i-1][0][0]]
    v1 = uExpanded[degrees[i-1][0][1]]
    u2 = uExpanded[degrees[i-1][1][0]]
    v2 = uExpanded[degrees[i-1][1][1]]
    uV = np.array([[u1], [v1], [u2], [v2]])
    factor = np.array([-c, -s, c, s])
    deformations.append(float((np.dot(factor, uV))*(1/lenghts[i-1])))


#Tensions
tensions = []
for i in range(len(info)):
    tensions.append(deformations[i]*Inc[i][2])

#Normal internal efforts
normals = []
for i in range(len(info)):
    normals.append(tensions[i]*Inc[i][3])

geraSaida('result', reactions, uExpanded, deformations, normals, tensions)

#Check if structure collapses
def checkCollapse():
    if(converged):
        for i in range(len(tensions)):
            if(abs(tensions[i]) > 18E6):
                print("Stress greater than 18 MPa found: " +  str(tensions[i]))

        for i in range(len(deformations)):
            if(abs(deformations[i] > 0.05)):
                print("Deformation greater than 5% found: " + str(deformations[i]))

        for i in range(len(uExpanded)):
            if(abs(uExpanded[i]) > 0.02):
                print("Displacement greater than 20mm found: " + str(uExpanded[i]))

        for i in lenghts:
            if i > 0.11:
                print("Bar size greater than 110mm found: " + str(i))
    else:
        print("Cannot check if structure collapses. Iterative method not had converged. Please check structure size and choose the adequated method.")

checkCollapse()
plota(N, Inc)


