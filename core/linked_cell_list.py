from .sphere import Sphere

class LinkedCellList:

    Lbox = [10., 10., 10.]
    totAtoms = int(1)
    maxCellsPerAxis = [int(200), int(200), int(200)]       #Max cells per axis to limit RAM usage
    minCellWidth = 1.                        #Min cell width (at least max sigma of the spheres in the system)
    cellsPerAxis = [int(10), int(10), int(10)]
    cellWidth = [1., 1., 1.]
    totCells = int(1000)
    neighboringIndexes = [int(2), int(2), int(2)]   #If the volume is too small to have 3 cells per axis

    headObj = [-1]         #head of the cell (contains first sphere)
    listObj = [-1]         #list of all the other spheres

    def __init__(self, Lbox=[10., 10., 10.], minCellWidth=1., totAtoms=1):

        self.Lbox = Lbox
        self.minCellWidth = minCellWidth
        self.totAtoms = totAtoms
        self.cellsPerAxis = [int(l / minCellWidth) for l in Lbox]
        for ax in range(3):
            self.cellsPerAxis[ax] = int(self.Lbox[ax] / self.minCellWidth)
            if self.cellsPerAxis[ax] > self.maxCellsPerAxis[ax]: 
                self.cellsPerAxis[ax] = self.maxCellsPerAxis[ax]
            self.cellWidth[ax] = self.Lbox[ax] / float(self.cellsPerAxis[ax])
        
        self.totCells = int(self.cellsPerAxis[0] * self.cellsPerAxis[1] * self.cellsPerAxis[2])

        self.headObj = [-1]*self.totCells
        self.listObj = [-1]*self.totAtoms

        for ax in range(3):
            if self.cellsPerAxis[ax] < 3: 
                self.neighboringIndexes[ax] = self.cellsPerAxis[ax] - 1

    def calculateObjectCell(self, cm=[0., 0., 0.]) -> list[int]:
        cellIndex = [0, 0, 0, 0]
        for ax in range(3):
            cellIndex[ax] = int((cm[ax] + 0.5 * self.Lbox[ax]) / self.cellWidth[ax])
        cellIndex[3] = cellIndex[0] + cellIndex[1] * self.cellsPerAxis[0] + cellIndex[2] * self.cellsPerAxis[0] * self.cellsPerAxis[1]
    

        return cellIndex

    def addObjectToList(self, objID: int, cm=[0., 0., 0.]):
        cellIndex = self.calculateObjectCell(cm)

        self.listObj[objID] = self.headObj[cellIndex[3]]
        self.headObj[cellIndex[3]] = objID


    def overlapCheck(self, sph : Sphere, spheres : list[Sphere]):
        objCellIndex = self.calculateObjectCell(sph.cm)

        cellIndex = [0, 0, 0, 0]
        for kx in range(-1, self.neighboringIndexes[0]):
            for ky in range(-1, self.neighboringIndexes[1]):
                for kz in range(-1, self.neighboringIndexes[2]):
        
                    cellIndex[0] = objCellIndex[0] + kx
                    cellIndex[1] = objCellIndex[1] + ky
                    cellIndex[2] = objCellIndex[2] + kz

                    for ax in range(3):
                        if cellIndex[ax] >= self.cellsPerAxis[ax]: 
                            cellIndex[ax] -= self.cellsPerAxis[ax]
                        elif cellIndex[ax] < 0:
                            cellIndex[ax] += self.cellsPerAxis[ax]
                    
                    cellIndex[3] = cellIndex[0] + cellIndex[1] * self.cellsPerAxis[0] + cellIndex[2] * self.cellsPerAxis[0] * self.cellsPerAxis[1]
                    cellIndex[3] = self.headObj[int(cellIndex[3])]

                    while cellIndex[3] >= 0 :
                        if sph.overlap(spheres[cellIndex[3]], self.Lbox):
                            return True
                        
                        cellIndex[3] = self.listObj[cellIndex[3]]
                    
        return False
