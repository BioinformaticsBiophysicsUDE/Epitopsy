'''
Created on Jun 29,  2011

@author: chris
'''

from numpy import array, sin, pi, deg2rad, dot, sqrt

'''
Imports at the bottom ... circular imports ... what a mess!
'''



class Result:
    '''
    classdocs
    '''
    X = 0
    Y = 1
    Z = 2
    def __init__(self):
        '''
        Constructor
        '''
        self.id = 1
        self.theta = 0.0
        self.phi = 0.0
        self.psi = 0.0
        self.coordinateTranslationX = 0.0
        self.coordinateTranslationY = 0.0
        self.coordinateTranslationZ = 0.0
        self.gridTranslationX = 0
        self.gridTranslationY = 0
        self.gridTranslationZ = 0
        self.thetaOffset = 0.0
        self.phiOffset = 0.0
        self.psiOffset = 0.0
        self.score = 0.0
    def getId(self):
        return self.id
    def setId(self, id):
        self.id = id
    def getPhi(self): 
        return self.phi
    def setPhi(self, phi): 
        self.phi = phi
    def getPsi(self): 
        return self.psi
    def setPsi(self, psi): 
        self.psi = psi
    def getTheta(self): 
        return self.theta
    def setTheta(self, theta): 
        self.theta = theta
    def getRotation(self): 
        return  [self.theta, self.phi, self.psi]
    def setRotation(self, theta, phi, psi): 
        self.theta = theta
        self.phi = phi
        self.psi = psi
    def scaleTranslationByMeshSize(self, meshSizeX, meshSizeY = 0, meshSizeZ = 0):
        if isinstance(meshSizeX, list):
            meshSizeY = meshSizeX[self.Y]
            meshSizeZ = meshSizeX[self.Z]
            meshSizeX = meshSizeX[self.X]
        self.setCoordinateTranslationX(meshSizeX * self.getGridTranslationX())
        self.setCoordinateTranslationY(meshSizeY * self.getGridTranslationY())
        self.setCoordinateTranslationZ(meshSizeZ * self.getGridTranslationZ())
    def appendToFile(self, filename):
        ResultWriter().append(filename, self)
    def writeToFile(self, filename):
        ResultWriter().write(filename, self)
    def getCoordinateTranslationX(self): 
        return self.coordinateTranslationX
    def setCoordinateTranslationX(self, coordinateTranslationX): 
        self.coordinateTranslationX = coordinateTranslationX
    def getCoordinateTranslationY(self): 
        return self.coordinateTranslationY
    def setCoordinateTranslationY(self, coordinateTranslationY):
        self.coordinateTranslationY = coordinateTranslationY
    def getCoordinateTranslationZ(self):
        return self.coordinateTranslationZ
    def setCoordinateTranslationZ(self, coordinateTranslationZ):
        self.coordinateTranslationZ = coordinateTranslationZ
    def getGridTranslationX(self):
        return self.gridTranslationX
    def setGridTranslationX(self, gridTranslationX):
        self.gridTranslationX = gridTranslationX
    def getGridTranslationY(self):
        return self.gridTranslationY
    def setGridTranslationY(self, gridTranslationY):
        self.gridTranslationY = gridTranslationY
    def getGridTranslationZ(self):
        return self.gridTranslationZ
    def setGridTranslationZ(self, gridTranslationZ):
        self.gridTranslationZ = gridTranslationZ
    def getGridTranslation(self):
        return  [self.getGridTranslationX(), self.getGridTranslationY(), self.getGridTranslationZ()]
    def setGridTranslation(self, gt):
        self.setGridTranslationX(gt[self.X])
        self.setGridTranslationY(gt[self.Y])
        self.setGridTranslationZ(gt[self.Z])
    def getCoordinateTranslation(self):
        return  [self.getCoordinateTranslationX(), self.getCoordinateTranslationY(), self.getCoordinateTranslationZ()]
    def setCoordinateTranslation(self, ct):
        self.setCoordinateTranslationX(ct[self.X])
        self.setCoordinateTranslationY(ct[self.Y])
        self.setCoordinateTranslationZ(ct[self.Z])
    def getPhiOffset(self):
        return self.phiOffset
    def setPhiOffset(self, phiOffset):
        self.phiOffset = phiOffset
    def getPsiOffset(self):
        return self.psiOffset
    def setPsiOffset(self, psiOffset):
        self.psiOffset = psiOffset
    def getThetaOffset(self):
        return self.thetaOffset
    def setThetaOffset(self, thetaOffset):
        self.thetaOffset = thetaOffset
    def getScore(self):
        return self.score
    def setScore(self, score):
        self.score = score


class ResultSet:
    '''
    classdocs
    '''
    def __init__(self, rs = None):
        if rs is None:
            self.results = []
        elif isinstance(rs, ResultSet):
            self.results = rs.results.clone()
        elif isinstance(rs, Result):
            self.results = rs
    '''*
     * Method to append a Result object to ResultSet
     * @param result Result object to be appended to ResultSet
     '''
    def append(self, result): 
        self.results.append(result)
    def remove(self, n): 
        self.results.remove(n)
    def get(self, n): 
        if(n >= self.size()):
            raise NameError("Index not in ResultSet.")
        else:
            return self.results[n]
    def getLast(self): 
        if(self.size() == 0):
            raise NameError("ResultSet is empty.")
        else:
            return self.get(self.size() - 1)
    def getBestScore(self): 
        if(self.size() == 0):
            raise NameError("ERROR: ResultSet is empty.")
        best = 0
        for i in range(self.size()): 
            if(self.results[i].getScore() > self.results[best].getScore()): 
                best = i
        return self.get(best)
    def getNBestScores(self, n): 
        if(self.size() < n):
            raise NameError("ResultSet shorter than requested number of items.")
        self.sortByScore()
        r = ResultSet()
        for i in range(n): 
            r.append(self.get(i))
        return r
    def size(self): 
        return len(self.results)
    def clear(self): 
        self.results = []
    def merge(self, rs): 
        for i in range(rs.size()): 
            if(rs.get(i) is not None): 
                self.append(rs.get(i))
    def sortByScore(self):
        # Implemented the comparator as a function below 
        self.results.sort(Comparator)
    def writeToFile(self, filename):
        ResultWriter().write(filename, self)
    def appendToFile(self, filename):
        ResultWriter().append(filename, self)
    def normalize(self, n): 
        if(n == 0.0): 
            raise NameError("ERROR: Cannot divide by 0. Normalization of scores failed!")
        for r in self.results: 
            r.setScore(r.getScore() / n)
    def meanAngle(self): 
        ma = array([0.0, 0.0, 0.0])
        for r in self.results:
            newVec = array([sin(deg2rad(r.theta) / 2), sin(deg2rad(r.phi) / 2), sin(deg2rad(r.psi) / 2)])
            ma = ma + newVec
        ma = ma * (1.0 / len(self.results))
        return ma
    def clone(self): 
        return ResultSet(self.results[:])
    def cluster(self, rs, maxDist): 
        clusters = []
        for r in rs.results: 
            clusters.append(ResultSet())
            clusters[-1].append(r)
        while(len(clusters) > 1): 
            bestJ = 0
            bestI = 0
            bestDist = 9999
            for i in range(len(clusters)): 
                iCenter = array(clusters[i].meanAngle())
                for j in range(len(clusters)): 
                    if(i != j): 
                        dist = iCenter - array((clusters[j].meanAngle()))
                        dist = sqrt(dot(dist, dist))
                        if(dist < bestDist): 
                            bestDist = dist
                            bestI = i
                            bestJ = j
            if(bestDist > maxDist):
                break
            clusters[bestI].merge(clusters[bestJ])
            clusters.remove(bestJ)
            bestI = 0
            bestJ = 0
            bestDist = 9999
        return clusters
    def getResults(self): 
        return self.results
    def setResults(self, results): 
        self.results = results
        
def Comparator(res1, res2):
    diff = res2.getScore() - res1.getScore()
    if diff > 0:
        return 1
    if diff < 0:
        return -1
    else:
        return 0

from ResultWriter import ResultWriter

