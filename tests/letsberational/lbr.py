from typing import List
import csv
import matplotlib.pyplot as plt 
import numpy as np

class NumericalResultLetsBeRational: 

    def __init__(self, name:str):
        self.name = name 
        self.data = self.loadData() 
    
    def loadData(self) -> dict[str, List[float]]:
        logMoneyness, beta, normalizedSigma, estimatedSigma, relativeDiffNormalizedSigma, = [],[], [], [], []
        with open('build/LetsBeRationalResult/'+self.name+'.csv', newline="") as csvfile:
            reader = csv.DictReader(csvfile)  
            for row in reader:
                logMoneyness.append(float(row['logMoneyness']))
                beta.append(float(row['beta']))
                normalizedSigma.append(float(row['normalizedSigma']))
                estimatedSigma.append(float(row['estimatedSigma']))
                relativeDiffNormalizedSigma.append(float(row['relativeDiffNormalizedSigma']))
        return {
            'logMoneyness':logMoneyness, 
            'beta' : beta,
            'normalizedSigma' : normalizedSigma,
            'estimatedSigma' : estimatedSigma,
            'relativeDiffNormalizedSigma' : relativeDiffNormalizedSigma
        }
    
    

    def plotErrorSurface(self):
        # Extract data
        x = np.array(self.data['logMoneyness'])
        y = np.array(self.data['normalizedSigma'])
        z = np.array(self.data['relativeDiffNormalizedSigma'])

        # Create grid
        x_unique = np.unique(x)
        y_unique = np.unique(y)
        X, Y = np.meshgrid(x_unique, y_unique)

        # Initialize Z grid
        Z = np.full(X.shape, np.nan)

        # Fill grid
        for xi, yi, zi in zip(x, y, z):
            ix = np.where(x_unique == xi)[0][0]
            iy = np.where(y_unique == yi)[0][0]
            Z[iy, ix] = zi

        # Plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, Z)

        ax.set_xlabel('logMoneyness')
        ax.set_ylabel('normalizedSigma')
        ax.set_zlabel('relativeDiffNormalizedSigma')

        plt.show()

    
NumericalResultLetsBeRational("Figure7").plotErrorSurface()
NumericalResultLetsBeRational("Figure8").plotErrorSurface()
NumericalResultLetsBeRational("Figure9").plotErrorSurface()
