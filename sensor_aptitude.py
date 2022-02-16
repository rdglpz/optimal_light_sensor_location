import numpy as np


class SensorFitness():
    """
    
    """
    
    
    def __init__(self,NLTI,EAM,sensitivity,local_variograms,coordinates):
        """
        
        """
        print("Selct cost functions: \n 'xor','max'")
        self.NLTI = NLTI
        self.EAM = EAM
        self.sensitivity = sensitivity
        self.local_variograms = local_variograms
        self.coordinates = coordinates
        
    def selectFitnessFunction(self,s):
        
        if s=="xor":
            self.f = self.xor
        elif s=="max":
            self.f = self.maximum
        elif s=="cover":
            self.f = self.evaluateCoverage

        
    def validate_coordinates(self,iy,ix):
        """
        
        """
        
        iy = np.where(self.coordinates[:,0]==iy)
        ix = np.where(self.coordinates[:,1]==ix)
        
        return np.intersect1d(ix,iy)
        
        
    def coverage(self,X):
        """
        
        """
        n_sensors = int(len(X)/2)
        sensor_list = X.reshape(n_sensors,len(self.NLTI.shape))

        coverage = np.zeros((n_sensors,self.NLTI.shape[0],self.NLTI.shape[1]))
        
        for i,s in enumerate(sensor_list):

            sy,sx = s[0],s[1]
            ix = self.validate_coordinates(sy,sx)
            coverage[i] = np.zeros(self.NLTI.shape)
 
            if len(ix)>0: 
                map0to1 = 1/(1+self.local_variograms[ix]) 
                coverage[i] = map0to1*(map0to1<1)*self.sensitivity
        return coverage
    
    def coverage2(self,X):
        """
        
        """
        n_sensors = int(len(X)/2)
        sensor_list = X.reshape(n_sensors,len(self.NLTI.shape))

        coverage = np.zeros((n_sensors,self.NLTI.shape[0],self.NLTI.shape[1]))
        
        for i,s in enumerate(sensor_list):

            sy,sx = s[0],s[1]
            ix = self.validate_coordinates(sy,sx)
            coverage[i] = np.zeros(self.NLTI.shape)
            

 
            if len(ix)>0:
        
                coordinates = self.coordinates[ix][0]
                pi = self.NLTI[coordinates[0]][coordinates[1]]
                tvar = self.local_variograms[ix][0]
                
                outofrange = (tvar==0)*(pi**2/2)
                tvar[coordinates[0]][coordinates[1]]=0
                
                M = tvar+outofrange
                
                       
                
                lb = pi**2/2
                map0to1 = (-M+lb)/lb

                coverage[i] = map0to1*self.sensitivity
        return coverage
    
    def evaluateCoverage(self,X):
        """
        
        """
        n_sensors = int(len(X)/2)
        sensor_list = X.reshape(n_sensors,len(self.NLTI.shape))

        coverage = np.zeros((n_sensors,self.NLTI.shape[0],self.NLTI.shape[1]))
        
        for i,s in enumerate(sensor_list):

            sy,sx = s[0],s[1]
            ix = self.validate_coordinates(sy,sx)
            coverage[i] = np.zeros(self.NLTI.shape)
            
            if len(ix)>0:
                coverage[i] = self.local_variograms[ix][0]>0
                
        M = coverage
        
        
        return -np.sum(np.sum(M,axis=0)>0)
                
    def maximum(self,X):
        """
        
        
        """
    
        M = self.coverage2(X)
        #creamos n mapas de cobertutura de cada sensor        
        return -np.sum(np.max(M,axis=0))
    
    def xor(self,X):
        """
        
        
        """
    
        M = self.coverage(X)
        
        #XOR
        
        #generamos una mascara indicando los valores que nos interesa tomar en cuenta, los cuales son regiones donde no hay intersección de cobertura
        
        mask = np.sum(M>0,axis=0)==1

        return -(np.sum([mask]*len(M)*M))
    
    def showPositions(self,X):
        """
        show positions, coverage, histogram, covered sensitivity
        """
        
        n_sensors = int(len(X)/2)
        sensor_list = X.reshape(n_sensors,len(self.NLTI.shape))
        
        positions = np.zeros(np.shape(self.NLTI))
        for i,p in enumerate(sensor_list.astype(int)):
            positions[p[0]][p[1]] = i+1
        return positions
        
        
                
        
        
        
        

        
    
 