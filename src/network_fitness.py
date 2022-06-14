import numpy as np


class NetworkFitness():
    """
    
    """
    
    
    def __init__(self,NLTI,EAM,sensitivity,local_variograms,local_variograms_m,coordinates):
        """
        
        """
        print("Selct cost functions: \n 'xor','max' or 'cover'")
        self.NLTI = NLTI
        self.EAM = EAM
        self.sensitivity = sensitivity
        self.local_variograms = local_variograms
        self.local_variograms_m = local_variograms_m
        self.coordinates = coordinates
        
    def selectFitnessFunction(self,s):
        
        if s=="xor":
            self.f = self.xor
        elif s=="max":
            self.f = self.maximum
        elif s=="min":
            self.f = self.minimum
        elif s=="cover":
            self.f = self.bruteCoverage
        elif s=="explicability":
            self.f = self.explicability

        
    def validate_coordinates(self,iy,ix):
        """
        
        """
        
        iy = np.where(self.coordinates[:,0]==iy)
        ix = np.where(self.coordinates[:,1]==ix)
        
        return np.intersect1d(ix,iy)
        
        
    def coverage(self,X):
        """
        
        """
        
        #get the number of sensors
        X = np.array(X)
        n_sensors = int(len(X)/2)
        
        #get the coordinate sensor list
        dimensions = len(self.NLTI.shape)
        sensor_list = X.reshape(n_sensors,dimensions).astype(int)

        #empty coverage of each sensor
        coverage = np.zeros((n_sensors,self.NLTI.shape[0],self.NLTI.shape[1]))
        
        for i,s in enumerate(sensor_list):

            
            sy,sx = s[0],s[1]
            ix = self.validate_coordinates(sy,sx)
            
            
            coverage[i] = np.zeros(self.NLTI.shape)
            
            #descartamos aquella región con semivariograma >= 2
            #esto se puede mejorar estimando el variograma teórico y 
            #considerar  el 95% de las semivarianzas antes de tocar umbral
            # sill, a partir de la primera semivarianza C_0
            bounded_local_var = (self.local_variograms[ix])*(self.local_variograms[ix]<=50)
 
            if len(ix)>0: 
                map0to1 = 1/(1+bounded_local_var) 
                coverage[i] = map0to1*(map0to1<1)*self.sensitivity
        return coverage
    
    def coverage_explicability(self,X):
        """
        
        """
        
        #get the number of sensors
        n_sensors = int(len(X)/2)
        
        #get the coordinate sensor list
        sensor_list = X.reshape(n_sensors,len(self.NLTI.shape))

        #empty coverage of each sensor
        coverage = np.zeros((n_sensors,self.NLTI.shape[0],self.NLTI.shape[1]))
        
        for i,s in enumerate(sensor_list):

            
            sy,sx = s[0],s[1]
            ix = self.validate_coordinates(sy,sx)
            
            
            coverage[i] = np.zeros(self.NLTI.shape)
            
            #descartamos aquella región con semivariograma >= 2
            #esto se puede mejorar estimando el variograma teórico y 
            #considerar  el 95% de las semivarianzas antes de tocar umbral
            # sill, a partir de la primera semivarianza C_0
            v = self.local_variograms[ix]
            M = v > 0
 
            if len(ix)>0: 

                coverage[i] = M/(1+v)
            
        return coverage
    
    def coverage2(self,X):
        """
        
        """
        
        #get the numner of sensors.
        n_sensors = int(len(X)/2)
        
        #list of n elements with 2 
        sensor_list = X.reshape(n_sensors,len(self.NLTI.shape))

        #coverate_layers
        coverage = np.zeros((n_sensors,self.NLTI.shape[0],self.NLTI.shape[1]))
        
        for i,s in enumerate(sensor_list):

            sy,sx = s[0],s[1]
            ix = self.validate_coordinates(sy,sx)
            coverage[i] = np.zeros(self.NLTI.shape)
            
            if len(ix)>0:
                
                coordinates = self.coordinates[ix][0]
                pi = self.NLTI[coordinates[0]][coordinates[1]]
                tvar = self.local_variograms[ix][0]
                
                tvar_m = self.local_variograms_m[ix][0]
                
                #outofrange = (tvar==0)*(pi**2/2)
                outofrange = (tvar_m==0)*(np.max(self.NLTI)**2/2)
                inrange = (tvar_m==1)*tvar

                
                #M = tvar+outofrange
                M = inrange+outofrange
                
                lb = (np.max(self.NLTI))**2/2
                map0to1 = (lb-M)/(lb)

                coverage[i] = map0to1*self.sensitivity
        return coverage
    
    def coverMaps(self,X):
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
            
        
        
        return coverage
                
    def maximum(self,X):
        """
        
        
        """
    
        M = self.coverage2(X)
        #creamos n mapas de cobertutura de cada sensor   
        
        return -np.sum(np.max(M,axis=0))/np.sum(self.sensitivity)
    
    def minimum(self,X):
        """
        
        
        """
    
        M = self.coverage(X)
        #creamos n mapas de cobertutura de cada sensor   
        
        return -np.sum(np.min(M,axis=0))
    
    def xor(self,X):
        """
        
        
        """
    
        M = self.coverage(X)
        
        #XOR
        
        #generamos una mascara indicando los valores que nos interesa tomar en cuenta, los cuales son regiones donde no hay intersección de cobertura
        
        mask = np.sum(M>0,axis=0)==1

        return -(np.sum([mask]*len(M)*M))
    
    def bruteCoverage(self,X):
        """
        
        """
        
        M = self.coverMaps(X)
        
        
        return -np.sum(np.sum(M,axis=0)>0)
    
    def explicability(self,X):
        
        T = np.sum(self.sensitivity)
        E = self.coverage_explicability(X)
        
        C = np.max(E,axis=0)*self.sensitivity
        
        
        
        #creamos n mapas de cobertutura de cada sensor   
        
        return -100*np.sum(C)/T
        
        
        
        
    
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
    
    def showVariogram(self,X):
        """
        
        """
        
    def project(self,X):
        """
        
        
        """
        
        dim = self.NLTI.shape
        IMGP = np.copy(self.coverage2(X))
        R = np.zeros((len(IMGP)+1,dim[0],dim[1]))
        dummy = np.ones(dim)*-1
        R[0] = dummy

        for i in range(1,len(R)):
            outofrange = (IMGP[i-1]==0)*-1
            inrange = (IMGP[i-1]!=0)*IMGP[i-1]
            R[i] = outofrange+inrange

        return np.argmax(R,axis=0)*(self.NLTI >0)
        
        
        
                
        
        
        
        

        
    
 