import itertools as it


class rGrowing():
    
    
    
    
    
    def __init__(self,Img):
        '''
        
        '''
        import numpy as np
        
        self.Img = Img
        self.rng = list(range(-1,2))
        self.A = np.zeros(Img.shape)
        self.coords = list([])
        
    
    def setParams(self,params):
        
        if params["type"]=="semivar":
            
            self.si = params["si"]
            self.atol = params["atol"]
            self.dtol = params["dtol"]
            

        
    def getRegion(self,sj):
        
        '''
        
        '''
        import numpy as np
        
        
        si = self.si
        self.sj=sj
  
        

        self.stack = []        
        visit = np.zeros(self.Img.shape)
        visit[sj] = 1
        self.coords = list([])
        self.coords.append(sj)
        self.stack.append(sj)
        
        while True if self.stack else False:
#        for i in range(1):

            N = list([])
            
            p = self.stack.pop()
            
            pr = it.product(self.rng,self.rng)

            for px in pr:

                y = p[0]+px[0]
                y = 0 if y < 0 else y
                y = self.Img.shape[0]-1 if y > self.Img.shape[0]-1 else y

                x = p[1]+px[1]
                x = 0 if x < 0 else x
                x = self.Img.shape[1]-1 if x > self.Img.shape[1]-1 else x

                N.append((y,x))


            for n in N:
                
                if self.stopCondition2(visit,n):
                    visit[n] = 1
                    self.coords.append(n)
                    self.stack.append(n)
            
        
        return visit  
                
    def stopCondition(self,visit,n,p):
        '''
            
        '''
        
        per = 1
        reach = per*self.Img.max()/100
        return (((self.Img[n]-self.Img[p])>=0.0 ) or (self.Img[p]==0)) and (visit[n] == 0) and (self.Img[n]<=reach)
    
    def stopCondition2(self,visit,sk):
        '''
            
        '''
        
        ###
        
        anglek = self.angle(self.si,sk)
        dk = ((self.si[0]-sk[0])**2 + (self.si[1]-sk[1])**2)**(0.5)
        
        anglej = self.angle(self.si,self.sj)
        dj = ((self.si[0]-self.sj[0])**2 + (self.si[1]-self.sj[1])**2)**(0.5)
            
      
        
        # correction if angleinf is <0, the anglek are in [0,360)

        anglesup = anglej + self.atol
        angleinf = anglej - self.atol
        
        s1 = (angleinf < 0) and ((anglek > (angleinf%360)) or (anglek<anglesup))
        s2 = (anglesup>=360 and ((anglek > angleinf) or (anglek<(anglesup%360))))
        s3 = (anglek<anglesup) and (anglek > angleinf)

        
#        return (s1 or s2 or s3) and (visit[sk] == 0)
        return ( (dk<=(dj+self.dtol)) and (dk>=(dj-self.dtol)) ) and (s1 or s2 or s3) and (visit[sk] == 0)
    
    def angle(self,si,sj):
        
        import numpy as np

    
        s0 = (0,0)
        
        
   
        
        sja = (sj[0]-si[0],sj[1]-si[1])


        dy =  sja[0]-s0[0]
        dx =  sja[1]-s0[1]


        if dx!=0 and dy!=0:
            m1= dy/dx
            if dy > 0 and dx > 0:

                r = np.rad2deg(np.arctan(m1))

            if dy > 0 and dx < 0:

                r = 180+np.rad2deg(np.arctan(m1))

            if dy<=0 and dx < 0:

                r = np.rad2deg(np.arctan(m1))+180

            if dy<=0 and dx > 0:

                r = 360+np.rad2deg(np.arctan(m1))


        if dy == 0 and dx >0:
            r = 90*0

        if dy>0 and dx == 0:
            r = 90*1

        if dy == 0 and dx < 0:
            r = 90*2

        if dx == 0 and dy < 0:
            r = 90*3

        if dx == 0 and dy == 0:
            r = 0

        return r
    
    def makeSemivarianceMap(self,si):
        '''
        
        
        '''
        import numpy as np
        z = np.zeros(self.Img.shape)

        for i in range(self.Img.shape[0]):
            for j in range(self.Img.shape[1]):
                sj = (i,j)
                seed = sj
                set_ri = self.getRegion(seed)
                z[sj] = np.mean([(self.Img[si]-self.Img[c])**2   for c in self.coords])/2
        z[si] = 0
        return z
    
    def getLine(self,sj):
    
        import numpy as np
        line  = list([])
        
        si = self.si
        self.sj = sj
        
        s0 = si
        line.append(s0)
        
        
        d = ((si[0]-sj[0])**2 + (si[1]-sj[1])**2)**(0.5) 
        
        visit = np.zeros(self.Img.shape)
        
  
        
        
        while(d!=0):
        
            visit[s0] = 1
            #l.append(si)
        
        
            pr = it.product(self.rng,self.rng)
            p = si

            dmin = 10000000
            for px in pr:

                y = p[0]+px[0]
                y = 0 if y < 0 else y
                y = self.Img.shape[0]-1 if y > self.Img.shape[0]-1 else y

                x = p[1]+px[1]
                x = 0 if x < 0 else x
                x = self.Img.shape[1]-1 if x > self.Img.shape[1]-1 else x
            
                N.append((y,x))
            
            
            l.append(pmin)
        
        
        
        
        
    
    
        
        
        
        