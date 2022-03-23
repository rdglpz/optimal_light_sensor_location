import itertools as it


class rGrowing():
    
    
    
    
    
    def __init__(self,Img):
        '''
        
        '''

        self.Img = Img
        self.rng = list(range(-1,2))
    
    
    def getRegion(self,t):
        '''
        
        '''
        
        print(t)
    #    Img = self.Img
        self.stack = []
        
        visit = np.zeros(self.Img.shape)
        visit[t] = 1
        self.stack.append(t)
        
        while True if self.stack else False:




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
           

                if self.stopCondition(visit,n,p):
                    visit[n] = 1
                    self.stack.append(n)
            
        
        return visit  
                
    def stopCondition(self,visit,n,p):
        '''
            
        '''
        per = 1
        reach = per*self.Img.max()/100
        return (((self.Img[n]-self.Img[p])>=0.0 ) or (self.Img[p]==0)) and (visit[n] == 0) and (self.Img[n]<=reach)
        
        
        
        