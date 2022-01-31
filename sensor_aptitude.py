def aptitude(X,,W1=W1 = nonsat=nonsat,varmask=varmask,only_aptitude=True):
    
    
    C = mapa(X)

 
    return -np.sum(res>0)


def validate_coordinates(iy,ix,coords=coords):
    iy = np.where(coords[:,0]==iy)
    ix = np.where(coords[:,1]==ix)

    ix = np.intersect1d(ix,iy)
    return ix

def mapa(X,nonsat=nonsat,varmask=varmask,only_aptitude=True):
    
    
    n_sensors = int(len(X)/2)
    sensor_list = X.reshape(n_sensors,len(varmask.shape))
    

    coverage = np.zeros((n_sensors,varmask.shape[0],varmask.shape[1]))
    
    for i,s in enumerate(sensor_list):

        sy = s[0]
        sx = s[1]
        ix = validate_coordinates(sy,sx)
        
        coverage[i] = np.zeros(varmask.shape)
 
        if len(ix)>0: 
            map0to1 = 1/(1+varmask[ix]) 
            coverage[i] = map0to1*(map0to1<1)*W1

    res = np.zeros(nonsat.shape)
    
    for s in coverage:
        res = res+s*(nonsat>0)*W1

 
    return coverage