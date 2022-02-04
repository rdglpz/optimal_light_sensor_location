

def aptitude(X,nonsat,varmask,only_aptitude=True):
    
    #M es rs 
    M = mapa(X, nonsat, varmask)
    
    #creamos n mapas binarios donde la cobertura de los sensores no se traslape
    # sumamos todos los valores mayores que ceros de una coordenada.
    # coleccionamos todos esos valores en un mapa
    # tomamos en cuenta donde solo tenemos un 
    mask = np.sum(M>0,axis=0)==1
    
    #mask = -np.sum(np.max(M,axis=0))
    
 #   np.sum([mask]*len(M)*M)

    return -(np.sum([mask]*len(M)*M))


def validate_coordinates(iy,ix,coords):
    iy = np.where(coords[:,0]==iy)
    ix = np.where(coords[:,1]==ix)

    ix = np.intersect1d(ix,iy)
    return ix

def mapa(X,nonsat,varmask,W1):
    
    
    n_sensors = int(len(X)/2)
    
    sensor_list = X.reshape(n_sensors,len(nonsat.shape))
    
    coverage = np.zeros((n_sensors,nonsat.shape[0],nonsat.shape[1]))
    
    for i,s in enumerate(sensor_list):

        sy = s[0]
        sx = s[1]
        #sy,sx = s[0],s[1]
        ix = validate_coordinates(sy,sx)
        
        coverage[i] = np.zeros(nonsat.shape)
 
        if len(ix)>0: 
            map0to1 = 1/(1+varmask[ix]) 
            coverage[i] = map0to1*(map0to1<1)*W1

 
    return coverage   