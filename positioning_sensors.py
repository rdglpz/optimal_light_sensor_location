import numpy as np

def f5(i,p,e=4):
    '''
    Pesos exponenciales
    
    W = i*a^(p-r)
    '''
    return i*e**(p)

def getN(W,i,j):
    
    ln = ([])
    
    c = W[i,j]
    
    ln.append(W[i+1,j])
    ln.append(W[i-1,j])
    ln.append(W[i-1,j-1])
    ln.append(W[i-1,j+1])
    ln.append(W[i+1,j+1])
    ln.append(W[i+1,j-1])
    ln.append(W[i,j+1])
    ln.append(W[i,j-1])
    
    opt = False if W[i,j]<max(ln) else True
    
    return opt

def fuerza_puntual(x,y,im2,h=2,e=2):
    '''
    this does not make sense if we have corrected verions solving the saturarion image problem
    
    '''
    s = 0.0 
    g1 = float(im2[x][y])
    xlb = 0 if(x-h)<0 else x-h
    xub = im2.shape[0] if(x+h)> im2.shape[0] else x+h
    
    ylb = 0 if(y-h)<0 else y-h
    yub = im2.shape[1] if(y+h)> (im2.shape[1]) else y+h 
    
    ri = np.arange(xlb,xub)
    rj = np.arange(ylb,yub)
    for i in ri:
        for j in rj:
            d=((i-x)**2+(j-y)**2)**(0.5)
            g2 = float(im2[i][j])
            if d >0: s += (g1*g2)/(d**e)
    return s

def filtro_fuerza_puntual(A,h=3,e=2):
    pot1 = np.zeros(A.shape)
    for i in range(pot1.shape[0]):
        for j in range(pot1.shape[1]):
            pot1[i][j] = fuerza_puntual(i,j,A,h,e)
            
    return pot1


def plot_examples(data,colormaps,filename):
    """
    Helper function to plot data with associated colormap.
    """
    n = len(colormaps)
    #fig, axs = plt.subplots(1, n, figsize=(n * 2 + 2, 3),
    #                        constrained_layout=False, squeeze=False)
    fig, axs = plt.subplots(1, n, figsize=(6, 5.2),constrained_layout=False, squeeze=False)
   
    scalebar = ScaleBar(1000) # 1 pixel = 0.2 meter
    plt.gca().add_artist(scalebar)
    for [ax, cmap] in zip(axs.flat, colormaps):
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        psm = ax.pcolormesh(np.flip(data,0), cmap=cmap, rasterized=False, vmin=0, vmax=data.max())
        fig.colorbar(psm, ax=ax)
    
    

    plt.savefig(filename, dpi = 150)

def getMax(img):
    
    B = np.zeros(img.shape)
    
    for i in range(1,B.shape[0]-1):
        for j in range(1,B.shape[1]-1):
            B[i][j] = getN(img,i,j)
    
    positionsB = B*(img>0)
    locationsP = np.zeros(nonsat.shape)
    for i in range(0,nonsat.shape[0]):
        for j in range(0,nonsat.shape[1]):
            locationsP[i][j] = 8 if positionsB[i][j]==1 else P[i,j]
    
    
    return positionsB

def get_coords(seeds):
    coords = np.array([[0,0]])
    for i in range(0,seeds.shape[0]):
        for j in range(0,seeds.shape[1]):
            if seeds[i,j]!=0:
                coords = np.append(coords,[[i,j]],axis=0)
    return coords

def makeSpatialScatter(pc,mask,img):

    i = j = 0
    
    p = np.array([i,j])
    d = np.linalg.norm(pc-pc)
    v = (img[pc[0],pc[1]]-img[pc[0],pc[1]])**2

    ls = [[d,v]]
    
    for i in range(1,img.shape[0]):
        for j in range(1,img.shape[1]):
            if mask[i][j]==1:
                p = np.array([i,j])

                d = np.linalg.norm(p-pc)
                
                v = (img[p[0],p[1]]-img[pc[0],pc[1]])**2
        
                ls = np.append(ls,[[d,v]],axis = 0)
 
  
                
    return ls

def variogram(sc,h=30):
    variogram = np.array([])
    for h in range(1,30):


        ix = (sc[:,0]>(h-1))*(sc[:,0]<=h)
        n = np.sum(ix)
        if n >0:
            v = np.sum(sc[ix,1])
            g = v/(2*n)
            variogram = np.append(variogram,g)
        else:
            variogram = np.append(variogram,0)
    
    return variogram


def getMax(img):
    B = np.zeros(img.shape)
    for i in range(1,B.shape[0]-1):
        for j in range(1,B.shape[1]-1):
            B[i][j] = getN(img,i,j)
    
    positionsB = B*(img>0)
    locationsP = np.zeros(img.shape)
    for i in range(0,img.shape[0]):
        for j in range(0,img.shape[1]):
            locationsP[i][j] = 8 if positionsB[i][j]==1 else 0
            
    return locationsP
