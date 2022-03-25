import numpy as np
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
from scipy import ndimage as ndi
from skimage.segmentation import watershed
from PIL import Image
from IPython.display import display, clear_output

def f2(i,p):
    '''
    Asignación de pesos posicional
    
    La caraceterística mas significativa se coloca 
    a la izquierda y es multiplicada por una 
    base compuesta definida por número total de valores que 
    puede tener la caracteristica inmediata colocada una
    posición a la derecha.
    
    W(P,L)=P*|L|+L
    
    W = P*|L|+L
    '''
    L = np.max(i)+1
    return p*L+i

def f5(i,p,e=4):
    '''
    Pesos exponenciales
    
    W = i*a^(p-r)
    '''
    print(" W = i*a^(p)")
    return i*e**(p)

def f6(i,s,p=4):
    '''
    Pesos exponenciales
    
    W = i*s^(p)
    '''
    print(" W = i*s^(p)")
    return i*s**p

def getN(W,i,j):
    
    ln = ([])
    
#    c = W[i,j]
    
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
    

def get_coords(seeds):
    coords = np.array([[0,0]])
    for i in range(0,seeds.shape[0]):
        for j in range(0,seeds.shape[1]):
            if seeds[i,j]!=0:
                coords = np.append(coords,[[i,j]],axis=0)
                
        
    return coords[1:]

def makeSpatialScatter(pc,mask,img):
    """
    pc: the centroid
    mask: pixels of interet, 
    img nlt image
    
    returns [distance,variance]
    """

    i = j = 0
    
    p = np.array([i,j])
    d = np.linalg.norm(pc-pc)
    v = (img[pc[0],pc[1]]-img[pc[0],pc[1]])**2

    ls = np.array([[d,v]])
    
    for i in range(1,img.shape[0]):
        for j in range(1,img.shape[1]):
            if mask[i][j]==1:
                p = np.array([i,j])

                d = np.linalg.norm(p-pc)
                
                v = (img[p[0],p[1]]-img[pc[0],pc[1]])**2
        
                ls = np.append(ls,[[d,v]],axis = 0)
 
  
                
    return ls

def variogram(sc,h=30):
    """
    sc: spatial scatter
    
    """
    if np.array(sc).shape[0] > 0:
    
        variogram = np.array([])
        for h in range(1,h):
    
       
            ix = (sc[:,0]>(h-0.5))*(sc[:,0]<=(h+0.5))
            n = np.sum(ix)

            
            if n > 0:
                v = np.sum(sc[ix,1])
                g = v/(2*n)
                variogram = np.append(variogram,g)
            else:
                variogram = np.append(variogram,0)
    else:
        variogram = range(1,h)
    
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
            locationsP[i][j] = 1 if positionsB[i][j]==1 else 0
            
    return locationsP

def waterShedRegions(S,C):
    """
    S: Seeds
    C: Coordinates
    """
    
    mask = np.zeros(S.shape, dtype=bool)
    #we set to 1 all the coordinates C in the mask
    mask[tuple(C.T)] = True

    markers, _ = ndi.label(mask)
    labels = watershed(-S, markers, mask=S)
    
    L = np.zeros((len(C),S.shape[0],S.shape[1]))
    
    for i in range(len(L)):
        c = C[i]
        L[i] = (labels == labels[c[0]][c[1]])
    
    return L

def getMaximumVariance(S,th):
    """
    th ={0,100}
    """
    S_flatten = S.flatten()
    vmax = np.max(S_flatten)
    vmin = np.min(S_flatten)
    vrange = (vmax-vmin)**2/2
    max_var = th * vrange/100
    return max_var

def getOrientationMaskII(si,sj,atol,Img):
    """
    atol: angle tolerance
    """
    
    # area constrained by angle
    A = np.zeros(Img.shape)
    angle = ps.angle(si,sj)
    angleinf = (angle-atol)
    anglesup = (angle+atol)
    
    for i in range(nonsat.shape[0]):
        for j in range(nonsat.shape[1]):
            
            #sk is the test coordinate. is the sk point inside the cone?
            sk = np.array([i,j])
            
            
            #probar este angulo si esta en el cono
            # b va de 0 a 359.999 
            #is the angle of interest?
            anglek = ps.angle(si,sk)
            
            # correction if angleinf is <0, the anglek are in [0,360)
            if angleinf < 0 and ((anglek > (angleinf%360)) or (anglek<anglesup)) :
                A[i][j] = 1
                    
            if anglesup>=360 and ((anglek > angleinf) or (anglek<(anglesup%360))) :
#                anglesupp = anglesup%360
#                if (anglek > angleinf) or (anglek<(anglesup%360)):

                A[i][j] = 1          
            if (anglek<anglesup) and (anglek>angleinf):
                A[i][j] = 1
    
    return A

def getOrientationMask(d,dtol,S,dx,dy):
    """
    d: direction in grads (0-360)
    dtol: direction tolerance (range) in angle
    S: is a matrix
    """

    gm1 = d-dtol
    gm2 = d+dtol
    m1 = np.tan(np.radians(gm1))
    m2 = np.tan(np.radians(gm2))
    
    print(m1,m2)
    
    setA = np.zeros(S.shape)
    setB = np.zeros(S.shape)
    
    for y in range(S.shape[0]):
        for x in range(S.shape[1]):

            #if angle is in right or left side
            if gm1<90 or gm1>90*3:
                setA[y][x] = 1 if m1*(x-dx)+dy-y <0 else 0
            else:
                setA[y][x] = 0 if m1*(x-dx)+dy-y <0 else 1

            if gm2<90 or gm2>90*3:
                setB[y][x] = 0 if m2*(x-dx)+dy-y <0 else 1
            else:
                setB[y][x] = 1 if m2*(x-dx)+dy-y < 0 else 0

    return setA*setB

def angle(si,sj):
    
    s0 = si-si
    sja = sj-si
#    print(s0)
#    print(sja)
    
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
    


    
def getMaxRadio(ac,mv):
    """
    Obtenemos el radio maximo que satisface la maxima varianza (mv) permitida
    dado el acumulado de la varianza en funcion de la distancia
    """
    
    bs = ac <= mv
     
    
    
    nf = np.where(bs==False)

    
    if len(nf[0])>0:
        bs[nf[0][0]:] = False
        radio = np.sum(bs)
    else:
        radio = len(bs)
        
    return radio
                
                

def fillArea(p,accum,radio,direction,mshape):
    emptyZ = np.zeros(mshape)
    emptyMZ = np.zeros(mshape)
    
    dy = p[0]
    dx = p[1]
    
    for r in np.arange(0,radio,0.1):
                    
        #asumiendo que el angulo está centrado en el origen calculamos las coordenadas (y,x) dado la direccon dado por el ángulo en radianes y el tamaño del radio
        y = np.int(np.round(r*np.sin(np.radians(direction))))
        x = np.int(np.round(r*np.cos(np.radians(direction))))

        #we take care of the positive squared boundaries
        #trasladamos las coordenadas a su posición original            
        py = dy+y
        px = dx+x
                
        #validamos que (py,px) estén dentro de los limites de la matriz        
        if py>=mshape[0]:
            py = mshape[0]-1
        elif py < 0 :
            py = 0
                        
        if px>=mshape[1]:
            px = mshape[1]-1
        elif px < 0 :
            px = 0
                    
    
        emptyZ[py,px] = accum[int(r)]
        emptyMZ[py,px] = 1
        
    return emptyZ,emptyMZ
    
    

def computeRegions(S,coords,th = 0.6, atol=30, direction_delta = 2,verbose=False):
    """
    S: NLTI  map
    coords: maximum points
    th: valid variannce in %
    atol = angle tolerance 
    direction_delta: angle steps
    verbose: show prints
    """
    
    max_var = getMaximumVariance(S,th)

    setC = np.zeros((len(coords),S.shape[0],S.shape[1]))
    
    z    = np.zeros((len(coords),S.shape[0],S.shape[1]))
    mz   = np.zeros((len(coords),S.shape[0],S.shape[1]))
    
    for i,c in enumerate(coords):
        clear_output(wait=True)
        display( "{:.2f}%".format(100*(i/len(coords))) )
        if verbose == True: print("Coords ", c)
    
        dy = c[0]
        dx = c[1]

        for direction in range(0,360,direction_delta):
            
        
            mask = getOrientationMask(direction,atol,S,dx,dy)

            p = np.array([dy,dx])
            
            sc = makeSpatialScatter(p,mask,S)
            accum = variogram(sc)
            

            # detectamos hasta que indice se cumple el requerimiento de la varianza
            #getVar
            
            # lo primero que se cumpla, rebasa la toleracia de la varianza o encuentre un maximo
            accum = np.append([0],accum)
            radio = getMaxRadio(accum,max_var)  
            
            
            args = np.where((accum[1:]-accum[:-1])<=0)
            radio_first_local_max = radio+1
            
            if len(args[0])>0:
                radio_first_local_max = args[0][0]+1
            
            #print(direction,radio,radio_first_local_max)
            radio = np.min([radio,radio_first_local_max])
           
            
            
            z_aux    = np.zeros((S.shape[0],S.shape[1]))
            mz_aux   = np.zeros((S.shape[0],S.shape[1]))
            
            a_aux, mz_aux = fillArea(p,accum,radio,direction,S.shape)
            
            xor_mask = np.logical_xor(a_aux,mz[i])

            mz[i]+= xor_mask*a_aux
            z[i] += z_aux
            
            #mz[i] += mz_aux
        
        setC[i][c[0]][c[1]]=1
        
        #use a mask to point out that the 0 is for 0 variance associated to the sensor locations
        mz[i][c[0]][c[1]]=0.000000000001
        
        if verbose == True: print("--")
    
    return z,mz,setC
    

    

def readIMG(img,invert=False,null=255):
    im1 = np.array(Image.open(img))
    if invert == False:
        im1 = np.array(Image.open(img))
        im1 = np.where(im1==null, 0, im1) 
    #    print("categories:", set(im1.flatten()))
    else:
        
        nc = 5
        P = np.where(np.isnan(im1),nc, im1)-1 
        im1 = P.max()-P

    return im1

def plotMasks(mask,L,W):
    ngrid = np.int32(np.ceil(np.sqrt(len(mask))))
    fig, axs = plt.subplots(ngrid, ngrid,figsize = (30,30))
    c=0
    for i in range(ngrid):
        for j in range(ngrid):
            if i * ngrid + j < len(mask):
                axs[i, j].imshow(mask[c])
                axs[i, j].set_title("w: {:.1f}, c {}".format(np.sum(mask[c]*W*L[c]),str(coords[c])))
                c+=1
    plt.show()
    
    
def desaturate(img,th=62):
    image = img>=th
    distance = ndi.distance_transform_edt(image)
    nonsat = img+(distance)
    return nonsat


def saveRegions(varmask,locations,name = "allcoversnonsatat15percent.csv"):
    
    #flattenizing

    fvm = [varmask[i].flatten() for i in range(len(varmask))]
    df = pd.DataFrame(fvm)
    
    c  = np.array([ps.get_coords(locations[i]) for i in range(len(locations)) ])
    cf = c.flatten()
    coords = cf.reshape(len(locations),2)
    df.insert(0,"coordsy", coords[:,0])
    df.insert(1,"coordsx", coords[:,1])
    df.to_csv(name)
    
#def readRegions(name):


def gom(si,sj,atol,Img):
    """
    gom: get orienttion mask
    atol: angle tolerance
    """
    
    # area constrained by angle
    A = np.zeros(Img.shape)
    ang = angle(si,sj)
    angleinf = (ang-atol)
    anglesup = (ang+atol)
    
    # en vez de recorrer toda la imagen, usar crecimiento de regiones inicializado en si
    # considerar si los vecinos sk;
    # 1) están dentro de angle inf y angle sup
    for i in range(Img.shape[0]):
        for j in range(Img.shape[1]):
            
            #sk is the test coordinate. is the sk point inside the cone?
            sk = np.array([i,j])
            
            
            #probar este angulo si esta en el cono
            # b va de 0 a 359.999 
            #is the angle of interest?
            anglek = angle(si,sk)
            
            # correction if angleinf is <0, the anglek are in [0,360)
            
            s1 =  (angleinf < 0) and ((anglek > (angleinf%360)) or (anglek<anglesup))
            s2 = (anglesup>=360 and ((anglek > angleinf) or (anglek<(anglesup%360))))
            s3 = (anglek<anglesup) and (anglek>angleinf)
            
            
#            if (angleinf < 0) and ((anglek > (angleinf%360)) or (anglek<anglesup)) :
#                A[i][j] = 1
                    
#            if anglesup>=360 and ((anglek > angleinf) or (anglek<(anglesup%360))) :
#                A[i][j] = 1   
                
#            if (anglek<anglesup) and (anglek>angleinf):
#                A[i][j] = 1
                

            #if (s1 or s2 or s3) and (d(sk,si) inside the range): 
            #    A[i][j] = 1
            
            if (s1 or s2 or s3):
                A[i][j] = 1
    
    return A
    
    
    
    