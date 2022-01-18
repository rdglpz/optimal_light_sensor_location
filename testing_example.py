#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 19:42:22 2022

@author: rodrigo
"""


import importlib

from PIL import Image
import numpy as np

from noisyopt import minimizeCompass

from sortedcontainers import SortedList

import pandas as pd

import matplotlib.pyplot as plt


import matplotlib.cbook as cbook
from matplotlib.pylab import rcParams
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm
from matplotlib.pyplot import figure
import matplotlib.patches as mpatches

from skimage.feature import peak_local_max
from skimage.segmentation import watershed

from scipy import ndimage as ndi

plt.rcParams["figure.figsize"] = (6.0*1.5,5.2*1.5)

from skimage import data

from skimage.util import img_as_ubyte

orig_phantom = img_as_ubyte(data.shepp_logan_phantom())

from skimage.morphology import (erosion, dilation, opening, closing,  # noqa
                                white_tophat)
from skimage.morphology import black_tophat, skeletonize, convex_hull_image  # noqa
from skimage.morphology import disk  # noqa

footprint = disk(3)
footprint = np.ones((2,2))

import importlib

import positioning_sensors as ps
importlib.reload(ps)


# funciones

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

def desaturate(img,th=62):
    image = img>=th
    distance = ndi.distance_transform_edt(image)
    nonsat = img+(distance)
    return nonsat


def computeRegions(S,coords,th = 0.6, atol=30, direction_delta = 2):
    """
    S: is the importante map
    
    """
    
    S_flatten = S.flatten()
    vmax = np.max(S_flatten)
    vmin = np.min(S_flatten)
    
    #100% of variability reference
    vrange = (vmax-vmin)**2/2
    max_var = th * vrange/100

    rs=30
    
    
    setC = np.zeros((len(coords),S.shape[0],S.shape[1]))
    z    = np.zeros((len(coords),S.shape[0],S.shape[1]))
    mz   = np.zeros((len(coords),S.shape[0],S.shape[1]))



    for i,c in enumerate(coords):
        print("Coords ", c)
        print(i/len(coords)*100,"%")


    
        dy = c[0]
        dx = c[1]

        
        
        for direction in range(0,360,direction_delta):
#            print("d : ",direction)

            setA = np.zeros(S.shape)
            setB = np.zeros(S.shape)

            gm1 = direction-atol
            gm2 = direction+atol

            m1 = np.tan(np.radians(gm1))
            m2 = np.tan(np.radians(gm2))

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

 #           setB = generateMaks(dims,p,theta,theta_range)
            mask = setA*setB
            p = np.array([dy,dx])
            
            p = np.array([dy,dx])

            sc = ps.makeSpatialScatter(p,mask,S)
            print(sc,"--",direction)
            print(type(sc))
            accum = ps.variogram(sc)
            
            print(accum,"-accum")

            # detectamos hasta que indice se cumple el requerimiento de la varianza

            bs = accum <= max_var
            nf = np.where(bs==False)
            
   
                
            if len(nf[0])>0:

                bs[nf[0][0]:] = False
                radio = np.sum(bs)

            else:

                radio = len(bs)
            
            if True:
#                print("enter",radio)
                inspect_bs = bs
                inspect_nf = nf
                inspect_radio = radio


                for r in np.arange(0,radio,0.1):
                    

                    y = np.int32(np.round(r*np.sin(np.radians(direction))))
                    x = np.int32(np.round(r*np.cos(np.radians(direction))))

                    #we take care of the positive squared boundaries
                    
                    py = dy+y
                    px = dx+x
                    if py>=nonsat.shape[0]:
                        py = nonsat.shape[0]-1
                    elif py < 0 :
                        py = 0
                        
                    if px>=nonsat.shape[1]:
                        px = nonsat.shape[1]-1
                    elif px < 0 :
                        px = 0
                    
                        
                        
                    #py = nonsat.shape[0]-1 if dy+y>=nonsat.shape[0] else dy+y
                    #py = 0 if dy+y<0 else dy+y
                    #px = nonsat.shape[1]-1 if dx+x>=nonsat.shape[1] else dx+x
                    #px = 0 if dx+x<0 else dx+x
                    
                    
                    
                    z[i][py,px] = accum[int(r)]
                    mz[i][py,px]=1



        setC[i][c[0]][c[1]]=1
        mz[i][c[0]][c[1]]=1
        
        print("--")
    
    return z,mz,setC

#Original nocturnal image without values <21
ilumina = 'example/inputs/qro_light_th'
niveles = 'example/inputs/prioridades'

luminance = readIMG(ilumina)
importance = readIMG(niveles,invert=True)
nonsat = desaturate(luminance,th=62)
W1 = ps.f5(nonsat,importance)
W2 = ps.f2(nonsat,importance)
seeds = ps.getMax(W1)
seeds_nonsat = ps.getMax(nonsat)
coords = ps.get_coords(seeds)
coords_nonsat = ps.get_coords(seeds)

variogram, mask ,locations = computeRegions(W1,[coords[-2]],th = 0.25 ,atol=15,direction_delta = 2)
    