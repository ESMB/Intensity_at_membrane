#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 12:57:38 2021

@author: Mathew
"""


from skimage.io import imread
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from skimage.segmentation import find_boundaries
from skimage.filters import threshold_local
import napari 
from PIL import Image


filename_contains="cropped.tif"

# Folders to analyse:


path="/Users/Mathew/Documents/Current analysis/Images-for-quant-memb-fluorescence/DL vs LIVE-PAINT images/7/"


def load_image(toload):
        
    image=imread(toload)
    
    return image


def subtract_bg(image):
    background = threshold_local(image, 11, offset=np.percentile(image, 1), method='median')
    bg_corrected =image - background
    return bg_corrected



for root, dirs, files in os.walk(path):
        for name in files:
                if filename_contains in name:
                    if ".tif" in name:
                        resultsname = name
                        print(resultsname)

image_path=path+resultsname
print(image_path)
    
    
# Load the image
cell_image=load_image(image_path)

# Label the cells- requires filling in cells.
cells_in=np.zeros(np.shape(cell_image),dtype=int)


viewer = napari.Viewer()
viewer.add_image(cell_image,rgb=False)
viewer.add_labels(cells_in,name='Cells')
cells_label= viewer.layers["Cells"].data
napari.run()

# while cells_label.max() < 1:
#     i=1



def analysis():

    # Width of "peripjhery"
    width_to_look_at=3
    
    # Find boundary.
    bound=find_boundaries(cells_label, mode='thick').astype(np.uint8)
    
    # Thicken the boundary
    for i in range(0,width_to_look_at):
        
        bound_new=find_boundaries(bound, mode='thick').astype(np.uint8)
        
        bound=bound+bound_new
    
    bound=(bound)>0
    boundary=bound*cells_label
    centre=cells_label-boundary
    
    boundary_image=(boundary>0)*cell_image
    centre_image=(centre>0)*cell_image
    
    fig,ax=plt.subplots(1,2)
    ax[0].imshow(boundary_image)
    ax[1].imshow(centre_image)
    plt.savefig(path+"cells.tif")
    plt.show()
    
    ims = Image.fromarray(boundary_image)
    ims.save(path+'boundary_image.tif')
    
    ims = Image.fromarray(centre_image)
    ims.save(path+'centre_image.tif')
    
    
    filtered=subtract_bg(cell_image)
    
    
    number_of_cells=cells_label.max()
    
    boundary_mean=[]
    centre_mean=[]

    
    Output_all_cells = pd.DataFrame(columns=['Boundary','Centre','Ratio'])
    
    
    for i in range(0,number_of_cells):
        print(i+1)
        
        boundary_cell=boundary==(i+1)
        centre_cell=centre==(i+1)
    
        boundary_list=(boundary_cell*filtered).flatten()
        boundary_list_values=boundary_list[boundary_list!=0]
    
        centre_list=(centre_cell*filtered).flatten()
        centre_list_values=centre_list[centre_list!=0]
    
        plt.hist(centre_list_values, bins = 50,range=[0,400], rwidth=0.9,color='#ff0000')
        plt.xlabel('Intensity (ADU)',size=20)
        plt.ylabel('Number of pixels',size=20)
        plt.title('Centre',size=20)
        plt.savefig(path+"Cell_"+str(i+1)+"_Centre.pdf")
        plt.show()
        
        plt.hist(boundary_list_values, bins = 50,range=[0,400], rwidth=0.9,color='#ff0000')
        plt.xlabel('Intensity (ADU)',size=20)
        plt.ylabel('Number of pixels',size=20)
        plt.title('Boundary',size=20)
        plt.savefig(path+"Cell_"+str(i+1)+"_Boundary.pdf")
        plt.show()
    
        mean_centre=centre_list_values.mean()
        mean_boundary=boundary_list_values.mean()
    
        boundary_mean.append(mean_boundary)
        centre_mean.append(mean_centre)
        ratio=mean_boundary/mean_centre
        
        Output_all_cells = Output_all_cells.append({'Boundary':mean_boundary,'Centre':mean_centre,'Ratio':ratio},ignore_index=True)
    
        Output_all_cells.to_csv(path + 'all_metrics.csv', sep = '\t')
