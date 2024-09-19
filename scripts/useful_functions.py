# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import math
import numpy as np
#import time, sys
from IPython.display import clear_output
import matplotlib.pyplot as plt
from copy import deepcopy
#from vispy.color import Colormap
from IPython.display import display
from matplotlib import colors
from matplotlib.colors import to_hex
#from label_comparison import *

import pandas as pd

import dask.array as da
from dask import delayed, array
#import ipywidgets

#import ipywidgets
from tqdm import tqdm, trange
from copy import deepcopy
#from numba import njit
import itertools

def update_progress(progress):
    bar_length = 20
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1

    block = int(round(bar_length * progress))

    clear_output(wait = True)
    text = "Progress: [{0}] {1:.1f}%".format( "#" * block + "-" * (bar_length - block), progress * 100)
    print(text)


def update_progress_4dp(progress):
    bar_length = 20
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1

    block = int(round(bar_length * progress))

    clear_output(wait = True)
    text = "Progress: [{0}] {1:.4f}%".format( "#" * block + "-" * (bar_length - block), progress * 100)
    print(text)


def get_block_indices(xarraysize, yarraysize, myblocksize):
    #I think this function may not go all the way up to the edges
    N_x = math.ceil(xarraysize/myblocksize)
    N_y = math.ceil(yarraysize/myblocksize)
    N_total = N_x * N_y
    
    indices = np.zeros([N_total, 5], dtype=int)
    
    for k in range(0, N_total):
        indices[k,0]=k

    xstartvalue = int(0)
    ystartvalue = int(0)
    xfinishvalue = int(myblocksize)
    yfinishvalue = int(myblocksize)
    
    for i in range(0,N_total):
        
        indices[i,1] = xstartvalue
        indices[i,2] = xfinishvalue
        indices[i,3] = ystartvalue
        indices[i,4] = yfinishvalue
        
        
        xstartvalue += myblocksize
        xfinishvalue += myblocksize
        
        if xfinishvalue > xarraysize:
            xfinishvalue = xarraysize
        
        if xstartvalue > xarraysize:
            xstartvalue = 0
            xfinishvalue = myblocksize
            ystartvalue += myblocksize
            yfinishvalue += myblocksize
            if yfinishvalue > yarraysize:
                yfinishvalue = yarraysize
        
    
    return indices
#
    

def get_block_indices_without_edges(x0arraysize, x1arraysize, myblocksize):
    """
    Produces tile coordinates and ignores any incomplete tiles

    Parameters
    ----------
    x0arraysize : int
        The size of the array that you would like to tile in the 0th direction
    x1arraysize : Tint
        DESCRIPTION.
    myblocksize : TYPE
        DESCRIPTION.

    Returns
    -------
    indices : TYPE
        DESCRIPTION.

    """
    
    N_x0 = math.floor(x0arraysize/myblocksize)
    N_x1 = math.floor(x1arraysize/myblocksize)
    N_total = N_x0 * N_x1
    
    indices = np.zeros([N_total, 5], dtype=int)
    
    for k in range(0, N_total):
        indices[k,0]=k

    x0startvalue = int(0)
    x1startvalue = int(0)
    x0finishvalue = int(myblocksize)
    x1finishvalue = int(myblocksize)
    
    for i in range(0,N_total):
        
        indices[i,1] = x0startvalue
        indices[i,2] = x0finishvalue
        indices[i,3] = x1startvalue
        indices[i,4] = x1finishvalue
        
        
        x0startvalue += myblocksize
        x0finishvalue += myblocksize
        
        if x0finishvalue > x0arraysize:
            x0startvalue = 0
            x0finishvalue = myblocksize
            x1startvalue += myblocksize
            x1finishvalue += myblocksize
    
    return indices


def get_block_indices_with_buffer(xarraysize, yarraysize, myblocksize, buffersize):
    #this includes the edges
    #After this has been outputted, it is good to assign variables like this:
    #i_tile = tile[0] #this is just a counting index of which tile you are on
    #out = tile[1:5] #these are the coordinates of the tile without a buffer
    #buffered_tile = tile[5:9] #these are the coordinates of the tile including a buffer
    #internal_extract = tile[9:13] #this is how to extract the tile from buffered tile
    #stardist tiling is an example of this
    
    N_x = math.ceil(xarraysize/myblocksize)
    N_y = math.ceil(yarraysize/myblocksize)
    N_total = N_x * N_y
    
    indices = np.zeros([N_total, 13], dtype=int)
    
    for k in range(0, N_total):
        indices[k,0]=k

    xstartvalue = int(0)
    ystartvalue = int(0)
    xfinishvalue = int(myblocksize)
    yfinishvalue = int(myblocksize)
    
    for i in range(0,N_total):
        
        indices[i,1] = xstartvalue
        indices[i,2] = xfinishvalue
        indices[i,3] = ystartvalue
        indices[i,4] = yfinishvalue
        
        
        xstartvalue += myblocksize
        xfinishvalue += myblocksize
        
        if xfinishvalue > xarraysize:
            xfinishvalue = xarraysize
        
        if xstartvalue > xarraysize:
            xstartvalue = 0
            xfinishvalue = myblocksize
            ystartvalue += myblocksize
            yfinishvalue += myblocksize
            if yfinishvalue > yarraysize:
                yfinishvalue = yarraysize
        
    #this next bit adds the buffer
    
    for i in range(len(indices)):
        lhs = indices[i, 1]
        rhs = indices[i, 2]
        top = indices[i, 3]
        bottom = indices[i, 4]
        
        if lhs - buffersize < 0:
            indices[i,5]=0
            indices[i,9] = lhs
        else:
            indices[i,5]=lhs-buffersize
            indices[i,9] = buffersize
        
        if top - buffersize < 0:
            indices[i,7]=0
            indices[i,11] = top
        else:
            indices[i,7]=top-buffersize
            indices[i,11] = buffersize
    
        if rhs + buffersize > xarraysize:
            indices[i,6] = xarraysize
            indices[i,10] = indices[i,9] + (xarraysize % myblocksize)
            
        else:
            indices[i,6] = rhs+buffersize
            indices[i,10] = indices[i,9] + myblocksize
            
        if bottom + buffersize > yarraysize:
            indices[i,8] = yarraysize
            indices[i,12] = indices[i,11] + (yarraysize % myblocksize)
        else:
            indices[i,8] = bottom + buffersize
            indices[i,12] = indices[i,11] + myblocksize
        
    return indices


def get_5_channels(thisstring, firstchannel):

    if 'EPCAM_IHC' in thisstring:
        thisstring = thisstring.replace('EPCAM_IHC', 'EPCAMIHC')
        
    if 'EPCAM_(IHC)' in thisstring:
        thisstring = thisstring.replace('EPCAM_(IHC)', 'EPCAMIHC')
        
    if 'EPCAM_IH' in thisstring and not 'EPCAM_IHC' in thisstring:
        thisstring = thisstring.replace('EPCAM_IH', 'EPCAMIHC')
    
    n1 = thisstring.find(firstchannel)
    d1 = thisstring[n1:].find('_')
    n2 = n1+d1+1
    d2 = thisstring[n2:].find('_')
    n3 = n2 + d2 + 1
    d3 = thisstring[n3:].find('_')
    n4 = n3 + d3 + 1
    d4 = thisstring[n4:].find('_')
    n5 = n4 + d4 + 1
    d5 = thisstring[n5:].find('_')
    thisstring[n5:]
    channels = [thisstring[n1:n1+d1], thisstring[n2:n2+d2], thisstring[n3:n3+d3], thisstring[n4:n4+d4], thisstring[n5:n5+d5]]
    return channels

def get_output_single_channel_name(file, firstchannel, k):
    chs = get_5_channels(file, firstchannel)
    for i, eachch in enumerate(chs):
        if i == k:
            continue
        else:
            file = file.replace(eachch + '_', '')
    return file

def get_output_single_channel_name_from_single_filename_list(filelist, preceding):
    filelistout = []
    for eachfile in filelist:
        ix = eachfile.find(preceding) + len(preceding) + 1;
        iy = eachfile[ix:].find('_') + ix;
        name = eachfile[ix:iy]
        filelistout.append(name)
    return filelistout

def upsize(array, scaling, show_progress=True):
    scale = round(scaling)
    x0_test = array.shape[0]*scale
    x1_test = array.shape[1]*scale
    x0 = round(x0_test)
    x1 = round(x1_test)
    if  x0-x0_test > 0.0001:
        print('Array dimensions must divide by scaling')
        exit()
    if  x1-x1_test > 0.0001:
        print('Array dimensions must divide by scaling')
        exit()
        
    if len(array.shape) == 2:
        outarray = np.zeros((x0, x1), dtype = array.dtype)
        
        patch = np.ones((scale, scale), dtype = array.dtype)
        
        for i in range(0, array.shape[0]):
    #        update_progress(i/array.shape[0])
            for j in range(0, array.shape[1]):
                outarray[scale*i:scale*(i+1), scale*j:scale*(j+1)] = array[i,j] * patch
    
    elif len(array.shape) == 3:
        outarray = np.zeros((x0, x1, array.shape[2]), dtype = array.dtype)
        
        patch = np.ones((scale, scale, array.shape[2]), dtype = array.dtype)
        
        for i in range(0, array.shape[0]):
    #        update_progress(i/array.shape[0])
            for j in range(0, array.shape[1]):
                outarray[scale*i:scale*(i+1), scale*j:scale*(j+1)] = array[i,j] * patch
            
    return outarray

def downsize(array, scaling):
    scale = round(scaling)
    x0_test = array.shape[0]/scale
    x1_test = array.shape[1]/scale
    x0 = round(x0_test)
    x1 = round(x1_test)
    if  x0-x0_test > 0.0001:
        print('Array dimensions must divide by scaling')
        exit()
    if  x1-x1_test > 0.0001:
        print('Array dimensions must divide by scaling')
        exit()
    outarray = np.zeros((x0, x1), dtype = array.dtype)

    
    for i in range(0, outarray.shape[0]):
        update_progress(i/outarray.shape[0])
        for j in range(0, outarray.shape[1]):
            outarray[i,j] = array[scale*i,scale*j]
    return outarray

def downsize2(myarray, scaling):
    scale = round(scaling)
    x0_test = myarray.shape[0]/scale
    x1_test = myarray.shape[1]/scale
    x0 = round(x0_test)
    x1 = round(x1_test)
    if  x0-x0_test > 0.0001:
        print('Array dimensions must divide by scaling')
        exit()
    if  x1-x1_test > 0.0001:
        print('Array dimensions must divide by scaling')
        exit()
    outarray = np.zeros((x0, x1), dtype = myarray.dtype)

    outarray = myarray[::scale,::scale]
    return outarray

def downsize_from_middle(myarray, scaling):
    scale = round(scaling)
    x0_test = myarray.shape[0]/scale
    x1_test = myarray.shape[1]/scale
    x0 = round(x0_test)
    x1 = round(x1_test)
    if  x0-x0_test > 0.0001:
        print('Array dimensions must divide by scaling')
        exit()
    if  x1-x1_test > 0.0001:
        print('Array dimensions must divide by scaling')
        exit()
    outarray = np.zeros((x0, x1), dtype = myarray.dtype)

    outarray = myarray[math.floor(scale/2)::scale,
                       math.floor(scale/2)::scale]
    return outarray


def downsize_by_mode(myarray, scale):
    import scipy
    N0_steps = math.floor(myarray.shape[0]/scale)
    N1_steps = math.floor(myarray.shape[1]/scale)
    
    
    if len(myarray.shape) == 2:
        output = np.zeros((N0_steps, N1_steps), dtype = myarray.dtype)
        for i in trange(N0_steps):
            for j in range(N1_steps):
    #            print(scale*i, scale*i+scale, scale*j, scale*j+scale)
                this_array = myarray[scale*i:scale*i+scale, scale*j:scale*j+scale]
                counts = np.bincount(this_array.ravel())
                output[i,j] = np.argmax(counts)
        
    #This doesn't work - need to switch to uniqueness
    elif len(myarray.shape) == 3:
        from joblib import Parallel
        import joblib
        def set_mode_to_output(myarray, i, output):
            output_intermediate = np.zeros((N1_steps, myarray.shape[2]), dtype=myarray.dtype)
            for j in range(N1_steps):
                this_array = myarray[scale*i:scale*i+scale, scale*j:scale*j+scale]
                this_array_reshape = this_array.reshape(-1, myarray.shape[2])
                mode = scipy.stats.mode(this_array_reshape)
                output_intermediate[j] = mode[0]
            return output_intermediate
            
        output = np.zeros((N0_steps, N1_steps, myarray.shape[2]), dtype = myarray.dtype)        
        results = Parallel(n_jobs=20)(joblib.delayed(set_mode_to_output)(myarray, i, output) for i in range(N0_steps))
        #results can be trange
        output = np.array(results)
    
    
#    elif len(myarray.shape) == 3:
#        output = np.zeros((N0_steps, N1_steps, myarray.shape[2]), dtype = myarray.dtype)
#        for i in trange(N0_steps):
#            for j in range(N1_steps):
#                this_array = myarray[scale*i:scale*i+scale, scale*j:scale*j+scale]
#                this_array_reshape = this_array.reshape(-1, myarray.shape[2])
#                mode = scipy.stats.mode(this_array_reshape)
#                output[i,j] = mode[0]
    return output


def crop_to_superpixel(array, pixelsize):
    """
    get_cropped_superpixel_dims is another useful function like this
    """
    x0 = array.shape[0] - array.shape[0]%int(pixelsize)
    x1 = array.shape[1] - array.shape[1]%int(pixelsize)
    outarray = array[0:x0, 0:x1]
    return outarray

def getvaluefromstring(folder, variable):
    i = folder.index(variable)
    length = len(variable)
    i_ans_start = i + length + 1
    try:
        i_ans_end = folder[i_ans_start:].index('_') + i_ans_start
    except:
        i_ans_end = len(folder)
    return folder[i_ans_start:i_ans_end]

def getvaluefromstringslashend(folder, variable):
    i = folder.index(variable)
    length = len(variable)
    i_ans_start = i + length + 1
    try:
        i_ans_end = folder[i_ans_start:].index('/') + i_ans_start
    except:
        i_ans_end = len(folder)
    return folder[i_ans_start:i_ans_end]

def getvaluefromstringdotend(folder, variable):
    i = folder.index(variable)
    length = len(variable)
    i_ans_start = i + length + 1
    try:
        i_ans_end = folder[i_ans_start:].index('.') + i_ans_start
    except:
        i_ans_end = len(folder)
    return folder[i_ans_start:i_ans_end]

def getvaluefromstringdotendnoprecharacter(folder, variable):
    i = folder.index(variable)
    length = len(variable)
    i_ans_start = i + length
    try:
        i_ans_end = folder[i_ans_start:].index('.') + i_ans_start
    except:
        i_ans_end = len(folder)
    return folder[i_ans_start:i_ans_end]

def get_list_from_string(mystring, variable, delimiter):
    #chooses the list in a string that is separated by a delimiter and ended by next dot slash or underscore
    i_start = mystring.find(variable) + len(variable) + 1
    mystring_afterstart = mystring[i_start:]
    
    list_of_possible_stop_char = ['.', '/', '_']
    dist_of_next_chars = []
    for each in list_of_possible_stop_char:
        if each in mystring_afterstart:
            dist_of_next_chars.append(mystring_afterstart.find(each))
    
    if len(dist_of_next_chars) != 0:
        i_end = min(dist_of_next_chars) + i_start
        liststring = mystring[i_start:i_end]
    elif len(dist_of_next_chars) == 0:
        liststring = mystring[i_start:]
    
    return liststring.split(delimiter)

def getvaluefromstringbest(folder, variable, preceding='_', ending='_', mydtype=str):
    i = folder.index(variable)
    length = len(variable)
    i_ans_start = i + length + len(preceding)
    
    if ending in folder[i_ans_start:]:
        i_ans_end = folder[i_ans_start:].index(ending) + i_ans_start
    else:
        i_ans_end = len(folder)
    return mydtype(folder[i_ans_start:i_ans_end])

def make_Ramp(ramp_colors, show=True): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    if show:
        plt.figure( figsize = (15,3))
        plt.imshow( [list(np.arange(0, len( ramp_colors ) , 0.1)) ] , interpolation='nearest', origin='lower', cmap= color_ramp )
        plt.xticks([])
        plt.yticks([])
    return color_ramp

def make_Ramp_flexible(ramp_colors, show=True): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ c1 for c1 in ramp_colors ] )
    if show:
        plt.figure( figsize = (15,3))
        plt.imshow( [list(np.arange(0, len( ramp_colors ) , 0.1)) ] , interpolation='nearest', origin='lower', cmap= color_ramp )
        plt.xticks([])
        plt.yticks([])
    return color_ramp

def create_napari_colormap(mycolors, name='mycolormap', show=True):
    """
    Function that takes a list of colors (can be in names) and converts them
    to a napari colormap

    Parameters
    ----------
    mycolors : TYPE
        DESCRIPTION.
    name : TYPE, optional
        DESCRIPTION. The default is 'mycolormap'.
    show : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    napari_colormap : TYPE
        DESCRIPTION.

    """
    from vispy.color import Colormap
    from matplotlib.colors import to_rgba
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap
    napari_colormap = {'colors': [to_rgba(each) for each in mycolors],
                      'name':name,
                      'interpolation': 'linear'}
    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ c1 for c1 in mycolors ] )
    
    if show:
        plt.figure( figsize = (15,3))
        plt.imshow( [list(np.arange(0, len( mycolors ) , 0.1)) ] , interpolation='nearest', origin='lower', cmap= color_ramp )
        plt.xticks([])
        plt.yticks([])
    return napari_colormap

def cmap_napari(list_of_colors):
    #this function is designed to create vispy colormaps for napari that are labels, but made as images
    
    #Example
    #viewer = napari.Viewer() 
    #this_cmap = cmap_napari(['#000000', 'yellow', (1.0,1.0,1.0)])
    #viewer.add_image(cellclass_pyr, colormap=('this_cmap', this_cmap), contrast_limits=[0,4])
    from vispy.color import Colormap
    rgb_list_of_colors = []
    for each in list_of_colors:
        if type(each) == str:
            rgb_list_of_colors.append(colors.to_rgb(each))
        else:
            rgb_list_of_colors.append(each)
    
    my_colormap = Colormap(rgb_list_of_colors)
    display(my_colormap)
    return my_colormap

def get_dims_from_numpy_file(fullpath_in):
    myarray = np.load(fullpath_in, mmap_mode='r')
    return myarray.shape

def change_classes(array, list0, list1):
    outarray = deepcopy(array)
    for (each0, each1) in zip(list0,list1):
        outarray[array==each0] = each1
    return outarray



def get_pyramid_hybrid_loading(array, max_size_of_min_layer, max_actual_load_size):
    #max_actual_load_size is the threshold, below which every pyramid layer is loaded into memory
    #max_size_of_min_layer is the maximum size of the smallest layer of the pyramid
    #making the gap between these two parameters wider increases the number of layers to be loaded into memory
    #this means slower load, but smoother inital scrolling
    #currently only work if napari window area is not too large
    
    #max_size_of_min_layer <=8192
    #Only tested for the case where max_size_of_min_layer is >=2*max_size_of_min_layer
    
    mydimsiterator = array.shape
    pyramid_fn = [array]
    factor_it = 2

    #This part is for the non-loaded arrays (just returns views of original)
    while np.max(mydimsiterator) > max_actual_load_size:
#        print('Starting = ', mydimsiterator, factor_it, 'Lazy load')
        newarray = array[::factor_it, ::factor_it]

        pyramid_fn.append(newarray)
        mydimsiterator = newarray.shape
#        print(mydimsiterator, factor_it, 'Lazy load')
        factor_it = factor_it * 2

    #This part is for the arrays that are loaded into memory before napari runs (actual load)
    while np.max(mydimsiterator) > max_size_of_min_layer:
#        print(factor_it)
        newarray = np.copy(array[::factor_it, ::factor_it])
        pyramid_fn.append(newarray)

        mydimsiterator = newarray.shape
#        print(mydimsiterator, factor_it, 'Copy load')
        factor_it = factor_it * 2

    return pyramid_fn


def get_dask_pyramid_hybrid_loading(array, max_size_of_min_layer, max_actual_load_size):
    #max_actual_load_size is the threshold, below which every pyramid layer is loaded into memory
    #max_size_of_min_layer is the maximum size of the smallest layer of the pyramid
    #making the gap between these two parameters wider increases the number of layers to be loaded into memory
    #this means slower load, but smoother inital scrolling
    #currently only work if napari window area is not too large
    
    #max_size_of_min_layer <=8192
    #Only tested for the case where max_size_of_min_layer is >=2*max_size_of_min_layer
    
    mydimsiterator = array.shape
    pyramid_fn = [da.from_array(array)]
    factor_it = 2

    #This part is for the non-loaded arrays (just returns views of original)
    while np.max(mydimsiterator) > max_actual_load_size:
#        print('Starting = ', mydimsiterator, factor_it, 'Lazy load')
        newarray = array[::factor_it, ::factor_it]

        pyramid_fn.append(da.from_array(newarray))
        mydimsiterator = newarray.shape
#        print(mydimsiterator, factor_it, 'Lazy load')
        factor_it = factor_it * 2

    #This part is for the arrays that are loaded into memory before napari runs (actual load)
    while np.max(mydimsiterator) > max_size_of_min_layer:
#        print(factor_it)
        newarray = np.copy(array[::factor_it, ::factor_it])
        pyramid_fn.append(da.from_array(newarray))

        mydimsiterator = newarray.shape
#        print(mydimsiterator, factor_it, 'Copy load')
        factor_it = factor_it * 2

    return pyramid_fn


def make_dask_pyramid_from_dask_array(dask_array, smallest_loading_level, largest_loading_level,
                                      show=False):
    """
    dask_array is 2D, can be loading with dask_image.imread, can be .rechunk ed 
    smallest_loadest_level is the smallest level of the pyramid to create
    largest_loading_level is the largest level of the pyramid to load into memory
    All other levels will be left as dask arrays    
    
    """
    myshape = dask_array.shape
    output_pyr = []
    current_shape = myshape
    shrink_factor = 2
    current_shrink_factor = 1
    base_memory_array_exists = False
    base_memory_shrink_factor = 0
    
    while np.mean(current_shape) > smallest_loading_level:
        
        if show:
            print('Current shape = ', current_shape)
        
        if np.mean(current_shape) > largest_loading_level:
            output_pyr.append(dask_array[::current_shrink_factor, ::current_shrink_factor])
        
        
        elif np.mean(current_shape) > smallest_loading_level and not base_memory_array_exists:
            if show:
                print('Getting base memory array')
            base_memory_array = dask_array[::current_shrink_factor, ::current_shrink_factor].compute()
            output_pyr.append(base_memory_array)
            base_memory_array_exists = True
            base_memory_shrink_factor = 1
            
        elif np.mean(current_shape) > smallest_loading_level and base_memory_array_exists:
            this_array = base_memory_array[::base_memory_shrink_factor, ::base_memory_shrink_factor]
            output_pyr.append(deepcopy(this_array))
            
            if show:
                print('Subsampling base memory')
                print('Just appended this array ', this_array.shape)
        
        current_shrink_factor = current_shrink_factor*shrink_factor
        current_shape = (int(current_shape[0]/shrink_factor), int(current_shape[1]/shrink_factor))
        base_memory_shrink_factor = base_memory_shrink_factor*shrink_factor
        

        
    return output_pyr


def mycc(mask1, mask2):
    #boolean arrays
    CC = (2*np.sum(mask1 & mask2))/(np.sum(mask1)+np.sum(mask2))
    return CC

def get_cc_array(masks1, masks2):
    #masks must have the same number of masks
    uniques1 = np.sort(np.unique(masks1))
    uniques2 = np.sort(np.unique(masks2))
#    if not np.array_equal(uniques1, uniques2):
#        print(uniques1, uniques2)
#        print('These masks are different')
#        return -1
    if len(uniques1) != len(uniques2):
        print('Different number of masks')
        return -1
    
    N = len(uniques1)
    
    cc_array = np.zeros([N+1, N+1])

    for i in uniques1:
        for j in uniques2:
            thismask1 = (masks1==i)
            thismask2 = (masks2==j)
            CC = mycc(thismask1, thismask2)
            cc_array[i,j] = CC
    return cc_array

def get_cc_averaged_over_array(masks1, masks2):
    cc_array = get_cc_array(masks1, masks2)
    maxis = [np.max(x) for x in cc_array];
    return np.mean(maxis)


def scale_max_to_new_dtype(myarray, mydtype):
    dictionary_max = {np.uint8: 255, np.uint16: 65535}
    maxi = np.max(myarray)
    
    output = np.zeros(myarray.shape, dtype = mydtype)
    tmparray = myarray/maxi*dictionary_max[mydtype]
    output = tmparray.astype(mydtype)
    return output


def scale_max_to_new_dtype_low_memory(myarray, mydtype, set_max_threshold):
    """
    This only works for 4D array that is of a certain size so that a 2D slice can be processed in it
    
    Set max threshold is the number in the current units that gets
    converted into the maximum value of the new array. Everything above
    that is set to the same maximum.
    """
    dictionary_max = {np.uint8: 255, np.uint16: 65535}
    output = np.zeros(myarray.shape, dtype = mydtype)
    scale = dictionary_max[mydtype] / set_max_threshold
    print('Scale = ', scale)
    if len(myarray.shape) == 4:
        for i in range(myarray.shape[0]):
            print('i=',i)
            for j in tqdm(range(myarray.shape[1])):
                tmp_array = deepcopy(myarray[i, j, :, :])
                tmp_array[tmp_array > set_max_threshold] = set_max_threshold
                output[i, j, :, :] = myarray[i, j, :, :] * scale
                
    if len(myarray.shape) == 3:
        for i in tqdm(range(myarray.shape[0])):
            tmp_array = deepcopy(myarray[i, :, :])
            tmp_array[tmp_array > set_max_threshold] = set_max_threshold
            output[i, :, :] = myarray[i, :, :] * scale
            
    return output


def cmap_threshold(lower, mid, upper, color1, color2):
    #this function is designed to create vispy colormaps for napari
    #these colormaps are designed to show the effect of the threshold
    #color1 and color2 are designed to be RGB lists but can do conversion
    
    #lower, mid and upper and mid, must be ints
    #lower and upper must be the same as the contrast_limits set, otherwise the threshold won't be in the right place
    
    #Example
    #viewer = napari.Viewer()
    #my_cmap = uf.cmap_threshold(100, 1000, 2000, colors.to_rgb('blue'), colors.to_rgb('orange'))
    #viewer.add_image(myarray, colormap=('my_cmap', my_cmap), contrast_limits=[100,2000])
    from vispy.color import Colormap
    
    my_colormap_list = []
    
    m = 1/(upper-lower)
    c = (-lower) / (upper - lower)
    
    for i in range(lower, mid):
        lower_scale_constant = i * m + c
        my_colormap_list.append([lower_scale_constant*color1[0], lower_scale_constant*color1[1],lower_scale_constant*color1[2]])
    
    for i in range(mid, upper+1):
        upper_scale_constant = i * m + c
        my_colormap_list.append([upper_scale_constant*color2[0], upper_scale_constant*color2[1],upper_scale_constant*color2[2]])
    
    my_colormap = Colormap(my_colormap_list)
    display(my_colormap)
    return my_colormap






def get_euclidean_distance_python_wrapper_for_list(list_vec1, list_vec2):
    import numba_funcs as nf
#    print(type(list_vec1))
    #This is useless for FISHDBC because things are an array there
    np_vec1 = np.array(list_vec1)
    np_vec2 = np.array(list_vec2)
    return nf.get_euclidean_distance(np_vec1, np_vec2)


def bokeh_imshow(im, color_mapper=None, plot_height=400, length_units='pixels', 
                 interpixel_distance=1.0):
    """
    Display an image in a Bokeh figure.
    
    Parameters
    ----------
    im : 2-dimensional Numpy array
        Intensity image to be displayed.
    color_mapper : bokeh.models.LinearColorMapper instance, default None
        Mapping of intensity to color. Default is 256-level Viridis.
    plot_height : int
        Height of the plot in pixels. The width is scaled so that the 
        x and y distance between pixels is the same.
    length_units : str, default 'pixels'
        The units of length in the image.
    interpixel_distance : float, default 1.0
        Interpixel distance in units of `length_units`.
        
    Returns
    -------
    output : bokeh.plotting.figure instance
        Bokeh plot with image displayed.
    """
    import bokeh.io
    import bokeh.models
    import bokeh.palettes
    import bokeh.plotting
    import bokeh    
    
    # Get shape, dimensions
    n, m = im.shape
    dw = m * interpixel_distance
    dh = n * interpixel_distance
    
    # Set up figure with appropriate dimensions
    plot_width = int(m/n * plot_height)
    p = bokeh.plotting.figure(plot_height=plot_height, plot_width=plot_width, 
                              x_range=[0, dw], y_range=[0, dh], x_axis_label=length_units,
                              y_axis_label=length_units,
                              tools='pan,box_zoom,wheel_zoom,reset')

    # Set color mapper; we'll do Viridis with 256 levels by default
    if color_mapper is None:
        color_mapper = bokeh.models.LinearColorMapper(bokeh.palettes.viridis(256))

    # Display the image
    im_bokeh = p.image(image=[im[::-1,:]], x=0, y=0, dw=dw, dh=dh, 
                       color_mapper=color_mapper)
    
    return p

def create_df_cross_for_separate_dfs(df1, df2):
    """A function that does the outer product of two pandas dataframes
    when they have no common columns"""
    df1['tmp'] = 1
    df2['tmp'] = 1
    df = pd.merge(df1, df2, on=['tmp'])
    df = df.drop('tmp', axis=1)
    return df



def get_cropped_superpixel_dims(img_shape, superpx_size):
    """
    Takes dims as input and returns dims that have been cropped to superpixel size
    
    """
    x0_diff = img_shape[0]%superpx_size
    x1_diff = img_shape[1]%superpx_size
    
    out_shape = (img_shape[0]-x0_diff, img_shape[1]-x1_diff)
    return out_shape
    


def get_heatmap_as_rgb_numpy(data_points, pixel_size, xlim, ylim, vmax, cmap):
    """
    To be used with additive_white_background

    Parameters
    ----------
    data_points : Nx2 array of points to be put into 2D density
        DESCRIPTION.
    pixel_size : TYPE
        DESCRIPTION.
    xlim : TYPE
        DESCRIPTION.
    ylim : TYPE
        DESCRIPTION.
    vmax : TYPE
        DESCRIPTION.
    cmap : TYPE
        DESCRIPTION.

    Returns
    -------
    image_from_plot : RGB image of the plt.hist2d function
        To be used with additive_white_background

    """
    fig, ax = plt.subplots(1,1, dpi=1000)
    plt.hist2d(data_points[:,0], data_points[:,1], 
               (np.arange(xlim[0],xlim[1],pixel_size), 
                np.arange(ylim[0], ylim[1], pixel_size)),
               cmap=cmap, vmax=vmax)
    ax.axis('off')
    fig.tight_layout(pad=0)

    # To remove the huge white borders
    ax.margins(0)

    fig.canvas.draw()
    image_from_plot = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    image_from_plot = image_from_plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    plt.close()
    return image_from_plot




    
def additive_white_background(list_of_arrays, bg=255):
    """
    Function that takes a list of rgb numpy images (FROM)
    and 'adds' them together with white (255, 255, 255) as the 0
    Corrects for the case where the numbers go below zero
    

    Parameters
    ----------
    list_of_arrays : python list of numpy arrays in shape e.g. (640, 480, 3)
        DESCRIPTION.
    bg : number that define, optional
        DESCRIPTION. The default is 255.

    Returns
    -------
    RGB int image
        Can be used with plot
        ax[0].imshow(myarray0, extent = [xlim[0], xlim[1], ylim[0], ylim[1]])
        to make plt.plot or plt.scatter axes

    """
    list_of_float_arrays = [array.astype(float) for array in list_of_arrays]
    initial_array = bg*np.ones(list_of_float_arrays[0].shape, dtype=float)
    overall_array = initial_array
    
    for current_array in list_of_float_arrays:
        to_subtract = bg - current_array
        overall_array  = overall_array - to_subtract
    overall_array[overall_array<0] = 0
    return overall_array.astype(np.uint8)


def return_color_scale(key, show=True):
    """
    A place to save my colormaps that exist for different purposes
    
    Parameters
    ----------
    key : TYPE
        DESCRIPTION.
    show : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    out : TYPE
        DESCRIPTION.

    """
    if key == 'white_to_color_scale':
        out = ['purple', 'darkgreen', 'darkorange', 'blue', 'brown', 'darkcyan', 'yellow']
        if show:
            for each in out:
                make_Ramp(['white', each])
                
    if key=='white_to_faded_color_scale':
        """
        This scale is for plotting lines that fade to white
        Not colorblind compatible
        Intended for separate plots
        #yes: 'yellow', 
        'darkorange', 'springgreen', 'lawngreen', 
        'aqua', 'cyan', 'orange', 'black', 
        #no: ''green', 'brown', 'limegreen', 
        'forestgreen', 'seagreen', 'lightgreen', 
        'gold', 'darkkhaki', 'khaki', 'maroon',
        #'blueviolet'
        
        """
        out =  ['blue', 'orange','fuchsia', 'lawngreen', 'cyan', 'red',
                'dimgray', 'springgreen']
        if show:
            for each in out:
                make_Ramp(['white', each])
                
    if key=='block_colors_for_labels_against_gray':
        """
        This scale is for plotting labels onto a 
        grayscale image so that the grayscale can be seen 
        as separate, and the colors can be used with a uniform
        value of alpha that can vary between 0 and 1 simulataneously 
        for all colours
        #removed 'peachpuff', 'navajowhite', 'lightsteelblue', 
        'steelblue', 'sandybrown','darkolivegreen','lawngreen',
        'purple', (0.3,0,0.3)
        
        """
        out =  ['dodgerblue', 'gold', 'midnightblue', 
                'lightskyblue','red', 'brown', 'purple', 'palegreen',
                'pink'] #yellow?
        #not sure about midnightblue in some cases
        if show:
            for each in out:
                make_Ramp([each, each])
                
                
    if key=='block_colors_for_labels_against_gray2':
        """
        This scale is for plotting labels onto a 
        grayscale image so that the grayscale can be seen 
        as separate, and the colors can be used with a uniform
        value of alpha that can vary between 0 and 1 simulataneously 
        for all colours
        #removed 'peachpuff', 'navajowhite', 'lightsteelblue', 
        'steelblue', 'sandybrown','darkolivegreen','lawngreen',
        'purple', (0.3,0,0.3)
        
        """
        out =  ['purple', 'darkgreen', 'darkorange', 'blue', 'red', 'yellow'] #yellow?
        #not sure about midnightblue in some cases
        if show:
            for each in out:
                make_Ramp([each, each])
                
    if key=='block_colors_for_labels_general':
        """
        This scale is for plotting labels where they don't mix
        
        """
        out =  ['purple', 'darkgreen', 'darkorange', 'blue', 'deeppink', 'yellow',
                'black', 'gray', 'midnightblue', 'white', 'lightskyblue', 'pink',
                'peachpuff']
        
        if show:
            for i, each in enumerate(out):
                make_Ramp([each, each])
                plt.text(-1.2,0,str(i), fontsize=20, va='center', ha='left')
                plt.show()
                
    if key == 'block_colors_for_labels_against_white':
        """
        This is for labels that don't mix and takes into account to some extent
        Omer's iniital preferences
        """
        out =  ['black', 'blue', 'darkorange', 'deeppink', 'darkgreen', 
                'gray', 'lightskyblue', 'purple', 
                 'gold', 'silver', 'darkgoldenrod',  'deepskyblue', ]
        #close? 'midnightblue', 'royalblue', 'red', 'steelblue',
        #removed#'brown''saddlebrown''sienna', 'peru', 'chocolate',
        #'darkturquoise', 'cadetblue', 'navy', 'darkblue', 
        # 'dodgerblue','darkcyan'
        if show:
            for each in out:
                make_Ramp([each, each])
#                plt.text(-1.2,0,each, fontsize=20, va='center', ha='left')
                plt.show()
                
    if key == 'block_colors_for_labels_against_white_small_points':
        """
        This is for labels that don't mix and but look gray as the points are 
        small. Used for IVY GAP project as colorblind_optimized
        """
        yy = [(1.0, 0.6, 0.6, 1.0),
             (0.0, 0.39215686274509803, 0.0, 1.0),
             (1.0, 0.0784313725490196, 0.5764705882352941, 1.0),
             (0.65, 0.45, 0.04, 1.0), #(0.7215686274509804, 0.5254901960784314, 0.043137254901960784, 1.0)
             (1.0, 0.5490196078431373, 0.0, 1.0),
             (0.0, 0.7490196078431373, 1.0, 1.0),
             (0.0, 0.0, 1.0, 1.0),
             (1.0, 0.8431372549019608, 0.0, 1.0),
             (0.5019607843137255, 0.0, 0.5019607843137255, 1.0),
             (0.0, 0.0, 0.33, 1.0),
                (0.0, 0.0, 0.0, 1.0),
                 (0.33, 0.33, 0.33, 1.0),
                 (0.66, 0.66, 0.66, 1.0),
                 (1.0, 1.0, 1.0, 1.0)]
        yy_hex = [to_hex(each) for each in yy]
        out = yy_hex
        if show:
            for i, each in enumerate(out):
            #    print(i)
                make_Ramp_flexible([each, each])
                plt.show()

    if key == 'block_colors_for_labels_against_white_small_points_white_first':
        """
        This is for labels that don't mix and but look gray as the points are 
        small. Used for IVY GAP project as colorblind_optimized / 'white_first'
        """
        yy = [(1.0, 0.6, 0.6, 1.0),
             (0.0, 0.39215686274509803, 0.0, 1.0),
             (1.0, 0.0784313725490196, 0.5764705882352941, 1.0),
             (0.65, 0.45, 0.04, 1.0), #(0.7215686274509804, 0.5254901960784314, 0.043137254901960784, 1.0)
             (1.0, 0.5490196078431373, 0.0, 1.0),
             (0.0, 0.7490196078431373, 1.0, 1.0),
             (0.0, 0.0, 1.0, 1.0),
             (1.0, 0.8431372549019608, 0.0, 1.0),
             (0.5019607843137255, 0.0, 0.5019607843137255, 1.0),
             (0.0, 0.0, 0.33, 1.0),
                (0.0, 0.0, 0.0, 1.0),
                 (0.33, 0.33, 0.33, 1.0),
                 (0.66, 0.66, 0.66, 1.0),
                 (1.0, 1.0, 1.0, 1.0)]
        yy_hex = [to_hex(each) for each in yy]
        out = yy_hex
        if show:
            for i, each in enumerate(out):
            #    print(i)
                make_Ramp_flexible([each, each])
                plt.show()
                
    if key == 'two_colors_for_overlap_plotting':
        """This is for situations where you want two overlapping colors in 
        matplotlib plot. Overlap will be purple to black
        Purple seems to look black in colorblindness
        And both colours can be distinguished These 
        colors also fade out to faded versions of themselves
        with alpha, unlike some other colours"""
        out = ['blue', 'darkorange']
        if show:
            for each in out:
                make_Ramp(['white', each])
                
    if key == 'black_to_color_scale':
        out = ['blue', 'yellow']
        if show:
            for each in out:
                make_Ramp(['black', each])
    
    if key == 'Omer_chosen_colors':
        out = ['yellow', 'lightskyblue', 'peachpuff', 'darkorange', 
               'blue', 'deeppink']
        if show:
            for each in out:
#                print('\n', each)
#                print(colors.to_rgb(each))
                make_Ramp([each, each])
#                plt.show()

    if key == 'Omer_chosen_colors_dark_first':
        out = ['blue', 'darkorange', 'deeppink', 'yellow', 'lightskyblue', 
               'peachpuff']
        if show:
            for each in out:
#                print('\n', each)
#                print(colors.to_rgb(each))
                make_Ramp([each, each])
#                plt.show()

    if key == 'gray_first_Omer_chosen_colors_plus_my_color_blind_choices':
        out = ['lightgray', 'blue', 'darkorange', 'deeppink', 'yellow', 'lightskyblue', 
               'peachpuff', 'purple', 'darkgreen',
                'black', 'dimgray']
        if show:
            for each in out:
#                print('\n', each)
#                print(colors.to_rgb(each))
                make_Ramp([each, each])
#                plt.show()

    if key == 'gray_first_Omer_chosen_colors_plus_my_color_blind_choices_plus_mcolorsx1000':
        """This part is for the agglomerative clustering that has unlimited colors
        because it is clustering through texture space so we are not protected by the 
        theorem that says we only need 4 or 5 or 6 colors
        """
        import matplotlib.colors as mcolors
        colors3 = list(mcolors.CSS4_COLORS.keys())
        colors3.remove('rebeccapurple')
        colors3.remove('teal')
        colors3.remove('white')
        colors3.remove('whitesmoke')
        colors3.remove('ghostwhite')
        colors3.remove('black')
        colors3.remove('darkgrey')
        colors3.remove('dimgrey')
        colors3.remove('lightgray')
        colors3.remove('lightgrey')
        colors3.remove('grey')
        colors3.remove('lightslategrey')
        colors3.remove('slategrey')
        colors3.remove('darkslategrey')
        colors_to_repeat = deepcopy(colors3)
        out = ['lightgray', 'blue', 'darkorange', 'deeppink', 'yellow', 'lightskyblue', 
                   'peachpuff', 'purple', 'darkgreen',
                    'black', 'dimgray'] + colors_to_repeat*1000
        if show:
            for each in out[0:200]:
#                print('\n', each)
#                print(colors.to_rgb(each))
                make_Ramp([each, each])
                plt.show()
                plt.close()
                
        

    return out

def create_mask_from_list_of_coords(dims, list_of_coords, scale_up = 1):
    import numba_funcs as nf
    """
    This function takes a list of coords [[0,0], [1,1]]
    and makes a mask image out of them 

    Parameters
    ----------
    dims : tuple of ints
        DESCRIPTION.
    list_of_coords : TYPE
        DESCRIPTION.
    scale_up : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    outarray : numpy array
        Full size image from mask coordinates

    """
    myarray = np.zeros(dims, dtype=np.bool_)
    for each in list_of_coords:
        myarray[each[0],each[1]] = 1
    outarray = nf.upsize(myarray, scale_up)
    return outarray

def create_heatmap_from_list_of_coords(dims, list_of_X0, list_of_X1, 
                                       list_of_class, scale_up = 1,
                                       bg_value=0, dtype=np.uint8):
    import numba_funcs as nf
    """
    Function that adds class integers onto a heatmap

    Parameters
    ----------
    list_of_X0 : TYPE
        DESCRIPTION.
    list_of_X1 : TYPE
        DESCRIPTION.
    list_of_class : TYPE
        DESCRIPTION.
    scale_up : TYPE, optional
        DESCRIPTION. The default is 1.
    bg_value :
        The value of the background that is needed.
        Usually either 0 or 1

    Returns
    -------
    outarray : TYPE
        DESCRIPTION.

    """
#    myarray = np.zeros((np.max(list_of_X0)+1, np.max(list_of_X1)+1), dtype=np.uint8)
    myarray = bg_value * np.ones(dims, dtype=dtype)
    for X0, X1, myclass in zip(list_of_X0, list_of_X1, list_of_class):
        myarray[X0,X1] = myclass
    outarray = nf.upsize(myarray, scale_up)
    return outarray


    
    
def rgb2labelint(img, array_of_colors = None):
    import numba_funcs as nf
    if array_of_colors is None:
        print('Started getting colors')
        array_of_colors = np.unique(img.reshape(-1, img.shape[2]), axis=0)
        print('Finished getting colors')
    output = nf.rgb2labelint_iterator(img, array_of_colors)
    return output



def label2rgb_with_dict(myarray, dict_of_colors):
    """
    This function takes an image (or other numpy array) and maps the 
    original numbers to the colors specified in the dictionary.
    So a '4' in the original array will get mapped to dict_of_colors[4],
    even if dict_of_colors[3], or dict_of_colors[2] (etc.) do not exist
    in the image. This is useful to maintain consistent colors across 
    images where each image is missing some of the classes.

    Parameters
    ----------
    myarray : numpy array 
        Normally integers that represent classes across an image
    dict_of_colors : dictionary
        The keys of this dictionary should correspond to the values in 
        myarray. The output values of the dictionary should be colors,
        and I think can be matplotlib names, hexes or rgb values

    Returns
    -------
    x : numpy array, image output
        float64 of image output, with values between 0 and 1

    """
    from skimage.color import label2rgb
    from matplotlib.colors import to_rgb
    uniques = np.unique(myarray)
    list_of_colors = [to_rgb(dict_of_colors[each]) for each in uniques]
    bg = np.zeros(myarray.shape, np.float64)
    x = label2rgb(myarray, bg, list_of_colors, alpha=1, bg_label=-100000.5)
    return x


def force_image_into_common_colors(image, N_pixels):
    """
    This function finds the most common colors in the image
    and then forces the others into them. Can be used in a mode 
    that selects all colors above a certain number of pixels. 
    

    Parameters
    ----------
    image : TYPE
        DESCRIPTION.
    N_pixels : TYPE
        DESCRIPTION.

    Returns
    -------
    output : TYPE
        DESCRIPTION.

    """
    output = np.zeros(image.shape, dtype=image.dtype)
    array_of_all_colors, counts = np.unique(image.reshape(-1, image.shape[2]), axis=0, return_counts=True)
    array_of_common_colors = array_of_all_colors[counts > N_pixels]
    print('Number of common colors =', len(array_of_common_colors))
    output = rgb2commonrgbclosest_iterator(image, array_of_common_colors)
    return output



def convert_rgb_to_rgba_without_transparency(img):
    """
    Converts rgb to rgba without transparency
    There is a numba func version of this that does optional transparency
    Parameters
    ----------
    img : numpy array (of floats i think, rgb)
        DESCRIPTION.

    Returns
    -------
    rgba_image : numpy array like img 
        RGBA version of img

    """
    shape = img.shape
    rgba_image = np.ones((shape[0], shape[1], 4), img.dtype)
    rgba_image[:,:, 0:3] = img
    
    return rgba_image


def get_optimal_overlap_of_classes(observed, gt, 
                                   elements_to_remove_from_observed=[-1],
                                   elements_to_remove_from_gt = []):
    import numba_funcs as nf
    import fastremap
    """
    This function takes two pandas series / lists/ numpy arrays
    And figures out which of the different classes in each best pair
    with the other image to produce maximum overlap.
    
    The observed and gt objects do not have to have the same number of classes.
    
    -1 in the observed object is treated as a special case
    

    Parameters
    ----------
    observed : ANYTHING
        DESCRIPTION.
    gt : ANYTHING
        DESCRIPTION.

    Returns
    -------
    out_dict
    overlap score

    """
    #This helps get the data in a consistent format
    observed_np = np.array(observed).ravel()
    gt_np = np.array(gt).ravel()

    #This removes both observed and gt locations where observed is in elements_to_remove_from_observed
    observed_to_remove = np.isin(observed_np, elements_to_remove_from_observed)
    if observed_to_remove.any():
        observed_np = observed_np[~observed_to_remove]
        gt_np = gt_np[~observed_to_remove]

    #This removes both observed and gt locations where gt is in elements_to_remove_from_gt
    gt_to_remove = np.isin(gt_np, elements_to_remove_from_gt)
    if gt_to_remove.any():
        observed_np = observed_np[~gt_to_remove]
        gt_np = gt_np[~gt_to_remove]
    
    
    #This removes both observed and gt locations where observed is -1
    #Deprecated for section above
#    observed_np_minus_1 = observed_np == -1
#    if observed_np_minus_1.any():
#        observed_np = observed_np[~observed_np_minus_1]
#        gt_np = gt_np[~observed_np_minus_1]

    #This section gets the unique values, then assigns them to ints from 0, then remaps the original 
    OBS_cats = pd.unique(observed_np)
    #display(OBS_cats)

    GT_cats = pd.unique(gt_np)
    #display(GT_cats)

    OBS_int_list = np.array([i for i in range(len(OBS_cats))])
    #print(OBS_int_list)

    OBS_mapped_to_ints = np.vectorize(dict(zip(OBS_cats, OBS_int_list)).get)(observed_np)
    #display(OBS_mapped_to_ints)

    GT_int_list = np.array([i for i in range(len(GT_cats))])
    #print(GT_int_list)

    GT_mapped_to_ints = np.vectorize(dict(zip(GT_cats, GT_int_list)).get)(gt_np)
    #print(GT_mapped_to_ints)

    #This part gets the iterator depending on the different 
    if len(OBS_int_list) == len(GT_int_list):
        #then can iterate over OBS_int_list
        myiter = itertools.permutations(OBS_int_list, len(OBS_int_list))
    elif len(OBS_int_list) > len(GT_int_list):
        #then can iterate over OBS_int_list
        myiter = itertools.permutations(OBS_int_list, len(OBS_int_list))
    elif len(OBS_int_list) < len(GT_int_list):
        #then can iterate over subsampling of GT_int_list
        myiter = itertools.permutations(GT_int_list, len(OBS_int_list))

    #This part iterates over all possible options and remembers the best one
    overlap_score = 0.0
    for this_perm in myiter:
        temp = fastremap.remap(OBS_mapped_to_ints, dict(zip(OBS_int_list,this_perm)), preserve_missing_labels=False)
        score = nf.IoU_labels_images_1d(temp, GT_mapped_to_ints)
        if score > overlap_score:
            overlap_score = score
            best_perm = this_perm
#    print(overlap_score, best_perm)

    #This part converts (by inversion) the matched integers to the original pixel values or objects
    OBS_inv_dict = dict(zip(OBS_int_list, OBS_cats))
    GT_inv_dict = dict(zip(GT_int_list, GT_cats))

    out_dict = {}
    for a,b in zip(OBS_int_list, best_perm):
        if b in GT_inv_dict.keys():
            out_dict[OBS_inv_dict[a]] = GT_inv_dict[b]
        else:
            out_dict[OBS_inv_dict[a]] = np.nan

    return out_dict, overlap_score


def intlabels_and_rgblabels_to_dict(intlabels, rgblabels, list_to_ignore = [0], divide_by_255 = False):
    """
    A function that takes 2D intlabels and rgblabels and makes the link.
    Not optimised or compiled - iterating over a numpy array, so not 
    for use with large arrays as it will be slow. 

    Parameters
    ----------
    intlabels : numpy ints
        2D array of ints (can be XY or YX) that indicate classes on an image
    rgblabels : numpy anything
        XYC or YXC array of numpy that looks like a rgb image or channel image
        that matches the intlabels in XY coordinates, just looks differently
    list_to_ignore : list, optional
        This is the list of things to not add to the dictionary. 
        The default is [0].
    divide_by_255 : Bool, optional
        A decision of whether to by 255 for the output colors. 
        Napari needs colors in range 0-1 so if rgblabels is in 0-255 then
        this should be set to True
        The default is False.

    Returns
    -------
    int_to_rgb_dict : dict
        A dictionary that links the ints provided by intlabels to the colors
        provided by rgblabels of the same image
        

    """
    int_to_rgb_dict = {}
    if divide_by_255:
        for i in range(intlabels.shape[0]):
            for j in range(intlabels.shape[1]):
                if intlabels[i,j] in int_to_rgb_dict.keys():
                    continue
                else:
                    int_to_rgb_dict[intlabels[i,j]] = rgblabels[i,j]/255
    else:
        for i in range(intlabels.shape[0]):
            for j in range(intlabels.shape[1]):
                if intlabels[i,j] in int_to_rgb_dict.keys():
                    continue
                else:
                    int_to_rgb_dict[intlabels[i,j]] = rgblabels[i,j]
    for each in list_to_ignore:
        del int_to_rgb_dict[each]
    return int_to_rgb_dict