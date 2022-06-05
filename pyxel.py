#!/usr/bin/env python3
############################################################################
############################################################################
# This code generates density maps of sparse simulation data. It takes 2D  #
# coordinates of N-body particles as input, maps them on a grid, and then  #
# expands each non-empty grid cell until it contains a set minimum number  #
# of particles. See examples at https://rerrani.github.io/code.html#pyxel. # 
# RE 2022                                                                  #
############################################################################
############################################################################
import numpy as np
import multiprocessing as mp
from multiprocessing import shared_memory

# input ASCII file with X Y particle coordinates; one set of coordinates per row
input_filename  = "Plummer.dat"

# output ASCII file with X Y grid pixel coordinates and surface density at each pixel
output_filename = "map.dat"

Ncores    = 4          # Number of processor cores
bins      = 500        # output grid pixel per side
Xmin      = -10        # output min and max X, Y
Xmax      = +10
Ymin      = -10
Ymax      = +10
threshold = 5          # min number of particles per expanded cell

############################################################################
# reading input file
print ("* Reading input file '%s'"%input_filename)
data = np.loadtxt(input_filename)
xdata = data[:,0]
ydata = data[:,1]
input_histogram, xedges, yedges = np.histogram2d(xdata, ydata, bins=bins, range=[[Xmin, Xmax], [Ymin, Ymax]])
input_histogram = input_histogram.T
xdim, ydim = input_histogram.shape
xcentres = (xedges[1:] + xedges[:-1])/2.
ycentres = (yedges[1:] + yedges[:-1])/2.
x2D, y2D = np.meshgrid(xcentres, ycentres)

############################################################################
# setting up shared memory array. 
# shared_array_1_name and shared_array_2_name must be unique on a given computer
shared_array_1_name = "shared_array_1"
shared_array_2_name = "shared_array_2"
dummy = np.zeros((xdim,ydim))
shm1 = shared_memory.SharedMemory(create=True, size=dummy.nbytes, name=shared_array_1_name)
shm2 = shared_memory.SharedMemory(create=True, size=dummy.nbytes, name=shared_array_2_name)
merged_pixels       = np.ndarray( dummy.shape, dtype=dummy.dtype, buffer=shm1.buf  )
merged_pixels_count = np.ndarray( dummy.shape, dtype=dummy.dtype, buffer=shm2.buf  )
merged_pixels[:] = dummy[:]
merged_pixels_count[:] = dummy[:]
lock = mp.Lock()

############################################################################
# routine that sums the input histogram along cell edges
def makesum(x0, y0, dx, dy):
  avg = 0
  avg     += np.sum(        input_histogram[x0:x0+dx, y0] )
  avg     += np.sum(        input_histogram[x0:x0+dx, y0+dy] )
  avg     += np.sum(        input_histogram[x0,       y0+1:y0+dy-1])
  avg     += np.sum(        input_histogram[x0+dx,    y0+1:y0+dy-1])
  return avg

############################################################################
# routine that expands the grid cell centred in xi, yi
# until the expanded cell contains a min number of particles
def expand_cell(xi, yi):
  if (xi == xdim-1) or (yi == ydim -1): return 0  

  avg = 0.
  avg += input_histogram[xi, yi]
  for d in range(1,xdim):
    x0 = max(0, xi - d)
    y0 = max(0, yi - d)
    dx = min(2*d+1, xdim - x0- 1)
    dy = min(2*d+1, ydim - y0- 1)
    
    avg += makesum(x0, y0, dx, dy)
    if avg > threshold: break    
  avg /= (dx*dy) 

  shm1 = shared_memory.SharedMemory(name=shared_array_1_name)
  shm2 = shared_memory.SharedMemory(name=shared_array_2_name)
  lock.acquire()  
  merged_pixels       = np.ndarray( dummy.shape, dtype=dummy.dtype, buffer=shm1.buf)
  merged_pixels_count = np.ndarray( dummy.shape, dtype=dummy.dtype, buffer=shm2.buf)
  merged_pixels[x0:x0+dx, y0:y0+dy] += avg
  merged_pixels_count[x0:x0+dx, y0:y0+dy] += 1
  lock.release()
  shm1.close()
  shm2.close()
  
  counts = dx*dy 
  return counts

############################################################################
# multiprocessing: calling expand_cell on all non-zero grid cells
print ("* Merging pixels on %i cores"%Ncores)

xarray, yarray =  np.nonzero(input_histogram)
pool = mp.Pool(processes=Ncores)    
input_array = zip( xarray.flatten() , yarray.flatten() )
counts = pool.starmap(expand_cell, input_array )
pool.close()

counts                     = np.array(counts)
merged_pixels_count        = np.array(merged_pixels_count)
merged_pixels              = np.array(merged_pixels)
nonzeroIdx                 = np.nonzero(merged_pixels_count)
merged_pixels[nonzeroIdx] /= merged_pixels_count[nonzeroIdx]

shm1.close()
shm1.unlink()
shm2.close()
shm2.unlink()
print ("* Merging pixels done.")

############################################################################
# Output: ASCII file, format X Y Sigma, where X Y are pixel coordinates and 
#         Sigma is the corresponding surface brightness at that pixel
print ("* Writing output file '%s'"%output_filename)
output = open(output_filename, "w")
for xi in range(xdim):
  for yi in range(ydim):
    print (x2D[xi,yi], y2D[xi,yi], merged_pixels[xi,yi], file=output)
  print ( "", file=output)
    

