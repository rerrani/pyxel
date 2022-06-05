# pyxel.py
This bit of code can be used to generate surface density maps of sparse simulation data. It takes 2D coordinates of N-body particles, maps them on a grid, and subsequently expands each non-empty grid cell until it contains a set minimum number of particles. Each pixel of the output grid is assigned the average surface density of all expanded grid cells that overlap on it. The figure below illustrates how the non-empty grid cells are expanded until each of them contains at least two particles. The resolution in different regions of the output grid can be intuitively read off from the expanded grid cell size.

![illustration](https://github.com/rerrani/rerrani.github.io/blob/master/code/pyxel.png?raw=true)

The code runs a minimum working example upon launch: The ASCII file Plummer.dat is read, which contains x,y coordinates of an N-body Plummer model (shown below in the left-hand panel). The output is written to the ASCII file map.dat, which contains x,y pixel coordinates, as well as the average surface brightness at each pixel (plotted below in the right-hand panel).

![application](https://github.com/rerrani/rerrani.github.io/blob/master/code/pyxel_plummer.png?raw=true)
