from environment import *
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

class reactormesh:
    """Class to arange reactors in a mesh for spatially resolving them for
    e.g. oceans, biofilms.

    """

    xvals = []
    yvals = []
    zvals = []

    def __init__(self, sphere=False, spheredim=[], cuboid=False, cuboiddim=[]):
        """For a sphere, pass dimensions [[inner,outer],[r, theta, phi]]
        where inner and outer are the inner and outer radii, and r, theta, phi
        the number of cells you want in each sperical direction.


        For a cuboid, pass dimensions [[x, y, z],[nx, ny, nz]] where n
        means send the nuber of cells you want for that parameter.

        """
        if sphere==True:

            # use spherical coordinates
            theta = np.linspace(0, np.pi, spheredim[1][1], endpoint=False)
            phi = np.linspace(0, 2*np.pi, spheredim[1][2], endpoint=False)
            rvals = np.linspace(spheredim[0][0], spheredim[0][1], spheredim[1][0])
            r, theta, phi = np.meshgrid(rvals, theta, phi)
            self.xvals = r * np.sin(theta) * np.cos(phi)
            self.yvals = r * np.sin(theta) * np.sin(phi)
            self.zvals = r * np.cos(theta)



        if cuboid==True:


            #length of mesh in each dimension, and number of cells
            lx, ly, lz = (cuboiddim[0][0], cuboiddim[0][1], cuboiddim[0][2])
            nx, ny, nz = (cuboiddim[1][0], cuboiddim[1][1], cuboiddim[1][2])
            # half-length of each cell in each dimension
            halfx, halfy, halfz = (lx/(2*nx), ly/(2*ny), lz/(2*nz))

            # make a meshgrid of the coordinates, with the orgin as the center.
            # x = np.linspace(-cuboiddim[0][0]/2, cuboiddim[0][0]/2., nx)
            # y = np.linspace(-cuboiddim[0][1]/2., cuboiddim[0][1]/2., ny)
            # z = np.linspace(-cuboiddim[0][2]/2., cuboiddim[0][2]/2., nz)

            # make a meshgrid of the center of each cell, centered on the origin.
            x = np.linspace(-halfx*(nx-1), halfx*(nx-1), nx)
            y = np.linspace(-halfy*(ny-1), halfy*(ny-1), ny)
            z = np.linspace(-halfz*(nz-1), halfz*(nz-1), nz)
            self.xvals, self.yvals, self.zvals = np.meshgrid(x, y, z)


    def plotcentres(self, ax=None, show=True):
        if ax==None:
            # initialise the figure first
            fig = plt.figure()
            ax = plt.axes(projection='3d')
        ax.scatter3D(self.xvals, self.yvals, self.zvals, marker='x')
        if show:
            plt.show()
        return ax
