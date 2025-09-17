import numpy as np
import matplotlib.pyplot as plt
import copy

from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.stats import norm
from sklearn.cluster import KMeans
from skimage import measure
from statsmodels.distributions.empirical_distribution import ECDF
from typing import Tuple

class CopulaSegmentation():
    def __init__(self, data: np.ndarray):
        """
        Parameters:
          data: N-d array for the datacube of shape (filters, y_pixel, x_pixel)
        """
        self.data: np.ndarray = data
        self.num_filters, self.height, self.width = self.data.shape
        self.N_pixels = self.height * self.width

        self.X: np.ndarray = self.__restructure_data__()
        self.U, self.F = self.extract_empirical_copula()

    def __restructure_data__(self) -> np.ndarray:
        """
        Restructure the data into an N-D array with dimensions (num_pixels, 2+num_channels),
        where the first two columns are the x and y coordinates of each pixel

        Returns:
          X: N-D array of shape (num_pixels, 2+num_channels)
        """
        # Create x and y coordinates for each pixel
        x_coords, y_coords = np.meshgrid(np.arange(self.width), np.arange(self.height))
        
        # Flatten x and y coordinates
        x_coords_flat = x_coords.flatten()
        y_coords_flat = y_coords.flatten()
        
        # Flatten the channels of image data
        z_coords_flat = self.data.reshape(self.num_filters, -1).T
        
        # Combine x, y, and channels into a single dataset
        X: np.ndarray = np.column_stack((x_coords_flat, y_coords_flat, z_coords_flat))
        return X

    def extract_empirical_copula(self):
        """
        Obtain the empirical marginal distribution function and copula from the data
        """
        F = []
        U = []
        for dim in range(self.X.shape[1]):
            F.append(ECDF(self.X[:, dim]))
            # Transform each dimension marginally to a uniform using the ECDF
            U.append(F[-1](self.X[:, dim]))
        return np.array(U), np.array(F)

    def copula_KDE(self, U_data=None):
        if U_data is not None:
            U = np.expand_dims(np.array(copy.copy(U_data)).T, axis=2)
        else:
            U = np.expand_dims(np.array(copy.copy(self.U)).T, axis=2)
        
        # 3-d grid on which to evaluate the KDE
        u1, u2, u3 = np.meshgrid(np.linspace(0,1,20), np.linspace(0,1,20), np.linspace(0,1,20))
        u_grid = np.expand_dims(np.stack((u1,u2,u3)), axis=3)
        
        # number of data points
        N = U.shape[0]
        nrm_fac = 1.0/N
        
        # estimate bandwidth -- this may need tuning
        h = np.std(U, axis=0).squeeze() * N**(-1/5)
        
        # stores KDE
        kde_grid = np.zeros(u1.shape)
        
        # one-dimensional Gaussian kernel
        f = lambda u, U, h :  norm.pdf((u-U)/h)/h
    
        # loop over two dimesions, parallelize over the remaining
        for i in range(u1.shape[0]):
            print(i)
            for j in range(u1.shape[1]):
                u = u_grid[:,i,j,:,:]
                kde_grid[i,j,:] = nrm_fac * np.sum(f(u[0,:,:], U[:,0,:], h[0]) \
                                                   * f(u[1,:,:], U[:,1,:], h[1]) \
                                                   * f(u[2,:,:], U[:,2,:], h[2]),
                                                   axis=0)
        
        kde_data = nrm_fac * np.sum(f(U[:,0,:].T, U[:,0,:], h[0]) \
                                    * f(U[:,1,:].T, U[:,1,:], h[1]) \
                                    * f(U[:,2,:].T, U[:,2,:], h[2]),
                                    axis=0)
    
        return u_grid.squeeze(axis=3), kde_grid, kde_data

    @classmethod
    def gradient_field(cls, data: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute the gradient field of a 2D image array.
        Parameters:
          data: 2D array
        Returns:
          grad_y, grad_x: gradients along y and x axes
        """
        ny, nx = data.shape
        y = np.arange(ny)
        x = np.arange(nx)
        x_mesh, y_mesh = np.meshgrid(x, y)  # coordinate grid
        
        grad_y, grad_x = np.gradient(data, y, x)  # gradients wrt pixel indices
        grad_mag = np.sqrt(grad_x**2 + grad_y**2)
        return grad_y, grad_x, grad_mag
    
    def segmentation_map(self, type='k-means', num_clusters=20):
        """
        Generate a segmentation map
        """
        # U rows: [x, y, flux_1, ..., flux_N]
        fluxes = self.U[2:].T
        
        if type == 'k-means':
            kmeans = KMeans(n_clusters=num_clusters, random_state=0)
            labels = kmeans.fit_predict(fluxes)
            label_image = labels.reshape(self.height, self.width)
            
            contours_by_cluster = []
            for cluster_id in range(num_clusters):
                mask = (label_image == cluster_id).astype(float)
                contours = measure.find_contours(mask, 0.5)
                contours_by_cluster.append(contours)
        
        return label_image, contours_by_cluster
    

hdu = fits.open('datacube_reg1.fits')
data = hdu[0].data
cs = CopulaSegmentation(data)
label_image, contours_by_cluster = cs.segmentation_map(type='k-means', num_clusters=20)

fig, ax = plt.subplots(1, 1, figsize=(7,7))
ax.imshow(label_image, cmap='tab10', origin='lower')
for cluster_id, contours in enumerate(contours_by_cluster):
    for c in contours:
        x, y = c[:,1], c[:,0]
        ax.plot(x, y, '-k', linewidth=1)
ax.set_title("KMeans Clustering + Contours in Copula Space")
ax.axis("off")
plt.savefig('segmentation_map.png', dpi=300)

