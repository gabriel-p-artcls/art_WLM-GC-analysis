

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from astropy.io import ascii
from scipy.spatial import cKDTree
from scipy.stats import gaussian_kde


def main(file, method='KDE'):
    """
    method: KDE, kNN, kNNInvSum
    """
    print(file)

    data = ascii.read(file + ".txt")
    xs, ys = data['RA(deg)'], data['Dec(deg)']
    dim = 2

    resolution = 200

    # # extent = [np.min(xs), np.max(xs), np.min(ys), np.max(ys)]
    # xmed, xstd = np.median(xs), np.std(xs)
    # ymed, ystd = np.median(ys), np.std(ys)
    # extent = [
    #     xmed - 3. * xstd, xmed + 3. * xstd, ymed - 3. * ystd, ymed + 3. * ystd]

    # Using fixed limits so all plots have the same size.
    extent = [0.39, 0.61, -15.62, -15.31]

    kde_pars, kNN_pars = [], []
    if method == 'KDE':
        # Parameters used by the KDE dens maps method.
        kde_vals, x_grid, y_grid, kde_pos, scotts_f = KDEPrep(
            xs, ys, resolution, extent)
        vals = [0, .1, .2, .5]
        kde_pars = (kde_vals, x_grid, y_grid, kde_pos, scotts_f)
    else:
        # Data in plot coordinates, given a resolution.
        xv = data_coord2view_coord(xs, resolution, extent[0], extent[1])
        yv = data_coord2view_coord(ys, resolution, extent[2], extent[3])
        # Grid with proper shape.
        grid = np.mgrid[0:resolution, 0:resolution].T.reshape(
            resolution**2, dim)
        vals = [0, 10, 20, 40]
        kNN_pars = (xv, yv, grid)

    makePlot(
        file, xs, ys, resolution, extent, method, vals, kde_pars, kNN_pars)


def kNNInvSumDist(xv, yv, grid, resolution, neighbors):
    """
    Assigns the inverse of the sum of the k distances to each grid point,
    where k is the number of neighbors.

    Original idea: https://stackoverflow.com/a/36515364/1391441
    """
    # Create the tree
    tree = cKDTree(np.array([xv, yv]).T)
    # Distance to the nearest neighbors for all the grid points.
    dists = tree.query(grid, neighbors)
    # Inverse of the sum of distances to each grid point.
    inv_sum_dists = 1. / dists[0].sum(1)

    # Reshape
    im = inv_sum_dists.reshape(resolution, resolution)
    return im


def kNNDens(xv, yv, grid, resolution, neighbors):
    """
    Assigns the normalized number of points divided by the radius (distance to
    the most distant neighbor) to each grid point.

    """
    # Create the tree
    tree = cKDTree(np.array([xv, yv]).T)
    # Find the closest nnmax-1 neighbors (first entry is the point itself)
    dists = tree.query(grid, neighbors)

    rads = dists[0].max(1)
    # Inverse of the sum of distances to each grid point.
    grid_dens = (neighbors / xv.size) / (np.pi * rads**2)

    # Reshape
    im = grid_dens.reshape(resolution, resolution)
    return im


def data_coord2view_coord(p, resolution, pmin, pmax):
    dp = pmax - pmin
    dv = (p - pmin) / dp * resolution
    return dv


def KDEPrep(xs, ys, resolution, extent):
    """
    """
    gd_c = complex(0, resolution)

    # Perform the kernel density estimate
    xx, yy = np.mgrid[extent[0]:extent[1]:gd_c, extent[2]:extent[3]:gd_c]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([xs, ys])

    kernel = gaussian_kde(values)
    scotts_f = kernel.covariance_factor() * np.max(values.std(axis=1))

    return values, xx, yy, positions, scotts_f


def KDEDens(resolution, kde_vals, kde_pos, scotts_f, bdw_f):
    """
    """
    bw_kde = (bdw_f * scotts_f) / np.max(kde_vals.std(axis=1))
    kernel = gaussian_kde(kde_vals, bw_method=bw_kde)
    kde = kernel(kde_pos)

    # Reshape
    kde = np.reshape(kde.T, (resolution, resolution))

    return kde


def makePlot(
        file, xs, ys, resolution, extent, method, vals, kde_pars, kNN_pars):
    """
    """
    fig = plt.figure(figsize=(15, 15))
    gs = gridspec.GridSpec(2, 2)

    for gs_i, bw_nghbrs in enumerate(vals):

        if gs_i == 0:
            ax = plt.subplot(gs[gs_i])
            ax.minorticks_on()
            ax.grid(b=True, which='both', color='gray', linestyle='--', lw=.5)
            ax.plot(xs, ys, 'k.', markersize=5)
            ax.set_aspect('equal')

        else:
            ax = plt.subplot(gs[gs_i])
            ax.minorticks_on()
            if method == 'KDE':
                kde_vals, x_grid, y_grid, kde_pos, scotts_f = kde_pars
                kde = KDEDens(
                    resolution, kde_vals, kde_pos, scotts_f, bw_nghbrs)
                ax.set_title(
                    "Bandwidth: {:.4f} deg".format(scotts_f * bw_nghbrs))
                ax.imshow(np.rot90(kde), cmap=cm.Blues, extent=extent)
                ax.contour(x_grid, y_grid, kde, colors='k', linewidths=.3)
            elif method in ('kNN'):
                xv, yv, grid = kNN_pars
                if method == 'kNN':
                    im = kNNDens(xv, yv, grid, resolution, bw_nghbrs)
                elif method == 'kNNInvSum':
                    im = kNNInvSumDist(xv, yv, grid, resolution, bw_nghbrs)
                ax.imshow(im, origin='lower', extent=extent, cmap=cm.Blues)
                ax.set_title("Smoothing over {} bw_nghbrs".format(bw_nghbrs))
                ax.contour(im, colors='k', linewidths=.3, extent=extent)

        plt.xlabel("ra")
        plt.ylabel("dec")
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.invert_xaxis()

    fig.tight_layout()
    plt.savefig('{}_{}.png'.format(file, method), dpi=150, bbox_inches='tight')
    plt.clf()
    plt.close("all")


if __name__ == '__main__':
    # files = (
    #     "rgb8.gc", "rgb8.gclikeall", "rgb8.old", "rgb8.poor", "rgb8.rich",
    #     "rgb8.young")
    files = ("rgb8.gc2", "rgb8.gclike2")
    for file in files:
        main(file)
