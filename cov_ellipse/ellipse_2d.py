
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import glob


def main(nstd=2.):
    """
    Generate an `nstd` sigma ellipse based on the mean and covariance of a
    point "cloud".
    """

    # # Generate some random, correlated data
    # cov = np.random.uniform(0., 1., (2, 2))
    # N = 50
    # points = np.random.uniform(0., 1., (N, 2))
    # # np.random.multivariate_normal(mean=(0, 0), cov=cov, size=N)

    files = glob.glob('*.txt')
    xy_cols = "RA(deg)", "Dec(deg)"
    color_col_name = 'LogAge'

    for file_name in files:
        print(file_name)

        data = ascii.read(file_name)
        x, y = data[xy_cols[0]], data[xy_cols[1]]
        color_col = data[color_col_name]
        points = np.array([x, y]).T

        # Bootstrap uncertainties in theta
        N_pts = len(points)
        N_btstrp = 10000
        theta_btstrp, ecc_btstrp = np.empty(N_btstrp), np.empty(N_btstrp)
        for i in range(N_btstrp):
            msk = np.random.choice(range(N_pts), size=N_pts, replace=True)
            vals_b, vecs_b = eigsorted(points[msk])
            theta_btstrp[i] = thetaAngle(vecs_b)
            w, h = cov_ellipse(vals_b, nstd)
            ecc_btstrp[i] = np.sqrt(1 - h**2 / w**2)
        theta_err = np.std(theta_btstrp)
        ecc_err = np.std(ecc_btstrp)

        # METHOD 1
        # Eigenvalues and eigenvectors of the covariance matrix.
        vals, vecs = eigsorted(points)
        # Obtain rotation angle from the eigenvectors
        theta = thetaAngle(vecs)
        # Obtain a, b, from the eigenvalues
        width1, height1 = cov_ellipse(vals, nstd)
        ecc = np.sqrt(1 - height1**2 / width1**2)

        # # METHOD 2
        # width2, height2 = cov_ellipse2(vals, nstd)
        # width2, height2 = 2. * width2, 2. * height2

        # Location of the center of the ellipse.
        mean_pos, stddev = points.mean(axis=0), points.std(axis=0)

        # Plot the raw points.
        x, y = points.T

        fig = plt.figure(figsize=(8, 8))
        ax = plt.gca()
        plt.title(
            r"$\theta={:.3f}\pm{:.3f}$, $e={:.3f}\pm{:.3f}$".format(
                theta, theta_err, ecc, ecc_err))

        ax.minorticks_on()
        ax.grid(b=True, which='both', color='gray', linestyle='--', lw=.5)

        # Width and height are "full" widths, hence the '2. *'
        width1, height1 = 2. * width1, 2. * height1
        im = plt.scatter(x, y, s=10, alpha=1, c=color_col)
        # First ellipse
        ellipse1 = Ellipse(
            xy=mean_pos, width=width1, height=height1, angle=theta,
            edgecolor='b', fc='None', ls='--', lw=2, zorder=4)
        ax.add_patch(ellipse1)

        # # Second ellipse
        # ellipse2 = Ellipse(
        #     xy=mean_pos, width=width2, height=height2, angle=theta,
        #     edgecolor='r', fc='None', lw=1, zorder=4)
        # ax.add_patch(ellipse2)

        plt.scatter(*mean_pos, marker='x', c='r')
        plt.xlabel("R.A. (deg)")
        plt.ylabel("Dec. (deg)")
        # plt.xlim(mean_pos[0] - 5. * stddev[0], mean_pos[0] + 5. * stddev[0])
        # plt.ylim(mean_pos[1] - 5. * stddev[1], mean_pos[1] + 5. * stddev[1])
        plt.xlim(1., 0.)
        plt.ylim(-15.8, -15.1)
        cbar = plt.colorbar(im, pad=.01, fraction=.02, aspect=40)
        cbar.ax.tick_params(labelsize=10)
        cbar.set_label(r'$Log_{10}[age]$', size=10)

        fig.tight_layout()
        plt.savefig(file_name[:-4] + '.png', dpi=300, bbox_inches='tight')


def eigsorted(points):
    '''
    Eigenvalues and eigenvectors of the covariance matrix.
    '''
    # The 2x2 covariance matrix to base the ellipse on.
    cov = np.cov(points, rowvar=False)

    # Eigenvalues and eigenvectors of the covariance matrix.
    vals, vecs = np.linalg.eigh(cov)
    # Sort putting the *largest* eigenvalue first.
    order = vals.argsort()[::-1]

    return vals[order], vecs[:, order]


def thetaAngle(vecs):
    """
    Both approaches give the same results
    """
    # theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    theta = np.degrees(np.arctan2(*vecs[::-1, 0]))

    return theta


def cov_ellipse(vals, nstd):
    """
    A N-sigma ellipse.

    Source: http://stackoverflow.com/a/12321306/1391441
    """

    # The eigenvalues are the variance of the [dimension associated to that
    # eigenvector]?.
    width, height = np.sqrt(vals) * nstd

    return width, height


def cov_ellipse2(vals, nstd):
    """
    Confidence interval?
    Source: https://stackoverflow.com/a/39749274/1391441
    """
    from scipy.stats import norm, chi2

    # Confidence level
    # Using this line, the ellipses are different
    # q = 2. * norm.cdf(nstd) - 1.
    # Using this line, the ellipses are equal
    q = 1. - np.exp(-(nstd ** 2 / 2.))

    k = np.sqrt(chi2.ppf(q, 2))
    # Gives the same result as the line above
    # k = np.sqrt(-2. * np.log(1. - q))

    width, height = np.sqrt(vals) * k

    return width, height


if __name__ == '__main__':
    main()
