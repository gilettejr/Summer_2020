import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from data_processor import data_processor
from selection_utils import selection_utils
import os

# import fitting and polygon modellers
from astropy.modeling import fitting
from astropy.modeling.models import Sersic1D

from shapely.geometry import Point, box, MultiPoint
from shapely.geometry.polygon import Polygon, LinearRing
import shapely.affinity

from matplotlib.patches import Ellipse, Rectangle

# class for defining agb, c and m background densities


class background_constructor(data_processor):

    # function to return constant background density from border around fov
    # only appropriate for ngc147 and 185

    def find_background_density_border(self, stars, marker='o', markersize='1',
                                       color='black', borderwidth=0.07,
                                       show_figure=False):
        # select star density to be returned
        # all agb stars
        if stars == 'agb':

            data = self.data
        # m stars
        elif stars == 'm':

            data = self.mdata
        # c stars
        elif stars == 'c':

            data = self.cdata
        # set selected set to variable, same process irrespective of stars
        border_data = data.copy()
        # define outer spatial corners of set
        corner_locs = self.data_corners(self.data)
        # convert to array
        for i in corner_locs:

            i = np.asarray(i)
        # unit vector to normalise direction of border width
        subtraction_identity = np.array([[-1, -1], [-1, 1], [1, 1], [1, -1]])
        # matrix created to construct border region
        spatial_subtraction_matrix = borderwidth * subtraction_identity
        # inner corners defined
        inner_corner_locs = corner_locs + spatial_subtraction_matrix
        # converted to tuple
        inner_corner_locs = map(tuple, inner_corner_locs)
        inner_corner_locs = tuple(inner_corner_locs)
        # method called to select all stars inside border region
        select_border = selection_utils()
        select_border.select_stars(
            border_data, inner_corner_locs[0], inner_corner_locs[2], unselect=True)
        # spatial area not in border region found
        inner_area = (inner_corner_locs[0][0]-inner_corner_locs[2]
                      [0]) * (inner_corner_locs[0][1]-inner_corner_locs[2][1])
        # total spatial area of set found
        outer_area = (corner_locs[0][0]-corner_locs[2]
                      [0]) * (corner_locs[0][1]-corner_locs[2][1])
        # find area of border region
        border_area = outer_area-inner_area
        # count number of stars in border region
        border_num = len(border_data.dropna())
        # find poisson error on count
        border_err = np.sqrt(border_num)
        # show border region on image if chosen
        if show_figure == True:

            borders = [Rectangle(xy=inner_corner_locs[2],
                                 width=inner_corner_locs[0][0] -
                                 inner_corner_locs[2]
                                 [0], height=inner_corner_locs[0][1]-inner_corner_locs[2][1],
                                 fill=False, color='red')]

            plt.rc('axes', labelsize=20)
            fig, ax = plt.subplots()
            ax.plot(data.xi, data.eta, linestyle='none', marker=marker,
                    markersize=markersize, color=color, zorder=1)
            for i in borders:

                ax.add_artist(i)
            ax.invert_xaxis()
            ax.set_ylabel(r'$\eta$')
            ax.set_xlabel(r'$\xi$')
        # returns numerical density value and error in n/deg^2
        return np.array([border_num/border_area,
                         border_err/border_area, border_num])

    # method to produce background density function for ngc205 or m32
    def find_background_grad(self, stars, ellipse_a=0.15, ellipticity=0.43,
                             clockrot=169.2, marker='o', markersize='1',
                             color='black', show_figure=False, binwidth=0.02):
        # same process carried out for c, m or agb stars, defined as argument
        self.stars = stars
        # default arguments for ngc205, reset here if m32 chosen
        if self.galaxy == 'm32':

            ellipticity = 0.14
            clockrot = 157.9
            ellipse_a = 0.1
        # choose appropriate dataset
        if stars == 'agb':

            data = self.data

        elif stars == 'm':

            data = self.mdata

        elif stars == 'c':

            data = self.cdata

        # galaxy stars inside defined ellipse masked
        e = selection_utils()

        dE_stars = e.select_ellipse(
            data, afl=ellipse_a, eccentricity=ellipticity, clockrot=clockrot, unselect=True)
        dE_ellipse = dE_stars[2]
        # only backgroud stars chosen
        data = dE_stars[0]

        # convert m31 coords to local tangent coords

        m31ra = 10.68470833
        m31dec = 41.26875

        if self.galaxy == 'ngc205':
            tra = 10.09189356
            tdec = 41.68541564
        elif self.galaxy == 'm32':
            tra = 10.6742708
            tdec = 40.8651694
        # m31 centre defined from centre of chosen galaxy
        m31xi, m31eta = self.eq_to_tan(m31ra, m31dec, tra, tdec)

        # find angle of rotation so that m31 gradient is on x axis
        theta = np.arctan((m31eta)/(m31xi))
        # xi and eta coordinates retrieved from data files
        xi = data.xi.copy()
        eta = data.eta.copy()

        # rotated coordinates constructed
        data['alpha'] = self.rotate_coords(xi, eta, theta)[0]
        data['beta'] = self.rotate_coords(xi, eta, theta)[1]
        # m31 coordinates found in terms of alpha and beta
        # conversion from radians required for imported coordinates
        m31xi = np.degrees(m31xi)
        m31eta = np.degrees(m31eta)

        m31_alpha = self.rotate_coords(m31xi, m31eta, theta)[0]
        m31_beta = self.rotate_coords(m31xi, m31eta, theta)[1]

        # seaborn formatting set
        sns.set_context('paper')

        params = {'legend.fontsize': '12', 'axes.labelsize': '15',
                  'axes.titlesize': '12', 'xtick.labelsize': '10',
                  'ytick.labelsize': '10', 'lines.linewidth': 2, 'axes.linewidth': 2, 'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})
        # plot image demonstrating rotation and new coordinate system
        # if chosen
        if show_figure == True:

            plt.plot(data.alpha, data.beta, marker=marker, markersize=markersize,
                     linestyle='none', color=color, label=self.galaxy.upper())
            if self.galaxy == 'm32':

                plt.xlim(0.3, -0.3)
                plt.ylim(-0.3, 0.3)

            else:

                plt.gca().invert_xaxis()
            plt.gca().set_ylabel(r'$\beta$ (degrees)')
            plt.gca().set_xlabel(r'$\alpha$ (degrees)')
        # define vertices of rotated set
        corner1 = np.array([np.min(data.xi), np.min(data.eta)])
        corner2 = np.array([np.min(data.xi), np.max(data.eta)])
        corner3 = np.array([np.max(data.xi), np.max(data.eta)])
        corner4 = np.array([np.max(data.xi), np.min(data.eta)])

        corners = [corner1, corner2, corner3, corner4]
        corners_rot = []
        for i in corners:

            corners_rot.append(self.rotate_coords(i[0], i[1], theta))
        # construct ring encircling dataset
        border = LinearRing(corners_rot)
        # create polygon containing all stars
        s = Polygon(border)
        corners_x = np.array(
            [corners_rot[0][0], corners_rot[1][0], corners_rot[2][0], corners_rot[3][0]])
        # define lower and upper bounds of set in alpha coordinates
        upper_bound = np.max(corners_x)-binwidth
        lower_bound = np.min(corners_x)

        # construct thin vertical bins

        bin_nums = []
        bin_areas = []
        bin_locs = []
        # define masked ellipse
        mask = dE_ellipse
        theta_deg = np.degrees(theta)
        # rotate galaxy masking ellipse for new coordinates
        # required to subtract bin and mask overlap regions
        mask = shapely.affinity.rotate(dE_ellipse, -theta_deg)
        # outer ellipse ring defined
        x, y = mask.exterior.xy
        # list to contain bin shapes
        bin_shapes = []
        # while loop to construct bins across whole alpha axis
        while upper_bound-lower_bound > binwidth:
            # bin width set by extended polygon of correct width
            bin_slice_fill = box(upper_bound-binwidth, -1.0, upper_bound, 1.0)
            # polygon changed to 1D line
            bin_slice = bin_slice_fill.exterior

            # find intercept to trim the vertical extent of bins
            bin_corners = border.intersection(bin_slice)
            # intersection points reordered for processing
            swap0 = bin_corners[0]
            swap1 = bin_corners[1]
            swap2 = bin_corners[3]
            swap3 = bin_corners[2]

            corners = [swap0, swap1, swap2, swap3]
            # intersection points represent vertices of final bin
            bin_corners = MultiPoint([swap0, swap1, swap2, swap3])
            # bin shape set as polygon
            final_bin = Polygon(bin_corners)
            # polygon added to list
            bin_shapes.append(final_bin)
            # bin exterior coordinates found
            x, y = final_bin.exterior.xy
        # if chosen, figure showing bins displayed
            if show_figure == True:

                plt.plot(x, y, color='orange')

                leg = plt.legend(handlelength=0, handletextpad=0,
                                 frameon=False, loc='upper right', markerscale=0.001)
                for item in leg.legendHandles:
                    item.set_visible(False)

                plt.savefig(self.galaxy + '_background_bins.png')
                plt.savefig(self.galaxy + '_background_bins.pdf')
            # define area
            final_bin_area = final_bin.area

            # count number of stars in each bin
            box_data = data.copy()

            for i in box_data.index:

                if final_bin.contains(Point(box_data.alpha[i], box_data.beta[i])) == False:

                    box_data.loc[i] = np.nan

            bin_nums.append(len(box_data.dropna()))

            # subtract area overlapping with ellipse

            if bin_slice_fill.intersects(mask) == True:

                overlap_area = mask.intersection(bin_slice_fill).area

                final_bin_area = final_bin_area-overlap_area
            # define alpha coordinates of bin
            bin_locs.append(upper_bound-(binwidth/2))
            # define area of bin
            bin_areas.append(final_bin_area)
            # increment bound so loop creates adjacent bin
            upper_bound = upper_bound-binwidth

            # criterion for while loop slightly different for m32
            # loop broken if below criterion reached
            if self.galaxy == 'm32' and upper_bound-lower_bound < (binwidth * 2):

                break
        # convert lists into arrays
        bin_nums = np.array(bin_nums)
        bin_locs = np.array(bin_locs)
        bin_areas = np.array(bin_areas)
        # define poisson uncertainties
        bin_uncs = np.sqrt(bin_nums)
        # define densities in n/arcsec
        bin_densities = (bin_nums/bin_areas)/3600
        # convert uncertainties
        bin_uncs = (bin_uncs/bin_areas)/3600
        # convert locations to arcmins
        bin_locs = bin_locs * 60

        m31_alpha = m31_alpha*60
        # set bin locations with horizontal translation to m31 basis
        # m31 at alpha=0
        m31_bin_locs = bin_locs-m31_alpha
        # create dataframe with useful background data
        background = pd.DataFrame({'bin_densities': bin_densities, 'bin_uncs': bin_uncs,
                                  'bin_locs': bin_locs, 'm31_bin_locs': m31_bin_locs})
        # set as class attributes for retrieval
        self.background = background

        self.bin_shapes = bin_shapes
        self.stars = stars
        # fit with Sersic profile
    # method to fit previously defined background with sersic profile

    def fit_close_background(self, marker='o', markersize='3', color='black'):
        # retrieve background dataframe
        background = self.background.copy()
        # x and y data set as bin locations in alpha, density respectively
        xdata = -background.m31_bin_locs

        ydata = background.bin_densities
        # error also found
        yerr = background.bin_uncs
        # models set with appropriate initial estimates depending on
        # stars and galaxies
        # overcrowded bins skipped for better fit
        # fit carried out
        model = Sersic1D(r_eff=2, n=10.00)

        fit = fitting.LevMarLSQFitter()

        if self.galaxy == 'ngc205' and self.stars == 'c':

            model = Sersic1D(r_eff=10, n=20)
            s = fit(model, xdata[2:], ydata[2:])

        elif self.galaxy == 'm32' and self.stars == 'c':
            model = Sersic1D(r_eff=20, n=10)
            s = fit(model, xdata[3:], ydata[3:])
        elif self.galaxy == 'm32':

            model = Sersic1D(r_eff=20, n=2)
            s = fit(model, xdata[3:], ydata[3:], weights=1/yerr[3:])

        else:

            s = fit(model, xdata[2:], ydata[2:], weights=1/yerr[2:])

        # fit and data plotted with uncertainties
        plt.figure()

        sat_xdata = background.bin_locs

        plt.gca().set_xlabel(r'$\alpha$ (arcmins)')
        plt.gca().set_ylabel('Density (N/arcmin$^2$)')
        plt.gca().invert_xaxis()
        plt.errorbar(sat_xdata, ydata, yerr=yerr, capsize=2,
                     linestyle='none', marker=marker, color=color)

        plt.plot(sat_xdata, s(xdata), label=self.galaxy.upper())

        leg = plt.legend(handlelength=0, handletextpad=0,
                         frameon=False, loc='upper right', markerscale=0.001)
        for item in leg.legendHandles:
            item.set_visible(False)

        plt.savefig(self.galaxy + self.stars + '_background_fit.png')
        plt.savefig(self.galaxy + self.stars + '_background_fit.pdf')

        # sersic fit added to background dataframe
        background['density_fit'] = s(xdata)

        self.background = background

    # function subtracts sersic fitted background to produce corrected
    # density profile for chosen stars

    def find_close_slice_profile(self):

        stars = self.stars

        print(stars)

        print('Make sure youve got a background dataframe for subtraction!')
        print('You must have run the make_close_background function to create this')
        # remake alpha coordinates
        slice_shapes = self.slice_shapes
        background = self.background
        rotate_coords = self.rotate_coords

        if stars == 'agb':

            data = self.data

        elif stars == 'm':

            data = self.mdata

        elif stars == 'c':

            data = self.cdata

        areas = self.areas
        slices = self.slices
        bin_shapes = self.bin_shapes

        outer_rad = self.outer_rad
        a_width = self.a_width

        m31ra = 10.68470833
        m31dec = 41.26875

        if self.galaxy == 'ngc205':
            tra = 10.09189356
            tdec = 41.68541564

        elif self.galaxy == 'm32':
            tra = 10.6742708
            tdec = 40.8651694

        m31xi, m31eta = self.eq_to_tan(m31ra, m31dec, tra, tdec)

        # find angle of rotation
        theta = np.arctan((m31eta)/(m31xi))

        xi = data.xi.copy()
        eta = data.eta.copy()
        theta_deg = np.degrees(theta)
        # rotate

        for i in slices:
            i['alpha'] = rotate_coords(xi, eta, theta)[0]
            i['beta'] = rotate_coords(xi, eta, theta)[1]
        slice_shapes_rot = []
        for i in slice_shapes:

            slice_shapes_rot.append(shapely.affinity.rotate(i, -theta_deg))

        slice_nums = []

        for i in slices:

            slice_nums.append(len(i.dropna()))

        slice_nums = np.array(slice_nums)
        areas = np.array(areas)

        slice_densities_deg = slice_nums/areas
        slice_densities_unc = np.sqrt(slice_nums)
        slice_densities_min = slice_densities_deg/3600
        slice_densities_unc = (slice_densities_unc/areas)/3600
        # background columns: bin_densities, bin_uncs, bin_locs
        # areas in degrees for ease, background in mins
        new_density = []
        new_uncs = []
        for i in range(len(slice_shapes_rot)):
            slice_bin_densities = []
            slice_bin_densities_uncs = []
            slice_bin_areas = []
            for j in background.index:

                if slice_shapes_rot[i].intersects(bin_shapes[j]) == True:

                    overlap_shape = slice_shapes_rot[i].intersection(
                        bin_shapes[j])
                    overlap_area = overlap_shape.area

                    # subtract background
                    overlap_density = slice_densities_min[i] - \
                        background.density_fit[j]
                    overlap_density_unc = np.sqrt(
                        slice_densities_unc[i]**2 + background.bin_uncs[j]**2)

                    slice_bin_densities.append(overlap_density)
                    slice_bin_densities_uncs.append(overlap_density_unc)
                    slice_bin_areas.append(overlap_area)

            slice_bin_areas = np.array(slice_bin_areas)
            slice_bin_densities = np.array(slice_bin_densities)

            weightings = slice_bin_areas/slice_shapes_rot[i].area

            avg_density = np.average(slice_bin_densities, weights=weightings)
            density_err = (
                np.sqrt(np.sum(np.square(slice_bin_densities_uncs))))/len(slice_bin_densities)
            new_density.append(avg_density)
            new_uncs.append(density_err)

        xdata = np.linspace(outer_rad-a_width/2, 0 +
                            a_width/2, num=len(new_density))
        xdata = xdata*60

        ydata = new_density
        yerr = new_uncs

        density_distribution = pd.DataFrame(
            {'a': xdata, 'density': ydata, 'density_err': yerr, 'slice_area_deg': areas, 'slice_nums': slice_nums})

        outfilename = self.galaxy + self.stars

        print(outfilename)

        try:
            density_distribution.to_parquet(
                'unfit_background_corrected_profiles/' + outfilename)

        except:

            os.system('mkdir unfit_background_corrected_profiles')
            density_distribution.to_parquet(
                'unfit_background_corrected_profiles/' + outfilename)
