from data_read import data_reader
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import seaborn as sns

# class for plotting graphs based on source density in spatial profiles
# takes galaxy and stage ('cm' for just AGB stars) as inputs


class density_graphs(data_reader):

    # create contour plot
    def make_spatial_contour(self):
        # read in numerical data
        # if you want just m or c stars, add the appropriate letter to 'data'
        # e.g to load data for just m stars, replace below line with galaxy_data=self.mdata
        galaxy_data = self.data

        # make plot from tangent coordinates
        # choose level spacing and number in logspace arguments
        sns.kdeplot(galaxy_data.xi, galaxy_data.eta,
                    levels=np.logspace(-1, 3, 30))
        # set labels

        plt.ylabel(r'$\eta$')
        plt.xlabel(r'$\xi$')
        plt.gca().invert_xaxis()

    # this is going to be a nightmare
    def make_sersic_graphs(self):
        data = self.data
        xi = data['xi']
        eta = data['eta']
        fitter = fitting.LevMarLSQFitter()
        model = models.Sersic1D()


# initialise class with chosen galaxy and stage of processing (use cm for just agb stars)
grapher = density_graphs(galaxy='ngc147', stage='cm')

grapher.make_spatial_contour()
