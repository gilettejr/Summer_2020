from data_read import data_reader
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import seaborn as sns


class density_graphs(data_reader):
    def make_spatial_contour(self):

        galaxy_data = self.data
        # set graphing visuals

        # make plot
        sns.kdeplot(galaxy_data.xi, galaxy_data.eta,
                    levels=np.logspace(-1, 3, 30))
        # don't reinvert x axis and relabel if being used as an overlay

        plt.ylabel(r'$\eta$')
        plt.xlabel(r'$\xi$')
        plt.gca().invert_xaxis()

    def make_sersic_graphs(self):
        data = self.data
        xi = data['xi']
        eta = data['eta']
        fitter = fitting.LevMarLSQFitter()
        model = models.Sersic1D()


grapher = density_graphs(galaxy='ngc147', stage='cm')
grapher.make_spatial_contour()
