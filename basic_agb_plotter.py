import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from data_read import data_reader


def set_graphing_visuals():

    sns.set_context('paper')

    params = {'legend.fontsize': '12', 'axes.labelsize': '18',
              'axes.titlesize': '14', 'xtick.labelsize': '12',
              'ytick.labelsize': '12', 'lines.linewidth': 2, 'axes.linewidth': 2, 'animation.html': 'html5'}
    plt.rcParams.update(params)
    plt.rcParams.update({'figure.max_open_warning': 0})


class basic_agb_plotter(data_reader):

    def plot_kj_cmd(self, marker='o', markersize=1, color='black', newfig=False):

        galaxy_data = self.data
        # conditional statements plot only c,m, or both sets depending on
        # optional stars argument

        # axes, figure set, CMD plotted

        if newfig == True:

            plt.figure()
        plt.rc('axes', labelsize=15)
        plt.plot(galaxy_data.jmag - galaxy_data.kmag, galaxy_data.kmag,
                 linestyle='none', markersize=markersize, marker=marker, color=color)

        if color != 'red':
            plt.gca().invert_yaxis()

        plt.ylabel('$K_0$')
        plt.xlabel('$J_0$-$K_0$')

    def plot_hk_cmd(self, marker='o', markersize=1, color='black', newfig=False):
        galaxy_data = self.data

        # conditional statements plot only c,m, or both sets depending on
        # optional stars argument

        # axes, figure set, CMD plotted

        # plt.figure()
        plt.rc('axes', labelsize=15)
        plt.plot(galaxy_data.hmag - galaxy_data.kmag, galaxy_data.kmag,
                 linestyle='none', markersize=markersize, marker=marker, color=color)

        if color != 'red':
            plt.gca().invert_yaxis()

        plt.ylabel('$K_0$')
        plt.xlabel('$H_0$-$K_0$')

    def plot_spatial(self, marker='o', markersize=1, color='black', show_background=False):
        galaxy_data = self.data
    # conditional statement based on stars defines data appropriately

        # if c+m cnot chosen for stars argument, data plotted in single plot

        plt.rc('axes', labelsize=20)
        plt.plot(galaxy_data.xi, galaxy_data.eta, linestyle='none',
                 marker=marker, markersize=markersize, color=color)
        plt.gca().set_ylabel(r'$\eta$')
        plt.gca().set_xlabel(r'$\xi$')

        if color != 'red':
            plt.gca().invert_xaxis()

    def plot_contour(self, overlay=False):
        galaxy_data = self.data
        # set graphing visuals

        # make plot
        sns.kdeplot(galaxy_data.xi, galaxy_data.eta,
                    levels=np.logspace(-1, 1, 50))
        # don't reinvert x axis and relabel if being used as an overlay
        if overlay == False:

            plt.ylabel(r'$\eta$')
            plt.xlabel(r'$\xi$')
            plt.gca().invert_xaxis()

    def plot_lum(self):
        galaxy_data = self.data

        lum_func = self.data.kmag

        plt.figure()

        sns.kdeplot(lum_func.dropna(), color='black')

        plt.xlabel('$K_0$')
        plt.gca().invert_xaxis()

    # graphing method for plotting h-k/j-h 2CD


def set_graphing_visuals():

    sns.set_context('paper')

    params = {'legend.fontsize': '12', 'axes.labelsize': '18',
              'axes.titlesize': '14', 'xtick.labelsize': '12',
              'ytick.labelsize': '12', 'lines.linewidth': 2, 'axes.linewidth': 2, 'animation.html': 'html5'}
    plt.rcParams.update(params)
    plt.rcParams.update({'figure.max_open_warning': 0})
