import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from data_read import data_reader

from subhess import subhess

from selection_utils import selection_utils


class four_agb_plotter:

    def __init__(self, stage='agb_crossed', arcmin_radius=False):

        self.n147 = data_reader(stage=stage, galaxy='ngc147')
        self.n185 = data_reader(stage=stage, galaxy='ngc185')
        self.n205 = data_reader(stage=stage, galaxy='ngc205')
        self.m32 = data_reader(stage=stage, galaxy='m32')
        self.dEs = np.array([self.n147, self.n185, self.n205, self.m32])

        self.n147_cm = data_reader(stage='cm', galaxy='ngc147')
        self.n185_cm = data_reader(stage='cm', galaxy='ngc185')
        self.n205_cm = data_reader(stage='cm', galaxy='ngc205')
        self.m32_cm = data_reader(stage='cm', galaxy='m32')

        if arcmin_radius != False:
            e = selection_utils()
            for i in self.dEs:
                i.data = e.select_ellipse(i.data, afl=arcmin_radius/60)[0]

    def plot_kj_cmd_hess(self, save=False):

        def jk_errors(data, binsize, start_loc):

            bin_locs = []
            binsy = []
            binsx = []

            while start_loc > np.min(data.kmag.dropna()):
                biny = []
                binx = []
                for i in data.index:

                    if start_loc > data.kmag[i] > start_loc - binsize:
                        biny.append(data.kerr[i])
                        binx.append(np.sqrt(data.kerr[i]**2 + data.jerr[i]**2))
                biny = np.array(biny)
                binx = np.array(binx)
                avgbiny = np.average(biny)
                avgbinx = np.average(binx)
                binsy.append(avgbiny)
                binsx.append(avgbinx)
                bin_locs.append(start_loc-binsize/2)
                start_loc = start_loc-binsize

            return bin_locs, binsx, binsy

        n147 = self.n147
        n185 = self.n185
        n205 = self.n205
        m32 = self.m32

        n147_c = self.n147_cm.cdata
        n185_c = self.n185_cm.cdata
        n205_c = self.n205_cm.cdata
        m32_c = self.m32_cm.cdata

        n147_m = self.n147_cm.mdata
        n185_m = self.n185_cm.mdata
        n205_m = self.n205_cm.mdata
        m32_m = self.m32_cm.mdata

        dim_data = [[n147.data, n185.data], [n205.data, m32.data]]

        sns.set_context('paper')

        params = {'legend.fontsize': '12', 'axes.labelsize': '15',
                  'axes.titlesize': '12', 'xtick.labelsize': '10',
                  'ytick.labelsize': '10', 'lines.linewidth': 2, 'axes.linewidth': 2, 'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})

        markersize = 0.5
        marker = 'o'

        trgb_locs = [18.1, 17.86, 17.93, 17.8]
        trgb_lims = (0, 4.5)
        fore_locs = [0.992, 0.964, 1.00, 0.92]

        fore_lims = (19, 11)

        trgb_locs_dim = [[trgb_locs[0], trgb_locs[1]],
                         [trgb_locs[2], trgb_locs[3]]]
        fore_locs_dim = [[fore_locs[0], fore_locs[1]],
                         [fore_locs[2], fore_locs[3]]]

        fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=[6, 6])

        subhess.plot_kj_cmd(n147, ax=axs[0, 0], label='NGC147')
        subhess.plot_kj_cmd(n185, ax=axs[0, 1], label='NGC185')
        subhess.plot_kj_cmd(n205, ax=axs[1, 0], label='NGC205')
        subhess.plot_kj_cmd(m32, ax=axs[1, 1], label='M32')
        start_loc = 19
        binsize = 1
        err_xpos = -0.5
        for i in range(2):
            for j in range(2):

                xtrgb = trgb_lims
                ytrgb = (trgb_locs_dim[i][j], trgb_locs_dim[i][j])

                xfore = (fore_locs_dim[i][j], fore_locs_dim[i][j])
                yfore = fore_lims

                err_ypos, xerr, yerr = jk_errors(
                    dim_data[i][j], binsize=binsize, start_loc=start_loc)

                err_xpos = np.full(shape=len(err_ypos), fill_value=-0.5)

                axs[i, j].errorbar(err_xpos, err_ypos, xerr=xerr, yerr=yerr,
                                   linestyle='none', elinewidth=1, capsize=0, color='black')
                axs[i, j].plot(xtrgb, ytrgb, linestyle='--', color='black')
                axs[i, j].plot(
                    xfore, yfore, linestyle='dashdot', color='black')

        plt.subplots_adjust(wspace=0, hspace=0)

        axes = [axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]
        for i in axes:
            # i.label_outer()
            leg = i.legend(handlelength=0, handletextpad=0,
                           frameon=False, loc='upper right', markerscale=0.001)
            for item in leg.legendHandles:
                item.set_visible(False)

        for i in range(2):

            for j in range(2):

                # axs[i,j].invert_yaxis()
                axs[i, j].set_ylim(bottom=19, top=13)
                # axs[i,j].set_xlim(left=-1)
                axs[i, j].xaxis.set_major_locator(plt.MaxNLocator(5))
                axs[i, j].yaxis.set_major_locator(plt.MaxNLocator(5))
                axs[i, j].xaxis.set_minor_locator(MultipleLocator(0.5))
                axs[i, j].yaxis.set_minor_locator(MultipleLocator(0.5))
        axs[1, 0].set_xlabel('(J-K)$_0$')
        axs[1, 1].set_xlabel('(J-K)$_0$')
        axs[1, 0].set_ylabel('K$_0$')
        axs[0, 0].set_ylabel('K$_0$')

        if save == True:

            plt.savefig('4_kj_cmd_hess.png')
            plt.savefig('4_kj_cmd_hess.pdf')

        # axs[0,0].set_title('NGC147')
        # axs[0,1].set_title('NGC185')
        # axs[1,0].set_title('NGC205')
        # axs[1,1].set_title('M32')

    def plot_cc_hess(self, save=False):

        n147 = self.n147
        n185 = self.n185
        n205 = self.n205
        m32 = self.m32

        sns.set_context('paper')

        params = {'legend.fontsize': '12', 'axes.labelsize': '15',
                  'axes.titlesize': '12', 'xtick.labelsize': '10',
                  'ytick.labelsize': '10', 'lines.linewidth': 2, 'axes.linewidth': 2, 'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})

        fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=[6, 6])

        hkcut = [[0.337, 0.323], [0.407, 0.477]]
        hklim = 2.3
        jhcut = [[0.883, 0.857], [0.930, 0.913]]
        jhlim = 3.18

        subhess.plot_cc(n147, ax=axs[0, 0], label='NGC147')
        subhess.plot_cc(n185, ax=axs[0, 1], label='NGC185')
        subhess.plot_cc(n205, ax=axs[1, 0], label='NGC205')
        subhess.plot_cc(m32, ax=axs[1, 1], label='M32')

        plt.subplots_adjust(wspace=0, hspace=0)

        for i in range(2):

            for j in range(2):

                # i.label_outer()
                leg = axs[i, j].legend(
                    handlelength=0, handletextpad=0, frameon=False, loc='upper left', markerscale=0.001)
                for item in leg.legendHandles:
                    item.set_visible(False)
                #axs[i,j].xaxis.grid(True, which='minor')

                xhk = (jhcut[i][j], jhlim)
                yhk = (hkcut[i][j], hkcut[i][j])

                xjh = (jhcut[i][j], jhcut[i][j])
                yjh = (hkcut[i][j], hklim)

                axs[i, j].plot(xhk, yhk, linestyle='--', color='black')
                axs[i, j].plot(xjh, yjh, linestyle='--', color='black')

                axs[i, j].set_xlim(-0.93, jhlim)
                axs[i, j].set_ylim(0, hklim)

                axs[i, j].xaxis.set_minor_locator(MultipleLocator(0.2))
                axs[i, j].yaxis.set_minor_locator(MultipleLocator(0.2))

        if save == True:

            plt.savefig('4_cc_hess.png')
            plt.savefig('4_cc_hess.pdf')

    def plot_spatial(self, save=False):
        n147 = self.n147.data
        n185 = self.n185.data
        n205 = self.n205.data
        m32 = self.m32.data
        sns.set_context('paper')

        params = {'legend.fontsize': '12', 'axes.labelsize': '15',
                  'axes.titlesize': '12', 'xtick.labelsize': '10',
                  'ytick.labelsize': '10', 'lines.linewidth': 2, 'axes.linewidth': 2, 'animation.html': 'html5'}
        plt.rcParams.update(params)
        plt.rcParams.update({'figure.max_open_warning': 0})

        markersize = 3
        marker = 'o'

        fig, axs = plt.subplots(
            2, 2, sharex=False, sharey=False, figsize=[6, 6])

        sns.scatterplot(
            n147.xi, n147.eta, ax=axs[0, 0], marker=marker, color='k', label='NGC147')
        sns.scatterplot(
            n185.xi, n147.eta, ax=axs[0, 1], marker=marker, color='k', label='NGC185')
        sns.scatterplot(
            n205.xi, n205.eta, ax=axs[1, 0], marker=marker,  color='k', label='NGC205')
        sns.scatterplot(
            m32.xi, m32.eta, ax=axs[1, 1], marker=marker, color='k', label='M32')

        plt.subplots_adjust(wspace=0, hspace=0)

        axes = [axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]
        for i in axes:
            # i.label_outer()
            leg = i.legend(handlelength=0, handletextpad=0,
                           frameon=False, loc='upper right', markerscale=0.001)
            for item in leg.legendHandles:
                item.set_visible(False)

        for i in range(2):

            for j in range(2):

                axs[i, j].xaxis.set_major_locator(plt.MaxNLocator(5))
                axs[i, j].yaxis.set_major_locator(plt.MaxNLocator(5))
                axs[i, j].xaxis.set_minor_locator(MultipleLocator(0.5))
                axs[i, j].yaxis.set_minor_locator(MultipleLocator(0.5))
        axs[1, 0].set_xlabel(r'$\xi$')
        axs[1, 1].set_xlabel(r'$\xi$')
        axs[1, 0].set_ylabel(r'$\eta$')
        axs[0, 0].set_ylabel(r'$\eta$')

        if save == True:

            plt.savefig('4_spatial.png')
            plt.savefig('4_spatial.pdf')
