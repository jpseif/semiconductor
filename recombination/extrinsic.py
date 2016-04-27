
import numpy as np
import sys
import os
import scipy.constants as const
import matplotlib.pylab as plt

from semiconductor.helper.helper import HelperFunctions
from semiconductor.general_functions.carrierfunctions import get_carriers
from semiconductor.material.ni import IntrinsicCarrierDensity as ni

sys.path.append(
    os.path.abspath(os.path.join(os.getcwd(), os.pardir, os.pardir)))


class SRH(HelperFunctions):

    cal_dts = {
        'material': 'Si',
        'defect': None,
        'temp': 300.,
        'vth': 2.05E+7,  # taken from PVlighthouse at 300 K
        'Nt': 1e10,  # the number of traps
        'Nd': 0,
        'Na': 1e16,
        'nxc': 1e10,
        'ni_author': None
    }

    author_list = 'SRH.defects'

    def __init__(self, **kwargs):

        # update any values in cal_dts
        # that are passed
        self._update_dts(**kwargs)

        # get the address of the authors list
        author_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            self.cal_dts['material'],
            self.author_list)

        # get the models ready
        self._int_model(author_file)

        # initiate the a defect
        self.change_model(self.cal_dts['defect'])
        self._cal_taun_taup()
        self._update_links()

    def _update_links(self):

        self.ni = ni(material=self.cal_dts['material'],
                     author=self.cal_dts['ni_author'],
                     temp=self.cal_dts['temp'],
                     )
        self.ni.update()

    def _cal_taun_taup(self):
        '''
        Determines the SRH lifetime values
        from the capture cross sections
        '''

        self.vals['tau_e'] = 1. / self.cal_dts['Nt']\
            / self.vals['sigma_e'] / self.cal_dts['vth']  # s
        self.vals['tau_h'] = 1. / self.cal_dts['Nt']\
            / self.vals['sigma_h'] / self.cal_dts['vth']  # s

    def tau(self, **kwargs):
        '''
        reports the lifetime of the current defect for the given
        excess carrier density.
        '''
        self._update_dts(**kwargs)

        if 'defect' in kwargs:
            self.change_model(self.cal_dts['defect'])
            self._cal_taun_taup()

        # if change ni, calculate it
        if 'author' in ''.join(kwargs.keys()):
            self._update_links()

        return self._tau(self.cal_dts['nxc'],
                         self.vals['tau_e'],
                         self.vals['tau_h'],
                         self.vals['et'])

    def itau(self, **kwargs):
        return 1. / self.tau(**kwargs)

    def _tau(self, nxc, tau_e, tau_h, Et):
        """
        The Shockley Read Hall lifetime for a defect level
        in the band gap. This is calculated from the kinetics
        of repairing.
        """

        # the escape from defects
        nh1 = self.ni.ni * \
            np.exp(-Et / (const.k * self.cal_dts['temp'] / const.e))
        ne1 = self.ni.ni *\
            np.exp(Et / (const.k * self.cal_dts['temp'] / const.e))

        # get the number of carriers
        ne, nh = get_carriers(Na=self.cal_dts['Na'],
                              Nd=self.cal_dts['Nd'],
                              nxc=self.cal_dts['nxc'],
                              temp=self.cal_dts['temp'],
                              material=self.cal_dts['material'],
                              ni=self.ni.ni
                              )

        # calculate the recombination rate
        U = (ne * nh - self.ni.ni**2) / \
            (tau_h * (ne + ne1) + tau_e * (nh + nh1))

        return self.cal_dts['nxc'] / U

    def _plot_all(self):
        fig, ax = plt.subplots(1, 2, figsize=(16, 6))
        # ax = plt.add_subplot(111)
        counter = 0
        defects = self.AvailableModels()
        deltan = np.logspace(12, 17)
        for defect in defects:
            # ax.plot(np.inf,np.inf,'k-',label = 'Auger')
            self.select_defect(defect)
            self.cal_tau(1e10)
            ax[0].plot(counter, self.Et, 'o', label=defect)
            ax[1].plot(deltan, self.tau(deltan) * 1e6, label=defect)
            counter += 1
        # ax.legend(loc=0)
        # ax.loglog()
        x = np.arange(counter)
        # print defects
        ax[0].set_xticks(x)
        ax[0].set_xticklabels(defects, rotation=90)
        ax[0].set_xlabel('Defect')
        ax[0].set_ylabel('E$_t$ from Ei (eV)')
        ax[1].loglog()
        ax[1].set_xlabel('$\Delta$ n (cm$^{-3}$)')
        ax[1].set_ylabel('Lifetime (us)')
