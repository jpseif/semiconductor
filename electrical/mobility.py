#!/usr/local/bin/python
# UTF-8
# encoding=utf8

import numpy as np
import matplotlib.pylab as plt
import os

try:
    import ConfigParser as configparser
except:
    import configparser

import semiconductor.electrical.mobilitymodels as model
from semiconductor.general_functions import get_carriers

from semiconductor.helper.helper import HelperFunctions


class Mobility(HelperFunctions):
    '''
    A class to provide the mobility in a semiconductor
    carriers
    '''
    author_list = 'mobility.models'

    cal_dts = {
        'material': 'Si',
        'temp': 300,
        'author': None,
        'Na': 1e16,
        'Nd': 0,
        'nxc': 1e10,
    }

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

        # initiate the first model
        self.change_model(self.cal_dts['author'])

    def electron_mobility(self,  **kwargs):
        '''
        returns the electron mobility

        inputs:
            kwargs: (optinal)
                any value with cal_dts, for which the mobility depends on

        output:
            The electron mobility cm^2 V^-1 s^-1
        '''
        if bool(kwargs):
            self._update_dts(**kwargs)
        return getattr(model, self.model)(
            self.vals, Na=self.cal_dts['Na'], Nd=self.cal_dts['Nd'],
            nxc=self.cal_dts['nxc'], carrier='electron',
            temp=self.cal_dts['temp'])

    def hole_mobility(self, **kwargs):
        '''
        returns the hole mobility

        inputs:
            kwargs: (optinal)
                any value with cal_dts, for which the mobility depends on

        output:
            The hole mobility cm^2 V^-1 s^-1
        '''
        if bool(kwargs):
            self._update_dts(**kwargs)
        return getattr(model, self.model)(
            self.vals, Na=self.cal_dts['Na'], Nd=self.cal_dts['Nd'],
            nxc=self.cal_dts['nxc'], carrier='hole', temp=self.cal_dts['temp'])

    def mobility_sum(self,  **kwargs):
        '''
        returns the sum of the electron and hole mobilities

        inputs:
            kwargs: (optinal)
                any value with cal_dts, for which the mobility depends on

        output:
            The sum of the electron and hole mobilities cm^2 V^-1 s^-1
        '''
        if bool(kwargs):
            self._update_dts(**kwargs)
        return self.hole_mobility() +\
            self.electron_mobility()

    def ambipolar(self, ni_author=None, **kwargs):
        '''
        returns the ambipolar mobility

        inputs:
            ni_author: (optional, str)
                an author for the intrinsic carrier density
            kwargs: (optinal)
                any value with cal_dts, for which the mobility depends on

        output:
            ambipolar mobility in cm^2 V^-1 s^-1
        '''

        if bool(kwargs):
            self._update_dts(**kwargs)

        # get the number of carriers
        ne, nh = get_carriers(
            Na=self.cal_dts['Na'],
            Nd=self.cal_dts['Nd'],
            nxc=self.cal_dts['nxc'],
            temp=self.cal_dts['temp'],
            material=self.cal_dts['material'],
            ni_author=ni_author)

        mob_h = self.hole_mobility()
        mob_e = self.electron_mobility()
        mob_ambi = mob_h*mob_e/(mob_e*ne + mob_h*nh)*(ne+nh)

    def check_models(self):
        check_klaassen()
        check_dorkel()


# these checks should not be here, but rather be in the models class
def check_klaassen():

    '''compares to values taken from www.PVlighthouse.com.au'''
    a = Mobility('Si')
    a.change_model('klaassen1992')

    print('''The model disagrees at low tempeature owing to dopant ionisation'''
           '''I am unsure if mobility should take ionisated dopants or non ionisaed'''
           '''most likley it should take both, currently it only takes one''')

    dn = np.logspace(10, 20)
    # dn = np.array([1e14])
    Nd = 1e14
    Na = 0

    folder = os.path.join(
        os.path.dirname(__file__), 'Si', 'test_mobility_files')
    fnames = [r'Klassen_1e14_dopants.dat',
              r'Klassen_1e14_temp-450.dat']
    print(os.path.isdir(folder))

    for temp, f_name in zip([300, 450], fnames):

        plt.figure('Mobility - Klaassen: Deltan at ' + str(temp))

        plt.plot(dn, a.hole_mobility(dn, Na, Nd, temp=temp),
                 'r-',
                 label='hole-here')
        plt.plot(dn, a.electron_mobility(dn, Na, Nd, temp=temp),
                 'b-',
                 label='electron-here')

        print(f_name)
        data = np.genfromtxt(os.path.join(folder, f_name), names=True)

        plt.plot(data['deltan'], data['uh'], 'b--',
                 label='hole - PV-lighthouse')
        plt.plot(data['deltan'], data['ue'], 'r--',
                 label='electron - PV-lighthouse')
        plt.legend(loc=0, title='Mobility from')

        plt.semilogx()
        plt.xlabel(r'$\Delta$n (cm$^{-3}$)')
        plt.xlabel(r'Moblity  (cm$^2$V$^{-1}$s$^{-1}$)')


def check_dorkel():
    '''compares to values taken from www.PVlighthouse.com.au'''

    a = Mobility('Si')
    a.change_model(author='dorkel1981')

    dn = np.logspace(10, 20)
    # dn = np.array([1e14])
    Nd = 1e14
    Na = 0

    folder = os.path.join(
        os.path.dirname(__file__), 'Si', 'test_mobility_files')
    # file name and temp its at
    compare = [
        ['dorkel_1e14_carriers.dat', 300],
        ['dorkel_1e14_temp-450.dat', 450],
    ]

    for comp in compare:

        plt.figure('Mobility - Dorkel: Deltan at ' + str(comp[1]))

        plt.plot(dn, a.hole_mobility(dn, Na, Nd, temp=comp[1]),
                 'r-',
                 label='hole-here')
        plt.plot(dn, a.electron_mobility(dn, Na, Nd, temp=comp[1]),
                 'b-',
                 label='electron-here')

        data = np.genfromtxt(os.path.join(folder, comp[0]), names=True)

        plt.plot(data['deltan'], data['uh'], 'b--',
                 label='hole - PV-lighthouse')
        plt.plot(data['deltan'], data['ue'], 'r--',
                 label='electron - PV-lighthouse')
        plt.legend(loc=0, title='Mobility from')

        plt.semilogx()
        plt.xlabel(r'$\Delta$n (cm$^{-3}$)')
        plt.xlabel(r'Moblity  (cm$^2$V$^{-1}$s$^{-1}$)')

if __name__ == "__main__":
    check_klaassen()
    check_dorkel()
    plt.show()
