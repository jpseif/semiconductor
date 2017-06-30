#!/usr/local/bin/python
# UTF-8

import numpy as np
import matplotlib.pylab as plt
import os
import configparser
import scipy.constants as C

from semiconductor.helper.helper import BaseModelClass
from semiconductor.material import bandgap_narrowing_models as Bgn
from semiconductor.general_functions import carrierfunctions as GF


class BandGapNarrowing(BaseModelClass):

    '''
    Band gap narrowing accounts for a reduction in bandgap which
    occurs as a result of no thermal effects. These include:
        doping
        excess carrier density (non thermal distribution)

    This class allows calculation of an effective intrinsic carrier
    density. However, it only uses boltzman stastics.

    Inputs to this class are:

        1. material: (str)
            The elemental name for the material. Defualt (Si)
        2. temp: (float)
            The temperature of the material in Kelvin (300)
        3. author: (str)
            The author of the model to be used
        4. nxc: (array like cm^-3)
            The number of excess carriers
        5. Na: (array like cm^-3)
            The number of acceptor dopants
        6. Nd: (array like cm^-3)
            The number of donar dopants

    '''
    _cal_dts = {
        'material': 'Si',
        'temp': 300.,
        'author': None,
        'nxc': 1e10,
        'Na': 0,
        'Nd': 1e16,
    }

    author_list = 'bandgap_narrowing.yaml'

    def __init__(self, **kwargs):

        # update any values in cal_dts
        # that are passed
        self.calculationdetails = kwargs

        # get the address of the authors list
        author_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            self._cal_dts['material'],
            self.author_list)

        # get the models ready
        self._int_model(author_file)

        # initiate the first model
        self.change_model(self._cal_dts['author'])

    def update(self, **kwargs):
        '''
        Calculates the band gap narrowing

        Inputs:
        Na, Nd, delta n, temp, ni

        output:
            band gap narrowing in eV
        '''

        self.calculationdetails = kwargs

        # a check to make sure the model hasn't changed
        if 'author' in kwargs.keys():
            self.change_model(self._cal_dts['author'])

        # this should be change an outside function alter
        ne, nh = GF.get_carriers(Na=self._cal_dts['Na'],
                                 Nd=self._cal_dts['Nd'],
                                 nxc=self._cal_dts['nxc'],
                                 temp=self._cal_dts['temp'])

        doping = np.array(np.abs(self._cal_dts['Na'] - self._cal_dts['Nd']))

        return getattr(Bgn, self.model)(
            self.vals,
            Na=np.copy(self._cal_dts['Na']),
            Nd=np.copy(self._cal_dts['Nd']),
            ne=ne,
            nh=nh,
            temp=self._cal_dts['temp'],
            doping=doping)

    def ni_eff(self, ni, **kwargs):
        '''
        returns the effective intrinsic carrier densitiy
        '''
        if not isinstance(ni, np.ndarray):
            ni = np.asarray([ni])

        mult = self.ni_multiplier(**kwargs)

        if not isinstance(mult, np.ndarray):
            mult = np.asarray([mult])

        assert ni.shape == mult.shape or ni.shape[0] == 1

        return ni * mult

    def ni_multiplier(self, **kwargs):
        '''
        returns a multiplification factor that when applied to the intrinsic
        carrier concentration provides the effective intrinsic carrier
        concentraion.
        '''

        BGN = self.update(**kwargs)
        vt = C.k * self._cal_dts['temp'] / C.e
        return np.exp(BGN / vt / 2.)

    def check_models(self):
        plt.figure('Bandgap narrowing')
        Nd = 0.
        dn = 1e14
        temp = 300.

        test_file = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'Si', 'check data', 'Bgn.csv')

        data = np.genfromtxt(test_file, delimiter=',', names=True)

        for author in data.dtype.names[1:]:

            if author in self.available_models():
                Na = data['N']
                BGN = self.update(Na=Na, Nd=Nd, nxc=dn,
                                  author=author,
                                  temp=temp)

                plt.plot(
                    data['N'], data[author], '.',
                    label='PV-lighthouse\'s: ' + author)
                if not np.all(BGN == 0):
                    plt.plot(Na, BGN, label=author)

                if np.average(np.abs(data[author] - BGN) / (BGN + 0.1)) > 0.003:
                    print('{0} failed to match test data'.format(author),
                          np.average(np.abs(data[author] - BGN) / (BGN + 0.1)))

            else:
                print('{0} not an available model'.format(author))
                print('models availables are', self.available_models())

        plt.semilogx()
        plt.xlabel('Doping (cm$^{-3}$)')
        plt.ylabel('Bandgap narrowing (K)')

        plt.legend(loc=0)


def check_Schenk(fig, ax):
    '''compared to values taken from XXX'''
    BGN = BandGapNarrowing(material='Si', author='Schenk1988fer')

    folder = os.path.join(
        os.path.dirname(__file__), 'Si', r'check data')

    fnames = ['BGN_Schenk_asN-dn-1e14.csv']
    nxc = 1e14

    ax.set_color_cycle(['c', 'c', 'm', 'm', 'b', 'b', 'r', 'r', 'g', 'g'])

    for f_name in fnames:
        data = np.genfromtxt(os.path.join(folder, f_name),
                             names=True,
                             delimiter=',',
                             skip_header=1)
        ND = np.zeros(data.shape)
        for temp in data.dtype.names[1::2]:
            bgn = BGN.update(data['N'], ND, nxc, temp=float(temp))
            ax.plot(data['N'], bgn,
                    '.')
            ax.plot(data['N'], data[temp],
                    '--',
                    label=temp)

        ax.legend(loc=0, title='Temperature (K)')

    ax.set_title('BGN comparison to PV-lighthouse: $\Delta$n=1e14:')
    ax.set_ylabel('Bang gap narrowing (eV)')
    ax.set_xlabel('Ionised Doping (cm$^{-3}$)')
    ax.semilogx()

if __name__ == '__main__':

    bgn = BandGapNarrowing()
    bgn.check_models()
    plt.show()
