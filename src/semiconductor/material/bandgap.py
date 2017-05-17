#!/usr/local/bin/python
# UTF-8

import sys

from semiconductor.material.bandgap_intrinsic import IntrinsicBandGap
from semiconductor.material.bandgap_narrowing import BandGapNarrowing
from semiconductor.helper.helper import BaseModelClass


class BandGap(BaseModelClass):

    '''
    A simple class that combines the intrinsic band gap and
    band gap narrowing classes for easy access

    Inputs to this class are:

        1. material: (str)
            The elemental name for the material. Defualt (Si)
        2. temp: (float)
            The temperature of the material in Kelvin (300)
        3. iEg_author: (str)
            The author of the intrinsic band gap model to be used
        4. multiplier: (float)
            A multipler. This is a hack that people use to adjust the bandgap to achieve other desired values.
        5. BGN_author: (str)
            The author of the band gap narrowing.
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
        'iEg_author': None,
        'multiplier': 1.,
        'BGN_author': None,
        'dopant': 'boron',
        'nxc': 1,
        'Na': 1e16,
        'Nd': 0,
    }

    def __init__(self, **kwargs):
        # update any values in cal_dts
        # that are passed
        self.calculationdetails = kwargs

        # pass values to models
        self._update_links()

    def _update_links(self):

        self.iEg = IntrinsicBandGap(material=self._cal_dts['material'],
                                    author=self._cal_dts['iEg_author'],
                                    temp=self._cal_dts['temp'],
                                    multiplier=self._cal_dts['multiplier'],
                                    )
        self.BGN = BandGapNarrowing(material=self._cal_dts['material'],
                                    author=self._cal_dts['BGN_author'],
                                    temp=self._cal_dts['temp'],
                                    nxc=self._cal_dts['nxc'],
                                    )

    def plot_all_models(self):
        self.iEg.plot_all_models()
        self.BGN.plot_all_models()

    def update(self, **kwargs):
        '''
        Calculates the band gap
        '''
        self.calculationdetails = kwargs

        # just prints a warning if the model is for the incorrect
        # dopants
        dopant_model_list = self.BGN.available_models(
            'dopant', self._cal_dts['dopant'])

        # check dopant and model line up
        if self._cal_dts['BGN_author'] not in dopant_model_list:
            sys.exit(
                '''\nThe BGN author you have selected was not for your'''
                ''' selected dopant.\n'''
                '''Please try selecting one of the following authors:\n\t''' +
                str('\n\t'.join([i for i in dopant_model_list])) +
                '''\nFor the selected dopant: {0}\n'''.format(
                    self._cal_dts['dopant'])
            )

        Eg = self.iEg.update(material=self._cal_dts['material'],
                             author=self._cal_dts['iEg_author'],
                             temp=self._cal_dts['temp'],
                             multiplier=self._cal_dts['multiplier'],
                             ) - \
            self.BGN.update(material=self._cal_dts['material'],
                            author=self._cal_dts['BGN_author'],
                            temp=self._cal_dts['temp'],
                            nxc=self._cal_dts['nxc'],
                            Na=self._cal_dts['Na'],
                            Nd=self._cal_dts['Nd'],
                            )
        return Eg

    def check_models(self):
        self.iEg.check_models()
        self.BGN.check_models()
