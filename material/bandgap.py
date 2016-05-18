#!/usr/local/bin/python
# UTF-8

import sys

from semiconductor.material.bandgap_intrinsic import IntrinsicBandGap
from semiconductor.material.bandgap_narrowing import BandGapNarrowing
from semiconductor.helper.helper import HelperFunctions


class BandGap(HelperFunctions):

    '''
    A simple class to combine the intrinsic band gap and
    band gap narrowing classes for easy access
    '''
    cal_dts = {
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
        self._update_dts(**kwargs)

        # pass values to models
        self._update_links()

    def _update_links(self):

        self.iEg = IntrinsicBandGap(material=self.cal_dts['material'],
                                    author=self.cal_dts['iEg_author'],
                                    temp=self.cal_dts['temp'],
                                    multiplier=self.cal_dts['multiplier'],
                                    )
        self.BGN = BandGapNarrowing(material=self.cal_dts['material'],
                                    author=self.cal_dts['BGN_author'],
                                    temp=self.cal_dts['temp'],
                                    nxc=self.cal_dts['nxc'],
                                    )

    def plot_all_models(self):
        self.iEg.plot_all_models()
        self.BGN.plot_all_models()

    def update(self, **kwargs):
        '''
        Calculates the band gap
        '''
        self._update_dts(**kwargs)

        # just prints a warning if the model is for the incorrect
        # dopants
        dopant_model_list = self.BGN.available_models(
            'dopant', self.cal_dts['dopant'])

        # check dopant and model line up
        if self.cal_dts['BGN_author'] not in dopant_model_list:
            sys.exit(
                '''\nThe BGN author you have selected was not for your'''
                ''' selected dopant.\n'''
                '''Please try selecting one of the following authors:\n\t''' +
                str('\n\t'.join([i for i in dopant_model_list])) +
                '''\nFor the selected dopant: {0}\n'''.format(
                    self.cal_dts['dopant'])
            )

        Eg = self.iEg.update(material=self.cal_dts['material'],
                             author=self.cal_dts['iEg_author'],
                             temp=self.cal_dts['temp'],
                             multiplier=self.cal_dts['multiplier'],
                             ) - \
            self.BGN.update(material=self.cal_dts['material'],
                            author=self.cal_dts['BGN_author'],
                            temp=self.cal_dts['temp'],
                            nxc=self.cal_dts['nxc'],
                            )
        return Eg

    def check_models(self):
        self.iEg.check_models()
        self.BGN.check_models()
