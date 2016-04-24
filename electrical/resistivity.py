
import numpy as np
import scipy.constants as const

from semiconductor.general_functions.carrierfunctions import get_carriers
from semiconductor.matterial.ni import IntrinsicCarrierDensity as ni
from mobility import Mobility as Mob
from ionisation import Ionisation as Ion


class Resistivity():
    cal_dts = {
        'matterial': 'Si',
        'temp': 300,
        'mob_author': None,
        'nieff_author': None,
        'ionis_author': None,
        'dopant': 'boron',
    }

    def _update_dts(self, **kwargs):
        '''
        assignes the inputted values that are requrired,
        befor calling a function to pass it to the downstream
        classes
        '''

        items = [i for i in kwargs.keys() if i in self.cal_dts.keys()]
        for item in items:
            self.cal_dts[item] = kwargs[item]

        self._update_links()

    def __init__(self, matterial='Si', mob_author=None,
                 nieff_author=None, ionis_author=None, temp=300,
                 **kwargs):

        temp = locals().copy()
        del temp['self']

        self._update_dts(**temp)

    def _update_links(self):

        # setting downstream values, this should change from initalisation
        # to just updating through the update function
        self.Mob = Mob(matterial=self.cal_dts['matterial'],
                       author=self.cal_dts['mob_author'],
                       temp=self.cal_dts['temp'])
        self.ni = ni(matterial=self.cal_dts['matterial'],
                     author=self.cal_dts['nieff_author'],
                     temp=self.cal_dts['temp'])
        self.ion = Ion(matterial=self.cal_dts['matterial'],
                       author=self.cal_dts['ionis_author'],
                       temp=self.cal_dts['temp'])

    def query_used_authors(self):
        return self.Mob.model, self.ni.model, self.ion.model

    def _conductivity(self, Na, Nd, nxc, **kwargs):

        Nid, Nia = get_carriers(Na, Nd, 0,
                                temp=self.cal_dts['temp'],
                                ni_author=self.cal_dts['nieff_author']
                                )

        if np.all(Nid > Nia):
            Nid = self.ion.update_dopant_ionisation(
                Nid, nxc, self.cal_dts['dopant'])
        elif np.all(Nia > Nid):
            Nia = self.ion.update_dopant_ionisation(
                Nia, nxc, self.cal_dts['dopant'])

        ne, nh = get_carriers(Nid, Nia, nxc,
                              temp=self.cal_dts['temp'],
                              ni_author=self.cal_dts['nieff_author']
                              )

        mob_e = self.Mob.electron_mobility(nxc, Na, Nd,
                                           temp=self.cal_dts['temp'])
        mob_h = self.Mob.hole_mobility(nxc, Na, Nd,
                                       temp=self.cal_dts['temp'])

        print mob_h, mob_e, Na

        return const.e * (mob_e * ne + mob_h * nh)

    def caculate(self, Na, Nd, nxc, **kwargs):
        '''
        caculates the resistivity
        '''

        self._update_dts(**kwargs)
        res = 1. / self._conductivity(Na, Nd, nxc, **kwargs)

        return res
