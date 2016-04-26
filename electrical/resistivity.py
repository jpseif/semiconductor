
import numpy as np
import scipy.constants as const

from semiconductor.general_functions.carrierfunctions import get_carriers
from semiconductor.material.ni import IntrinsicCarrierDensity as ni
from mobility import Mobility as Mob
from ionisation import Ionisation as Ion
from semiconductor.helper.helper import HelperFunctions


class Resistivity(HelperFunctions):
    cal_dts = {
        'material': 'Si',
        'temp': 300,
        'mob_author': None,
        'nieff_author': None,
        'ionis_author': None,
        'dopant': 'boron',
        'Na': 1e16,
        'Nd': 0,
        'nxc': 1e10
    }

    def __init__(self, **kwargs):
        self._update_dts(**kwargs)

    def _update_links(self):

        # setting downstream values, this should change from initalisation
        # to just updating through the update function
        self.Mob = Mob(material=self.cal_dts['material'],
                       author=self.cal_dts['mob_author'],
                       temp=self.cal_dts['temp'])
        self.ni = ni(material=self.cal_dts['material'],
                     author=self.cal_dts['nieff_author'],
                     temp=self.cal_dts['temp'])
        self.ion = Ion(material=self.cal_dts['material'],
                       author=self.cal_dts['ionis_author'],
                       temp=self.cal_dts['temp'])

    def query_used_authors(self):
        return self.Mob.model, self.ni.model, self.ion.model

    def _conductivity(self, **kwargs):

        Nid, Nia = get_carriers(nxc=0,
                                Na=self.cal_dts['Na'],
                                Nd=self.cal_dts['Nd'],
                                temp=self.cal_dts['temp'],
                                ni_author=self.cal_dts['nieff_author']
                                )

        if np.all(Nid > Nia):
            Nid = self.ion.update_dopant_ionisation(
                Nid, nxc, impurity=self.cal_dts['dopant'])
        elif np.all(Nia > Nid):
            Nia = self.ion.update_dopant_ionisation(
                Nia,
                nxc=self.cal_dts['nxc'],
                impurity=self.cal_dts['dopant'])

        ne, nh = get_carriers(Nid, Nia,
                              nxc=self.cal_dts['nxc'],
                              temp=self.cal_dts['temp'],
                              ni_author=self.cal_dts['nieff_author']
                              )

        mob_e = self.Mob.electron_mobility(nxc=self.cal_dts['nxc'],
                                           Na=self.cal_dts['Na'],
                                           Nd=self.cal_dts['Nd'],
                                           temp=self.cal_dts['temp']
                                           )
        mob_h = self.Mob.hole_mobility(nxc=self.cal_dts['nxc'],
                                       Na=self.cal_dts['Na'],
                                       Nd=self.cal_dts['Nd'],
                                       temp=self.cal_dts['temp'])

        return const.e * (mob_e * ne + mob_h * nh)

    def caculate(self, **kwargs):
        '''
        calculates the resistivity
        '''

        self._update_dts(**kwargs)
        self._update_links()

        res = 1. / self._conductivity(**kwargs)

        return res
