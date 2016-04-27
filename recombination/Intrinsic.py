
import numpy as np
import matplotlib.pylab as plt
import os
import ConfigParser

from semiconductor.helper.helper import HelperFunctions, change_model
from semiconductor.general_functions.carrierfunctions import get_carriers
import radiative_models as radmdls
import auger_models as augmdls


class Intrinsic(HelperFunctions):

    cal_dts = {
        'material': 'Si',
        'temp': 300.,
        'ni_author': None,
        'rad_author': None,
        'aug_author': None,
    }

    def __init__(self, **kwargs):
        # update any values in cal_dts
        # that are passed
        self._update_dts(**kwargs)

        # pass values to models
        self._update_links()

    def _update_links(self):

        self.Radiative = Radiative(
            material=self.cal_dts['material'],
            author=self.cal_dts['rad_author'],
            temp=self.cal_dts['temp'],
            ni_author=self.cal_dts['ni_author']
        )

        self.Auger = Auger(
            material=self.cal_dts['material'],
            author=self.cal_dts['aug_author'],
            temp=self.cal_dts['temp'],
            ni_author=self.cal_dts['ni_author']
        )

    def tau(self, nxc, Na, Nd, **kwargs):
        '''
        Returns the intrinsic carrier lifetime
        '''
        return 1. / self.itau(nxc, Na, Nd, **kwargs)

    def itau(self, nxc, Na, Nd, **kwargs):
        '''
        Returns the inverse of the intrinsic carrier lifetime
        '''
        self._update_dts(**kwargs)
        if 'author' in ''.join(kwargs.keys()):
            self._update_links()

        itau = self.Radiative.itau(nxc, Na, Nd, **kwargs) +\
            self.Auger.itau(nxc, Na, Nd, **kwargs)

        return itau


class Radiative(HelperFunctions):

    author_list = 'radiative.model'

    cal_dts = {
        'material': 'Si',
        'temp': 300.,
        'author': None,
        'ni_author': None
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

    def tau(self, nxc, Na, Nd, **kwargs):
        self._update_dts(**kwargs)
        self.change_model(self.cal_dts['author'])

        ne0, nh0 = get_carriers(
            Na,
            Nd,
            nxc=0,
            ni_author=self.cal_dts['ni_author'],
            temp=self.cal_dts['temp']
        )

        B = self._get_B()

        return getattr(radmdls, self.model)(
            self.vals, nxc, nh0, ne0, B, temp=self.cal_dts['temp']
        )

    def itau(self, nxc, Na, Nd, **kwargs):
        return 1. / self.tau(nxc, Na, Nd, **kwargs)

    def _get_B(self):

        # if there is a model for blow, apply it
        if 'blow_model' in self.vals.keys():
            vals, model = change_model(self.Models, self.vals['blow_vals'])

            B = getattr(radmdls, self.vals['blow_model'])(
                vals, self.cal_dts['temp']
            )

        # else use the constant value
        else:
            B = self.vals['b']
        return B


class Auger(HelperFunctions):
    author_list = 'auger.model'

    cal_dts = {
        'material': 'Si',
        'temp': 300.,
        'author': None,
        'ni_author': None
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

    def tau(self, nxc, Na, Nd, **kwargs):
        self._update_dts(**kwargs)

        if 'author' in kwargs.keys():
            self.change_model(self.cal_dts['author'])

        ne0, nh0 = get_carriers(
            Na,
            Nd,
            nxc=0,
            ni_author=self.cal_dts['ni_author'],
            temp=self.cal_dts['temp']
        )

        return getattr(augmdls, self.model)(
            self.vals, nxc, ne0, nh0, temp=self.cal_dts['temp'])

    def itau(self, nxc, Na, Nd, **kwargs):
        return 1. / self.tau(nxc, Na, Nd, **kwargs)

    def check(self, author, fig=None, ax=None):
        if ax is None:
            fig, ax = plt.subplots(1)
        self.change_model(author, self.Models)

        func = getattr(augmdls, self.model)

        getattr(augmdls, author + '_check')(self.vals, func, fig, ax)
        ax.set_xlim(left=1e13)
