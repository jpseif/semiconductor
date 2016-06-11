import numpy as np
from matplotlib.pylab import *
import sys
import os
import scipy.constants as const

import semiconductor.material.ni as ni
import semiconductor.optical.opticalproperties as opticalproperties
import semiconductor.optical.absorptance as absorptance
from semiconductor.helper.helper import HelperFunctions

import inspect


class SpontaneousRadiativeEmission(HelperFunctions):

    """
    This class calculates the spectral spontaneous radiative emisison from
    the genralised planks law per wavelength or per photon interval
    Currents it takes in silicons properties by defult_QF_split
    """

    def __init__(self, **kwargs):
        """
        Can provide a specific instance of:
            A material. then it will attempt to look up the optical constants
                intrinsic carrier concentration

            if one of
            optical properties  or  ni module is provided
            These will then be used
        """
        self.material = 'Si'
        self.ni_author = None
        self.temp = 300.
        self.optics_abs_author = 'Green2008'
        self.optics_ref_author = 'Green2008'

        self._update_vars(**kwargs)
        self._update_links()

    def _update_links(self):
        '''
        updates the models used, or takes a passed class
        '''

        self._optics = opticalproperties.TabulatedOpticalProperties(
            material=self.material, temp=self.temp,
            abs_author=self.optics_abs_author,
            ref_author=self.optics_ref_author)
        self._ni = ni.IntrinsicCarrierDensity(material=self.material,
                                              author=self.ni_author,
                                              temp=self.temp)

    def blackbody_photon_per_wavelength(self, emn_wavelegnth=None, **kwargs):
        """
        Returns photon emission per wavelength interval per solid angle for a
        black body emission.

        sometimes is is obsrved is with a factor of 4 pi,
        this is when integrated for emission in all direction

        input:
            wavelength: (np array)
                    wavelength in nm
            kwargs: (optional)
                anything returned by the cal_dts function

        Currently:
         1. I divided by 10000, for some reaon? I have temp commented it out,
         2. I have no just multiplied by 1000 for NO reason
        """

        if emn_wavelegnth is None:
            emn_wavelegnth = self._optics.wavelength
        # wavelength to meters
        emn_wavelegnth = emn_wavelegnth * 1e-9

        return 2 * const.c / emn_wavelegnth**4 * 1. / (
            np.exp(const.h * const.c / emn_wavelegnth / const.k / self.temp) -
            1.) * 1000

    def genralised_planks_PerWavelength_Carriers(self, np=1e16, **kwargs):
        """
        generalised planks law.
        Provides emitted photons per wavelength interval
        Uses the format outlined by green

        inputs:
                ni (float, optional):
                    the product of carriers in cm^6
        returns:
                The black body emission spectrum
        """
        if bool(kwargs):
            self._update_vars(**kwargs)
            self._update_links()

        # black body here is per solid angle
        BB = self.blackbody_photon_per_wavelength(temp=self.temp)

        # The PL spectrum with no QF splitting
        rsp_thermal = (
            BB * self._optics.abs_cof_bb) / self._optics.ref_ind**2

        return rsp_thermal * ((np) / self._ni.update()**2)

    def genralised_planks_PerEnergy(self, QF_split=0.1, **kwargs):
        """
        generalised planks law.
        Provides emitted photons per energy interval
        Uses the traditional form

        inputs:
                QF_split (float, optional):
                    Quasi fermi energly level splitting in eV
        returns:
                The black body emission spectrum
        """

        if bool(kwargs):
            self._update_vars(**kwargs)
            self._update_links()

        QF_split *= const.e

        E = const.h * const.c / (self._optics.wavelength * 1e-9)

        # speed of light in medium
        try:
            c = const.c / self._optics.ref_ind
        except:
            c = const.c

        # Density of state of phtons
        D = E**2 / c**3 * 2**2 * const.pi / const.h**3

        # Note that in Gesikers phd he droped from the denumerator
        # The spectrum with QF splitting
        return self._optics.abs_cof_bb * c * D / (
            np.exp(E / const.k / self.temp) *
            np.exp(-QF_split / const.k / self.temp) - 1
        )

    def genralised_planks_PerWavelength(self, **kwargs):
        """
        generalised planks law.
        Provides emitted photons per wavelength interval
        Is just an adjustedment to the energy interval expression
        """

        if bool(kwargs):
            self._update_vars(**kwargs)
            self._update_links()

        # we just need to multip the per energy by the derivative below
        dEdwl = const.h * const.c / (self._optics.wavelength * 1e-9)**2

        # Adjust the values
        return self.genralised_planks_PerEnergy(**kwargs) * dEdwl


class luminescence_emission(HelperFunctions):

    """
    A class that simualted the PL emitted by a device
    i.e. tries to account for reabsorption and reflections
    """

    # Given a deltan V x profile  provides the PL out the system
    # Currently can not adjust for dector

    # Dictionaries
    wafer_optics_dic = {'polished': 'double_side_polished',
                        'textured': 'double_side_lambertian'}
    PL_Dection_side_depth = {'rear': 'Escape_rear',
                             'front': 'Escape_front'}

    def __init__(self, **kwargs):

        self.wafer_opitcs = 'polished'
        self.DetectionSide = 'front'
        # self.alpha_version = 'Schinke2015'
        self.material = 'Si'
        self.temp = 300.  # temp in kelvin
        self.width = 0.018  # width in cm
        self.ni_author = None  # author of intrinsic carrier density
        self.optics_abs_author = 'Green2008'
        self.optics_ref_author = 'Green2008'
        self.nxc = None  # the number of excess carrier with depth
        self.doping = 1e16  # the doping in cm^-3
        self._index = None

        self._update_vars(**kwargs)
        self._update_x_dist()
        self._update_links()

    def _update_x_dist(self):
        '''
        updates the distance
        '''
        self.nxc = np.ones(10)
        self._x = np.linspace(0, self.width, self.nxc.shape[0])

    def _update_links(self):
        '''
        A function where the links to other
        instances's is completely refreshed.
        '''
        self._sre = SpontaneousRadiativeEmission(
            temp=self.temp,
            optics_abs_author=self.optics_abs_author,
            optics_ref_author=self.optics_ref_author,
            material=self.material,
            ni_author=self.ni_author,
            )

        # I got lasy, so i'm using the previous classes stuff
        self._optics = self._sre._optics

        if self._index is None:
            self._index = self._optics.wavelength > 0

        self._optics.wavelength = self._optics.wavelength[self._index]
        self._optics.abs_cof_bb = self._optics.abs_cof_bb[self._index]
        self._optics.ref_ind = self._optics.ref_ind[self._index]

        self._sre._optics = self._optics

        self._esc = absorptance.EscapeProbability(
            material=self.material,
            x=self._x)

        self._esc._optics = self._optics

        self._update_escape()

    def update_carrierdensity(self, deltan, doping=None):
        """
        inputs for carrier density
        If not doping is given, used the doping in self.doping
        """
        if doping is not None:
            self.doping = doping

        if deltan.shape != self.x.shape:
            print('number of x-values not equal to delta n values')

        self.np = self.doping * deltan

    def limit_wavelegnths(self, wl_min=None, wl_max=None):
        """
        Used for obtained the basic values required for this caculation
        i.e. optical cosntants, ni, an escape fraction
        """

        self._index = self._optics.wavelength > wl_min
        self._index *= self._optics.wavelength < wl_max

    def _update_escape(self):
        """
        Can be used to update the escape fraction, no inputs
        """

        getattr(self._esc, self.wafer_optics_dic[self.wafer_opitcs])()

        self._escapeprob = getattr(
            self._esc,
            self.PL_Dection_side_depth[self.DetectionSide])

    def calculate_spectral(self, **kwargs):
        """
        deteries the spectral PL emitted from a sample

        don't think this is correct at the moment
        """
        # ensure inputs are good
        if bool(kwargs):
            self._update_vars(**kwargs)
            self._update_x_dist()
            self._update_links()

        # cacualte the generated PL
        sre = self._sre.genralised_planks_PerWavelength_Carriers(
            self._sre._ni.update()**2)

        # this is the spectral distribution from each point
        # Normalised to deltan = 1, so we can just multi this by deltan
        assert self.nxc.shape == self._x.shape, (
            "nxc is different length to x spacing")

        Spectral_PL = np.trapz(
            (sre * self._escapeprob).T * self.nxc * self.doping /
            self._sre._ni.update()**2,
            self._x,
            axis=1)

        return Spectral_PL

    def calculate_emitted(self, **kwargs):
        """
        multiples the detected PL by an EQE
        currently this does NOTHING
        """
        spectral = self.calculate_spectral(**kwargs)
        return np.trapz(spectral, self._optics.wavelength)


class Alpha_from_PL():

    """
    This class is for given a PL spectrum trying to
    determine the absorption coefficents
    """

    Temp = 300
    known_alpha = 1
    known_wavelength = 1
    wavelength_measured = np.array([1])
    PL = np.array([1])

    W = 0.018
    # Dictionaries
    wafer_optics_dic = {'polished': 'escprob_polished_schick1992',
                        'textured': 'escprob_textured_rudiger2007'}
    PL_Dection_side_depth = {'rear': 'Escape_rear',
                             'front': 'Escape_front'}
    wafer_opitcs = 'polished'
    DetectionSide = 'front'

    # this class should start by doing a point wise fit to the
    def PL_normtoBB(self):
        """
        Just a little function to remove black body stuff
        """

        BB = BlackBody.PhotonFlux(self.wavelength_measured, self.temp)

        self.PLnBB = self.PL / BB

    def Guess_alpha(self):
        """
        Used to align the realtive PL measurement to calibrated alpha at
        a known wavelength. Then does a first guess at to what alpha is.
        """

        self.PL_normtoBB()
        norm_PL = np.interp(
            self.known_wavelength, self.wavelength_measured, self.PLnBB)

        # Find the proportionality in PL for the known alpha
        self.wavelength = np.array(
            (self.known_wavelength, self.known_wavelength))
        self.optics_abs_cof_bb = np.array((self.known_alpha, self.known_alpha))
        self.update_escape()

        # For checking the escape probability
        # plt.figure('test')
        # plt.plot(self.x, self.escapeprob[:,0])
        # plt.show()
        self.Calibrationconstant = norm_PL / self.known_alpha / \
            np.trapz(self.np * self.escapeprob[:, 0], self.x)

        self.PL_alpha = self.PLnBB / norm_PL * self.known_alpha

    def update_escape(self):
        """
        calculates the escape probability given alpha
        """

        getattr(self, self.wafer_optics_dic[self.wafer_opitcs])()

        self.escapeprob = getattr(
            self, self.PL_Dection_side_depth[self.DetectionSide])

    def Iterate_alpha(self, n=10):
        """
        A function (not verified) to determine determine alpha from a PL
        spectrum itterativly

        Does the following itteration n times:

        1. It assumes alpha, calcs the escape probaility, calculates PL
        2. The calc PL is compared to the real PL and alpha is updated


        a note:
        Kramers Kronig could be used here to provide a relationship for how
        alpha should behave
        "http://www.doria.fi/bitstream/handle/10024/96800/Conventional%20and%20nonconventional%20Kramers-Kronig%20analysis%20in%20optical%20spectroscopy.pdf?sequence=3"
        """

        for i in range(n):
            # print i
            self.wavelength = self.wavelength_measured
            self.optics_abs_cof_bb = self.PL_alpha

            self.update_escape()

            self.PL_alpha = self.PLnBB / self.Calibrationconstant /\
                np.trapz(self.np * self.escapeprob.T,
                         self.x,
                         axis=1)


if __name__ == "__main__":
    a = Simulated_PL_emission()
    a.initalise_EmittedPL()
    a.calculate_spectral_PL()
    a.update_escape()
    plt.plot(a.optics.wavelength, a.Spectral_PL / np.amax(a.Spectral_PL))

    a.wafer_opitcs = 'textured'
    a.update_escape()
    a.calculate_spectral_PL()
    plt.plot(a.optics.wavelength, a.Spectral_PL / np.amax(a.Spectral_PL))
    # plt.plot(a.optics.wavelength, a.optics.abs_cof_bb)

    a.genralised_planks_PerWavelength()
    plt.plot(a.optics.wavelength, a.rsp / np.amax(a.rsp), '--')

    a.genralised_planks_PerEnergy()
    plt.plot(a.optics.wavelength, a.rsp / np.amax(a.rsp), '--')

    a.genralised_planks_PerWavelength_Carriers()
    plt.plot(a.optics.wavelength, a.rsp / np.amax(a.rsp), ':')

    # plt.plot(a.optics.wavelength,
    #          a.blackbody_photon_per_wavelength(
    #          ) / np.amax(a.blackbody_photon_per_wavelength()),
    #          '--')
    # plt.semilogy()
    plt.show()
