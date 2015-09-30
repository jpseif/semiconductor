import numpy
from pylab import *
import sys
import os
import scipy.constants as Const

sys.path.append('../matterial')
sys.path.append('./silicon')

import semiconductor.matterial.ni as ni
import semiconductor.optical.silicon.opticalproperties as opticalproperties
import semiconductor.optical.absorptance as absorptance


class SpontaneousRadiativeMeission(object):

    """
    This class caculates the spectral spontaneous radiative emisison from the genralised planks law per wavelength or per photon interval
    Currents it takes in silicons properties by defult_QF_split
    """

    
    defult_QF_split = .1 * Const.e

    def __init__(self, matterial='Si',
                 optical_properties=None,
                 intrinsic_carrier_concentration=None,
                 temp = None):
        """
        Can provide a specific instance of:
            A material. then it will attempt to look up the optical constants
                intrinsic carrier concentration

            if one of
            optical properties  or  ni module is provided
            These will then be used
        """
        self.matterial = matterial

        if optical_properties is None:
            self.optics = opticalproperties.OpticalProperties(self.matterial)
        else:
            self.optics = optical_properties

        self.optics.initalise_optical_constants()

        if intrinsic_carrier_concentration is None:
            self.ni = ni.IntrinsicCarrierDensity(self.matterial)
        else:
            self.ni = intrinsic_carrier_concentration

        if temp is None:
            self.temp = 300.
        else:
            self.temp = temp

        self.optics.initalise_optical_constants()

    def blackbody_photon_per_wavelength(self, emn_wavelegnth=None, temp=None):
        """
        Returns photon emission per wavelength interval per solid angle for a black body emission at:
         temperature T (can be provide)
         for wavelengths  (entered as nm)

        sometimes is is obsrved is with a factor of 4 pi,
        this is when integrated for emission in all direction

        Currently:
         1. I divided by 10000, for some reaon? I have temp commented it out,
         2.  I have no just multiplied by 1000 for NO reason
        """

        if emn_wavelegnth is None:
            emn_wavelegnth = self.optics.wavelength_emission
        emn_wavelegnth = emn_wavelegnth * 1e-9

        if temp is None:
            temp = self.temp

        return 2 * Const.c / emn_wavelegnth**4 * 1. / (np.exp(Const.h * Const.c / emn_wavelegnth / Const.k / temp) - 1.) * 1000

    def update_tempature(self, temp):
        self.temp = temp
        self.change_temp_Green2008()
        self.ni.temp = self.temp

    def genralised_planks_PerWavelength_Carriers(self, np=None, temp=None):
        """
        generalised planks law. 
        Provides emitted photons per wavelength interval
        Uses the format outlined by green`
        """
        if temp is None:
            temp = self.temp

        # black body here is per solid angle
        BB = self.blackbody_photon_per_wavelength(temp=temp)

        # The PL spectrum with no QF splitting
        self.rsp_thermal = (BB * self.optics.alpha_BB) / self.optics.n**2

        if np is None:
            self.rsp = self.rsp_thermal
        else:
            self.rsp = self.rsp_thermal * (np) / self.ni.update_ni()**2
        # The spectrum with QF splitting

    def genralised_planks_PerEnergy(self, QF_split=False, temp=None):
        """
        generalised planks law. 
        Provides emitted photons per energy interval
        Uses the traditional form
        """
        if temp is None:
            temp = self.temp

        E = Const.h * Const.c / (self.optics.wavelength_emission * 1e-9)

        if not QF_split:
            QF_split = self.defult_QF_split

        # speed of light in medium
        c = Const.c / self.optics.n

        # Density of state of phtons
        D = E**2 / c**3 * 2**2 * Const.pi / Const.h**3

        # A note that in Gesikers phd he droped from the denumerator
        # The spectrum with QF=0 splitting
        self.rsp_thermal = self.optics.alpha_BB * c * D / (
            numpy.exp(E / Const.k / self.temp) - 1)

        # The spectrum with QF splitting
        self.rsp = self.optics.alpha_BB * c * D / (
            numpy.exp(E / Const.k / self.temp) *
            numpy.exp(-QF_split / Const.k / self.temp) - 1
        )

    def genralised_planks_PerWavelength(self, QF_split=False, temp=None):
        """
        generalised planks law. 
        Provides emitted photons per wavelength interval
        Is just an adjustedment to the energy interval expression
        """
        if temp is None:
            temp = self.temp

        # we just need to multip the per energy by the derivative below

        dEdwl = Const.h * Const.c / (self.optics.wavelength_emission * 1e-9)**2

        self.genralised_planks_PerEnergy(QF_split)

        # Adjust the values
        self.rsp *= dEdwl
        self.rsp_thermal *= dEdwl


class Simulated_PL_emission(SpontaneousRadiativeMeission):

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
    wafer_opitcs = 'polished'
    DetectionSide = 'front'

    alpha_version = 'Schinke2015'

    def __init__(self, matterial=None, optical_constants=None, excess_carriers_array=None, width = None):
        super(Simulated_PL_emission, self).__init__()

        if width is None:
            width = 0.018  # in cm


        if excess_carriers_array is None:
            self.x = np.linspace(0, width, 100)  # cm
            self.doping = 1e16                   # cm^-3
            self.np = np.ones(self.x.shape) * 1e12 * self.doping  # cm^-6

        self.Esc = absorptance.EscapeProbability(
            matterial=matterial, optical_constants=optical_constants, x=self.x)

    #     self.initalise_EmittedPL()

    def update_carrierdensity(self, deltan, doping=None):
        """
        inputs for carrier density
        If not doping is given, used the doping in self.doping
        """
        if doping is not None:
            self.doping = doping

        if deltan.shape != self.x.shape:
            print 'number of x-values not equal to delta n values'

        self.np = self.doping * deltan

    def initalise_EmittedPL(self):
        """
        Used for obtained the basic values required for this caculation
        i.e. optical cosntants, ni, an escape fraction
        """

        self.optics.initalise_optical_constants()

        # index = self.wavelength_emission > 800
        # index *= self.wavelength_emission < 1400

        # self.wavelength_emission = self.wavelength_emission[index]
        # self.optics_alpha_BB = self.optics_alpha_BB[index]
        # self.optics.n = self.optics.n[index]

        self.ni.update_ni

        # The wafter thickiness is taken as the last value in the x-direction
        self.update_escape()

    def update_temperature(self, temp=False):
        """
        Used to change the sample termpature,
        and all constants with temrpature 
        if not provided used the defult value
        """
        if not temp:
            temp = self.temp

        self.temp = temp

        # This updates the other classes, need to update
        # The optics
        # ni
        # Escape

        self.optics.change_temp_Green2008(temp)
        self.Esc.optics = self.optics
        self.ni.update_ni()
        self.update_escape()

    def update_escape(self):
        """
        Can be used to udpate the escape fraction, no inputs
        """

        self.Esc.optics = self.optics
        getattr(self.Esc, self.wafer_optics_dic[self.wafer_opitcs])()

        self.escapeprob = getattr(self.Esc,
                                  self.PL_Dection_side_depth[self.DetectionSide])

    def caculate_spectral_PL(self):
        """
        deteries the spectral PL emitted from a sample
        """
        self.genralised_planks_PerWavelength_Carriers()

        # this is the spectral distribution from each point
        # Normalised to deltan = 1, so we can just multi this by deltan

        if self.np.shape == self.x.shape:

#            print self.rsp.shape, self.escapeprob.shape, self.ni.ni, self.np.shape, '\n\n'
            self.Spectral_PL = numpy.trapz((self.rsp * self.escapeprob).T
                                           * self.np / self.ni.ni**2,
                                           self.x,
                                           axis=1)

        else:
            print 'x and np differnt lengths'

    def caculate_detected_PL(self):
        """
        multiples the detected PL by an EQE
        currently this does NOTHING
        """
        self.caculate_spectral_PL()


class Alpha_from_PL():

    """
    This class is for given a PL spectrum trying to determine the absorption coefficents
    """

    Temp = 300
    known_alpha = 1
    known_wavelength = 1
    wavelength_measured = numpy.array([1])
    PL = numpy.array([1])

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
        Used to align the realtive PL measurement to calibrated alpha at a known wavelength
        Then does a first guess at to what alpha is
        """

        self.PL_normtoBB()
        norm_PL = numpy.interp(
            self.known_wavelength, self.wavelength_measured, self.PLnBB)

        # Find the proportionality in PL for the known alpha
        self.wavelength_emission = np.array(
            (self.known_wavelength, self.known_wavelength))
        self.optics_alpha_BB = np.array((self.known_alpha, self.known_alpha))
        self.update_escape()

        # For checking the escape probability
        # plt.figure('test')
        # plt.plot(self.x, self.escapeprob[:,0])
        # plt.show()
        self.CalibrationConstant = norm_PL / self.known_alpha / \
            numpy.trapz(self.np * self.escapeprob[:, 0], self.x)

        self.PL_alpha = self.PLnBB / norm_PL * self.known_alpha

    def update_escape(self):
        """
        Caculates the escape probability given alpha
        """

        getattr(self, self.wafer_optics_dic[self.wafer_opitcs])()

        self.escapeprob = getattr(
            self, self.PL_Dection_side_depth[self.DetectionSide])

    def Iterate_alpha(self, n=10):
        """
        A function (not verified) to determine determine alpha from a PL spectrum itterativly

        Does the following itteration n times:

        1. It assumes alpha, calcs the escape probaility, caculates PL
        2. The calc PL is compared to the real PL and alpha is updated 


        a note:
        Kramers Kronig could be used here to provide a relationship for how alpha should behave
        "http://www.doria.fi/bitstream/handle/10024/96800/Conventional%20and%20nonconventional%20Kramers-Kronig%20analysis%20in%20optical%20spectroscopy.pdf?sequence=3"
        """

        for i in range(n):
            # print i
            self.wavelength_emission = self.wavelength_measured
            self.optics_alpha_BB = self.PL_alpha

            self.update_escape()
            # print self.PLnBB.shape, self.PL.shape, numpy.trapz(self.np*self.escapeprob.T,
            #                                     self.x,
            #                                     axis=1).shape
            self.PL_alpha = self.PLnBB / self.CalibrationConstant / numpy.trapz(self.np * self.escapeprob.T,
                                                                                self.x,
                                                                                axis=1)


if __name__ == "__main__":
    a = Simulated_PL_emission()
    a.initalise_EmittedPL()
    a.caculate_spectral_PL()
    a.update_escape()
    plt.plot(a.optics.wavelength_emission, a.Spectral_PL/ np.amax(a.Spectral_PL))

    a.wafer_opitcs = 'textured'
    a.update_escape()
    a.caculate_spectral_PL()
    plt.plot(a.optics.wavelength_emission, a.Spectral_PL/ np.amax(a.Spectral_PL))
    # plt.plot(a.optics.wavelength_emission, a.optics.alpha_BB)

    a.genralised_planks_PerWavelength()
    plt.plot(a.optics.wavelength_emission, a.rsp / np.amax(a.rsp), '--')

    a.genralised_planks_PerEnergy()
    plt.plot(a.optics.wavelength_emission, a.rsp / np.amax(a.rsp), '--')

    a.genralised_planks_PerWavelength_Carriers()
    plt.plot(a.optics.wavelength_emission, a.rsp / np.amax(a.rsp), ':')

    # plt.plot(a.optics.wavelength_emission,
    #          a.blackbody_photon_per_wavelength(
    #          ) / np.amax(a.blackbody_photon_per_wavelength()),
    #          '--')
    # plt.semilogy()
    plt.show()