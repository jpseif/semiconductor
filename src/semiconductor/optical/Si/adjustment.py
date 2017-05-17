
import numpy as np
import matplotlib.pylab as plt
import scipy.constants as const


def EnergyToWavelength(data):
    '''
    appends array with the wavelength data
    data is required to have the headings:
    'energy' and 'alpha'
    '''

    assert 'wavelength' not in data.dtype.names
    assert 'energy' in data.dtype.names
    # created the names and type for the named array
    names = list(data.dtype.names + ('wavelength',))
    fmat = ['f4' for i in names]

    datanew = np.zeros((data.shape[0]),
                       dtype={'names': names, 'formats': fmat}
                       )

    # assignes the values
    for name in data.dtype.names:
        datanew[name] = data[name]

     # cacls and appends the wavelength values
    datanew['wavelength'] = const.c * const.h / data['energy'] / const.e * 1e9

    # returns the new array
    return datanew


def wl2nrg(data):
    '''
    appends array with the wavelength data
    data is required to have the headings:
    'wavelength' and 'alpha'
    '''

    assert 'energy' not in data.dtype.names
    # created the names and type for the named array
    names = list(data.dtype.names + ('energy',))
    fmat = ['f8' for i in names]
    print(names)

    datanew = np.zeros((data.shape[0]),
                       dtype={'names': names, 'formats': fmat}
                       )

    # assignes the values
    for name in data.dtype.names:
        datanew[name] = data[name]

     # cacls and appends the wavelength values
    datanew['energy'] = (const.c * const.h) / \
        data['wavelength'] * 1e9 / const.e

    # returns the new array
    return datanew


def rmfield(a, *fieldnames_to_remove):
    return a[[name for name in a.dtype.names if name not in fieldnames_to_remove]]


def add_wavelength_oren_ergy_to_file(fname):

    data = np.genfromtxt(
        fname, names=True, delimiter=',', filling_values=np.nan)

    data = rmfield(data, 'f0')

    if 'wavelength' in data.dtype.names and 'energy' not in data.dtype.names:
        print('added energy')
        data = wl2nrg(data)
    elif 'wavelength' not in data.dtype.names and 'energy' in data.dtype.names:
        print('added wavelength')
        data = EnergyToWavelength(data)

    # a = tuple(b.split(';'))
    # data.dtype.names = a

    np.savetxt(fname, data, delimiter=',', header=','.join(
        data.dtype.names), comments='', fmt=('%0.4e'), )


def add_alpha_to_file(fname):

    data = np.genfromtxt(
        fname, names=True, delimiter=',', filling_values=np.nan)

    data = rmfield(data, 'alpha')

    # if 'wavelength' in data.dtype.names and 'energy' not in data.dtype.names:
    #     data = wl2nrg(data)
    # elif 'wavelength' not in data.dtype.names and 'energy' in data.dtype.names:
    #     data = EnergyToWavelength(data)

    data = k2alpha(data)

    # a = tuple(b.split(';'))
    # data.dtype.names = a

    np.savetxt(fname, data, delimiter=',', header=','.join(
        data.dtype.names), comments='', fmt=('%0.4e'), )


def k2alpha(data):
    '''
    appends array with the wavelength data
    data is required to have the headings:
    'energy' and 'alpha'
    '''

    assert 'alpha' not in data.dtype.names
    assert 'wavelength' in data.dtype.names

    # created the names and type for the named array
    names = list(data.dtype.names)
    for name in data.dtype.names:
        if 'k_' in name:
            suffix = name.strip('k')
            print(name, 'alpha' + suffix)
            names.append('alpha' + suffix)

    fmat = ['f4' for i in names]

    datanew = np.zeros((data.shape[0]),
                       dtype={'names': names, 'formats': fmat}
                       )

    # assignes the values
    for name in data.dtype.names:
        datanew[name] = data[name]

    # fill in the new alpha values
    for name in datanew.dtype.names:
        if 'alpha_' in name:
            print(name)
            suffix = name.strip('alpha')
            datanew[name] = 4 * np.pi * data['k' + suffix] / \
                data['wavelength'] / 1e-9 / 100

    # returns the new array
    return datanew

# fnames = ['Si_Schinke_2014','Si_Green_2008','Si_Green_1995','Si_Daub_1995']
fname = 'Si_Vuye_1993'
add_wavelength_oren_ergy_to_file(fname)
add_alpha_to_file(fname)

# plt.semilogy()
# plt.show()
