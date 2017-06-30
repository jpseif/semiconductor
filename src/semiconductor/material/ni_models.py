
import numpy as np
import scipy.constants as Const


def ni_temp(vals, temp, **kargs):
    """
     This form comes from Bludau, Onton, and
     Heinke3 and Macfarlane et a1.31 as cited by Green,3 is
     given by
    """

    if not isinstance(temp, np.ndarray):
        temp = np.asarray([temp])

    if np.all(temp == 0):
        ni = 0
    else:
        ni = vals['a'] * (temp)**vals['power'] * \
            np.exp(- vals['eg'] / temp)

    return ni


def ni_temp_eg(vals, temp,  Eg, *args):
    """
     This form comes from Bludau, Onton, and
     Heinke3 and Macfarlane et a1.31 as cited by Green,3
    """

    if not isinstance(temp, np.ndarray):
        temp = np.asarray([temp])

    if np.all(temp == 0):
        ni = 0
    else:
        print('\n\n', type(vals['A']), vals['A'])
        print(type(vals['power']), vals['power'])
        ni = vals['A'] * temp**vals['power'] * \
            np.exp(- Eg * Const.e / 2. / Const.k / temp)

    return ni
