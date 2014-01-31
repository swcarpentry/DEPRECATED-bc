import numpy as np
from scipy.optimize import fmin

#
# Attempt 1
#

def point_on_line1(x, p1, p2):
    y = p1[1] + (x - p1[0])*(p2[1] - p1[1]) / (p2[0] - p1[0])
    return np.array([x, y])


def dist_from_line1(x, pdata, p1, p2):
    pline = point_on_line1(x, p1, p2)
    return np.sqrt(np.sum((pline - pdata)**2))


def closest_data_to_line1(data, p1, p2):
    dists = np.empty(len(data), dtype=float)
    for i, pdata in enumerate(data):
        x = fmin(dist_from_line1, p1[0], (pdata, p1, p2), disp=False)[0]
        dists[i] = dist_from_line1(x, pdata, p1, p2)
    imin = np.argmin(dists)
    return imin, data[imin]


#
# Attempt 2
#

def dist_from_line2(pdata, p1, p2):
    a = np.sqrt(np.sum((p1 - pdata)**2))
    b = np.sqrt(np.sum((p2 - pdata)**2))
    c = np.sqrt(np.sum((p2 - p1)**2))
    h = a * np.sqrt(1.0 - ((a**2 + c**2 - b**2) / (2.0 * a * c))**2)
    return h

def closest_data_to_line2(data, p1, p2):
    dists = np.empty(len(data), dtype=float)
    for i, pdata in enumerate(data):
        dists[i] = dist_from_line2(pdata, p1, p2)
    imin = np.argmin(dists)
    return imin, data[imin]

#
# Attempt 3
#

def perimeter3(pdata, p1, p2):
    a = np.sqrt(np.sum((p1 - pdata)**2))
    b = np.sqrt(np.sum((p2 - pdata)**2))
    c = np.sqrt(np.sum((p2 - p1)**2))
    return (a + b + c)

def closest_data_to_line3(data, p1, p2):
    peris = np.empty(len(data), dtype=float)
    for i, pdata in enumerate(data):
        peris[i] = perimeter3(pdata, p1, p2)
    imin = np.argmin(peris)
    return imin, data[imin]

#
# Attempt 4
#

def closest_data_to_line4(data, p1, p2):
    return data[np.argmin(np.sqrt(np.sum((p1 - data)**2, axis=1)) + \
                np.sqrt(np.sum((p2 - data)**2, axis=1)))]

