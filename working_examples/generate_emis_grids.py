import pyneb as pn
import numpy as np
from scipy import interpolate
from scipy.ndimage.interpolation import map_coordinates
from bisect import bisect_left
from timeit import default_timer as timer

def bilinear_axis_interpolator(x, y, x_range, y_range, data_grid):
    i = bisect_left(x_range, x) - 1
    j = bisect_left(y_range, y) - 1

    x1, x2 = x_range[i:i + 2]
    y1, y2 = y_range[j:j + 2]

    z11, z12 = data_grid[j][i:i + 2]
    z21, z22 = data_grid[j + 1][i:i + 2]

    inter_value = (z11 * (x2 - x) * (y2 - y) +
                   z21 * (x - x1) * (y2 - y) +
                   z12 * (x2 - x) * (y - y1) +
                   z22 * (x - x1) * (y - y1)) / ((x2 - x1) * (y2 - y1))

    return inter_value

def generate_emiss_grid(atom, lines_array, tem_range, den_range):

    grid = np.empty((tem_range.shape[0], den_range.shape[0], lines_array.shape[0]))
    for i in range(lines.shape[0]):
        grid[:, :, i] = atom.getEmissivity(tem_range, den_range, wave=lines[i], product=True)

    return grid

#Declare pyneb objects
H1 = pn.RecAtom('H', 1)
He1 = pn.RecAtom('He', 1)

#Declare target Te, ne and lines
Te, ne = 10550.0, 150.0
Te, ne = 14567.0, 275.0

lines = np.array([4026.0, 4471.0, 5876.0, 6678.0])

#Generate emissivity grid
min_den, max_den, step_den = 0, 1000, 10  # Density grid
den_interval    = np.arange(min_den, max_den, step_den)
den_interval[0] = 1
min_tem, max_tem, step_tem = 7000, 20000, 50  # Temperature grid
tem_interval    = np.arange(min_tem, max_tem, step_tem)
He_emis_grid    = generate_emiss_grid(He1, lines, tem_interval, den_interval)

# ----New approach
start = timer()
bilinear_manual_axis = bilinear_axis_interpolator(ne, Te, den_interval, tem_interval, He_emis_grid)
end = timer()
time_new = end - start
print 'Interpolating together: {} - Time: {} s'.format(bilinear_manual_axis, time_new)

#----Previous approach approach
n_lines = len(lines)
n_range = np.arange(n_lines)
emis_vector = np.empty(len(lines))
start = timer()
for i in n_range:
    emis_vector[i] =  He1.getEmissivity(Te, ne, wave=lines[i], product=False)
end = timer()
time_old = end - start
print 'Pyneb interpolation: {} - Time: {} s'.format(emis_vector, time_old)



# import pyneb as pn
# import numpy as np
# from scipy import interpolate
# from scipy.ndimage.interpolation import map_coordinates
# from bisect import bisect_left
# from timeit import default_timer as timer
#
# def bilinear_interpolator_axis(x, y, x_range, y_range, data_grid):
#     i = bisect_left(x_range, x) - 1
#     j = bisect_left(y_range, y) - 1
#
#     x1, x2 = x_range[i:i + 2]
#     y1, y2 = y_range[j:j + 2]
#
#     z11, z12 = data_grid[j][i:i + 2]
#     z21, z22 = data_grid[j + 1][i:i + 2]
#
#     inter_value = (z11 * (x2 - x) * (y2 - y) +
#                    z21 * (x - x1) * (y2 - y) +
#                    z12 * (x2 - x) * (y - y1) +
#                    z22 * (x - x1) * (y - y1)) / ((x2 - x1) * (y2 - y1))
#
#     return inter_value
#
# H1 = pn.RecAtom('H', 1)
# He1 = pn.RecAtom('He', 1)
# Te, ne = 10552.0, 157.0
# lines_ions = np.array(['S1', 'S2', 'O2', 'O3'])
# hot_ions = np.array(['S1', 'Ar2'])
# lines = np.array([4026.0, 4471.0, 5876.0, 6678.0])
#
# min_den, max_den, step_den = 0, 1000, 5
# den_interval = np.arange(min_den, max_den, step_den)
# den_interval[0] = 1
# min_tem, max_tem, step_tem = 7000, 20000, 10
# tem_interval = np.arange(min_tem, max_tem, step_tem)
#
# He_emis_grid = np.empty((tem_interval.shape[0], den_interval.shape[0], lines.shape[0]))
#
# for i in range(lines.shape[0]):
#     He_emis_grid[:, :, i] = He1.getEmissivity(tem_interval, den_interval, wave=lines[i], product=True)
#
#
# np.save('/home/vital/Desktop/testing_save.npy',He_emis_grid)
#
# He_emis_grid2 = np.load('/home/vital/Desktop/testing_save.npy')
#
# print 'Coso'
# idx_high = np.in1d(lines_ions, hot_ions)
#
# start = timer()
# bilinear_manual_axis = bilinear_interpolator_axis(ne, Te, den_interval, tem_interval, He_emis_grid)
# end = timer()
# time_new = end - start
# print 'new method', time_new
#
#
# bilinear_manual_axis2 = bilinear_interpolator_axis(ne, Te, den_interval, tem_interval, He_emis_grid2)
#
# n_lines = len(lines)
# n_range = np.arange(n_lines)
# emis_vector = np.empty(len(lines))
# start = timer()
# for i in n_range:
#     emis_vector[i] =  He1.getEmissivity(Te, ne, wave=lines[i], product=False)
# end = timer()
# time_old = end - start
# print 'Old method', time_old
# print 'ratio', time_new / time_old
#
# print emis_vector
# print bilinear_manual_axis
# print bilinear_manual_axis2