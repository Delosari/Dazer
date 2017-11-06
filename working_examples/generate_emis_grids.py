import pyneb as pn
import numpy as np
from scipy import interpolate
from scipy.ndimage.interpolation import map_coordinates
from bisect import bisect_left
from timeit import default_timer as timer


def bilinear_interpolator_axis(x, y, x_range, y_range, data_grid):
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


He1 = pn.RecAtom('He', 1)
Te, ne = 10550.0, 155.0
lines = np.array([4026.0, 4471.0, 5876.0, 6678.0])

min_den, max_den, step_den = 0, 1000, 5
den_interval = np.arange(min_den, max_den, step_den)
den_interval[0] = 1
min_tem, max_tem, step_tem = 7000, 20000, 50
tem_interval = np.arange(min_tem, max_tem, step_tem)

He_emis_grid = np.empty((tem_interval.shape[0], den_interval.shape[0], lines.shape[0]))

for i in range(lines.shape[0]):
    He_emis_grid[:, :, i] = He1.getEmissivity(tem_interval, den_interval, wave=lines[i], product=True)

start = timer()
bilinear_manual_axis = bilinear_interpolator_axis(ne, Te, den_interval, tem_interval, He_emis_grid)
end = timer()
time_new = end - start
print 'new method', time_new

start = timer()
for i in range(lines.shape[0]):
    He_emis_grid[:, :, i] = He1.getEmissivity(tem_interval, den_interval, wave=lines[i], product=True)
end = timer()
time_old = end - start
print 'Old method', time_old

print 'ratio', time_new / time_old

print bilinear_manual_axis
print He1.getEmissivity(Te, ne, wave=lines[0]), He1.getEmissivity(Te, ne, wave=lines[1]), He1.getEmissivity(Te, ne, wave=lines[2]), He1.getEmissivity(Te, ne, wave=lines[3])