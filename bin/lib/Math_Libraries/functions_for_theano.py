import theano.tensor as tt
from theano import function
from pymc3.math import abs_

class BilinearInterpTheano():

    def __init__(self):

        self.indexing_bilinear()

        self.operation_bilinear()

        return

    def neighbour_points(self, x0_idx, y0_idx, x_range, y_range, image_matrix):

        # x1, x2 = x_range[x0_idx:x0_idx + 2]
        # y1, y2 = y_range[y0_idx:y0_idx + 2]
        #
        # z11, z12 = image_matrix[y0_idx][x0_idx:x0_idx + 2]
        # z21, z22 = image_matrix[y0_idx + 1][x0_idx:x0_idx + 2]

        x1, x2 = x_range[x0_idx], x_range[x0_idx+1]
        y1, y2 = y_range[y0_idx], y_range[y0_idx+1]
        z11, z12 = image_matrix[x0_idx,y0_idx], image_matrix[x0_idx,y0_idx+1]
        z21, z22 = image_matrix[x0_idx+1,y0_idx], image_matrix[x0_idx+1,y0_idx+1]

        return x1, x2, y1, y2, z11, z12, z21, z22

    def indexing_bilinear(self):

        # Declare variables
        a = tt.dvector()
        b = tt.dscalar()

        # Build symbolic expression
        out = tt.argmin(tt.abs_(a - b))

        # Compile function
        self.tt_index_array = function([a, b], out)

    def indexing_bilinear_pymc3(self, range_array, in_cord):

        abs_dif = abs_(range_array - in_cord)

        out = tt.argmin(abs_dif)

        return out


    def operation_bilinear(self):

        # Declare variables
        te_1, te_2, ne_1, ne_2, emis_11, emis_12, emis_21, emis_22, x_in, y_in = tt.dscalars(10)

        # Output function
        out_bInterp = (emis_11 * (te_2 - x_in) * (ne_2 - y_in) +
                       emis_21 * (x_in - te_1) * (ne_2 - y_in) +
                       emis_12 * (te_2 - x_in) * (y_in - ne_1) +
                       emis_22 * (x_in - te_1) * (y_in - ne_1)) / ((te_2 - te_1) * (ne_2 - ne_1))

        # Compile function
        self.tt_interp_array = function(inputs=[te_1, te_2, ne_1, ne_2, emis_11, emis_12, emis_21, emis_22, x_in, y_in], outputs=out_bInterp)

    def bilinear_tt_interpolator(self, x, y, x0_idx, y0_idx, x_range, y_range, image_matrix):

        x1, x2, y1, y2, z11, z12, z21, z22 = self.neighbour_points(x0_idx, y0_idx, x_range, y_range, image_matrix)

        interpol_value = self.tt_interp_array(x1, x2, y1, y2, z11, z12, z21, z22, x, y)

        return interpol_value