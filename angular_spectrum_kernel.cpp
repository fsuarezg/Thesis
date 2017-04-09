//
// Created by Else Groen on 27/03/2017.
//

#include "angular_spectrum_kernel.h"

angular_spectrum_kernel::angular_spectrum_kernel(float wlen, float pp, float z, af::dim4 cdim) {
    float array1[cdim[1]];
    float array2[cdim[0]];
    int minvalue = -cdim[1]/2;
    for (int j = 0; j < cdim[1]; ++j) {
        array1[j] = (float) (pow((wlen/pp)*minvalue/cdim[1],2));
        array2[j] = (float) (pow((wlen/pp)*minvalue/cdim[0],2));
        minvalue++;
    }
    af::array X(cdim[1], array1);
    af::array Y(cdim[0], array2);
    X = af::tile(af::transpose(X), cdim[1]);
    Y = af::tile(Y, 1, cdim[0]);
    af::array final_matrix = 2 *  M_PI / wlen * z * af::sqrt(1 - X - Y);
    final_matrix = af::complex(af::constant(0, cdim[1], cdim[1]) ,final_matrix);
    final_matrix = af::exp(final_matrix);
    final_matrix = ifftshift(final_matrix);
    angSpecKernel = final_matrix;
}

array<double, 512> angular_spectrum_kernel::fill_array(int min, int size) {
    std::array<double, 512> array;
    int minvalue = min;
    for(int i = 0; i < size; i++) {
        array[i] = minvalue;
        minvalue += 1;
    }
    return array;
}

af::array angular_spectrum_kernel::fftshift(af::array a) {
    int dim1 = a.dims()[0];
    int dim2 = a.dims()[1];
    a = af::shift(a, (dim1+1)/2, (dim1+1)/2);
    return a;
}

af::array angular_spectrum_kernel::ifftshift(af::array a) {
    int dim1 = a.dims()[0];
    int dim2 = a.dims()[1];
    a = af::shift(a, (dim1+1)/2, (dim1+1)/2);
    return a;

    //NOTES:
    //loc_shiftdim[ii] look up
    //shiftdims
    //shift
    //arr.dims()+1/2
}

const af::array &angular_spectrum_kernel::getAngSpecKernel() const {
    return angSpecKernel;
}

void angular_spectrum_kernel::setAngSpecKernel(const af::array &angSpecKernel) {
    angular_spectrum_kernel::angSpecKernel = angSpecKernel;
}

af::dim4 angular_spectrum_kernel::getDims() {
    return angSpecKernel.dims();
}
