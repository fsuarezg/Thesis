//
// Created by Else Groen on 27/03/2017.
//

#ifndef THESIS_ANGULAR_SPECTRUM_KERNEL_H
#define THESIS_ANGULAR_SPECTRUM_KERNEL_H
#include "arrayfire.h"
#include <array>
#include "boost/multi_array.hpp"
using namespace std;

class angular_spectrum_kernel {
    af::array angSpecKernel;
public:

    angular_spectrum_kernel(float wlen, float pp, float z, af::dim4 cdim);
    std::array<double, 512> fill_array(int min, int size);
    af::array fftshift(af::array a);
    af::array ifftshift(af::array a);

    const af::array &getAngSpecKernel() const;
    void setAngSpecKernel(const af::array &angSpecKernel);
    af::dim4 getDims();

};


#endif //THESIS_ANGULAR_SPECTRUM_KERNEL_H
