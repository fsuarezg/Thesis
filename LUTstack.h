//
// Created by Else Groen on 27/03/2017.
//

#ifndef THESIS_LUTSTACK_H
#define THESIS_LUTSTACK_H
#include <arrayfire.h>
#include "angular_spectrum_kernel.h"


class LUTstack {
    public:
    std::vector<af::array> C;

    LUTstack(af::array U, float wlen, float pixelPitch, double depthZ, int N);


    af::array padArray(af::array a, int pad);
    af::array centercrop(af::array X, af::dim4 S);
    bool ispow2(int x);
    const std::vector<af::array> &getC() const;
    void setC(const std::vector<af::array> &C);
    void setLUT(int idx, af::array newLUT);
};


#endif //THESIS_LUTSTACK_H
