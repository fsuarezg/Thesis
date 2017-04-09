//
// Created by Else Groen on 28/03/2017.
//

#ifndef THESIS_HOLOGRAM_H
#define THESIS_HOLOGRAM_H

#include <arrayfire.h>
#include "pointcloud.h"
#include "LUTstack.h"

class hologram {

public:
    af::array plane;
    const af::array &getPlane() const;

    hologram(const af::array &plane);
    void padPlane(int padsiz);
    void removePadPlane(int padsiz);
    void ang_spec_prop(float wlen, int z, int padsiz, float pp);
    af::array ifftshift(af::array a);
    void occludePoint(int point, pointcloud &PC, int holoRes, float pp, float ocm, af::array &Mask);
    void applyPoint(int point, pointcloud &PC, int holoRes, float pp, int Q, LUTstack &LUT);
};


#endif //THESIS_HOLOGRAM_H
