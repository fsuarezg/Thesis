//
// Created by Else Groen on 28/03/2017.
//

#include "hologram.h"

hologram::hologram(const af::array &planeNew) : plane(planeNew) {
    plane.eval();
}

void hologram::padPlane(int padsiz) {
    af::array side = af::complex(af::constant(0, plane.dims(0), padsiz)); //complex values
    af::array out = af::join(1, side, plane, side);
    af::array top = af::complex(af::constant(0, padsiz, plane.dims(1)+2*padsiz)); //complex values
    plane = af::join(0, top, plane, top);
}

void hologram::removePadPlane(int padsiz) {
    plane = plane.rows(1+padsiz, plane.dims(0)-padsiz).cols(1+padsiz, plane.dims(1)-padsiz); //remove padding

}

const af::array &hologram::getPlane() const {
    return plane;
}

void hologram::ang_spec_prop(float wlen, int z, int padsiz, float pp) {

    padPlane(padsiz); //add padding
    plane = fft2(plane);
    int size = plane.dims(1)/2-1;

    float array1[size];
    float array2[size];
    int minvalue = -size/2;
    for (int j = 0; j < size; ++j) {
        array1[j] = (float) (pow((wlen/pp)*minvalue/size,2));
        array2[j] = (float) (pow((wlen/pp)*minvalue/size,2));
        minvalue++;
    }
    af::array X(size, array1);
    af::array Y(size, array2);
    X = af::tile(af::transpose(X), size);
    Y = af::tile(Y, 1, size);
    af::array final_matrix = 2 *  M_PI / wlen * z * af::sqrt(1 - X - Y);
    final_matrix = af::complex(af::constant(0, size, size) ,final_matrix);
    final_matrix = af::exp(final_matrix);
    final_matrix = ifftshift(final_matrix);

    plane = plane * final_matrix;
    plane = af::ifft2(plane);
    removePadPlane(padsiz); //remove padding
}

af::array hologram::ifftshift(af::array a) {
    int dim1 = a.dims()[0];
    int dim2 = a.dims()[1];
    a = af::shift(a, (dim1+1)/2, (dim1+1)/2);
    return a;
}

void hologram::occludePoint(int point, pointcloud &PC, int holoRes, float pp, float ocm, af::array &Mask) {
    float x = PC.pointcloudmatrix[point][0];
    float y = PC.pointcloudmatrix[point][1];
    float z = PC.pointcloudmatrix[point][2];
    float X = round(holoRes/2 + round(x/pp)); //find X & Y position in holo - indJ
    float Y = (holoRes - ceil(holoRes / 2 + round(y / pp)));

    int xfo = (int) std::fmax(1, X - ocm)-1; int xlo = (int) std::fmin(holoRes, X + ocm)-1; //1 moet 0 zijn
    int yfo = (int) std::fmax(1, Y - ocm)-1; int ylo = (int) std::fmin(holoRes, Y + ocm)-1;
    plane(af::seq(yfo,ylo), af::seq(xfo,xlo), 0, 0) *= Mask; // ocMaskN;
    plane.eval();

}

void hologram::applyPoint(int point, pointcloud &PC, int holoRes, float pp, int Q, LUTstack &LUT) {
    float x = PC.pointcloudmatrix[point][0];
    float y = PC.pointcloudmatrix[point][1];
    float z = PC.pointcloudmatrix[point][2];
    float X = round(holoRes/2 + round(x/pp)); //find X & Y position in holo - indJ
    float Y = (holoRes - ceil(holoRes / 2 + round(y / pp)));

    af::array PointField = LUT.getC()[Q-1]; //get correct table
    int SJ = (int) PointField.dims(0); //(16)
    int SI = (int) PointField.dims(1); //(16)
    //crop LUT in case the point is close to the boundary of the hologram
    int xf = (int) std::fmax(0, (X - SI / 2) - 1);
    int xl = (int) std::fmin(holoRes-1, (X + SI / 2 - 1) - 1) ;
    int yf = (int) std::fmax(0, (Y - SJ / 2) - 1);
    int yl = (int) std::fmin(holoRes-1, (Y + SJ / 2 - 1) - 1);
    int mincol = (int) std::fmax(0, (2 - Y + SJ / 2) - 1) ;
    int maxcol = (int) std::fmin(SJ-1, (SJ + holoRes - Y - SJ / 2 + 1) - 1);
    int minrow = (int) std::fmax(0, (2 - X + SI / 2) - 1);
    int maxrow = (int) std::fmin(SI-1, (SI + holoRes - X - SI / 2 + 1) - 1);
    PointField = PointField.cols(mincol, maxcol).rows(minrow, maxrow);
    PointField.eval();
    plane(af::seq(yf,yl), af::seq(xf,xl), 0, 0) += PointField;
    plane.eval();


}
