//
// Created by Else Groen on 27/03/2017.
//

#include "LUTstack.h"

LUTstack::LUTstack(af::array U, float wlen, float pixelPitch, double depthZ, int N) {

    if(ispow2(U.dims(0)) && ispow2(U.dims(1))) {
        if (N < 1 || N % 2 == 1) {
            std::cout << "N must be a positive and odd integer" << "N:" << N << std::endl;
        }
        std::vector<af::array> p(N);
        setC(p);
        af::array paddedU = padArray(U, U.dims(0)/2); //U.dims(1) / 2);


        af::array F = af::fft2(paddedU);
        angular_spectrum_kernel angkernel = angular_spectrum_kernel(wlen, pixelPitch, depthZ, F.dims()); //depthZ
        af::array H = angkernel.getAngSpecKernel();

        float mid;
        if(N%2 == 1){mid= (N/2)+1;} else{mid = N/2;} //ceil(N/2)
        mid = mid-1;
        C[mid] = U;

        af::array G = F;

        for (int i = mid - 1; i > -1; i--) {
            G = af::operator*(G,H);
            C[i] = centercrop(af::ifft2(G), U.dims());
        }

        af::array ones = af::constant(1, H.dims(0), H.dims(1));
        H = af::operator/(ones,H);// 1./H
        for (int i = mid + 1; i < N; i++) {
            F = af::operator*(F, H);
            C[i] = centercrop(af::ifft2(F), U.dims());
        }

    }
    else {
        std::cout << "U must be 2-dimensional, with power of 2 dimensions" << std::endl;
    }
    af::deviceGC();

}

af::array LUTstack::padArray(af::array a, int pad){ //pad a 2 dim array with zeroes, 'pad' is the number of zeroes added to on each side
    af::array side = af::constant(0, a.dims(0), pad); //complex values
    af::array out = af::join(1, side, a, side);
    af::array top = af::constant(0, pad, a.dims(1)+2*pad); // for complex values: af::complex(af::constant(0, pad, a.dims(1)+2*pad));
    out = af::join(0, top , out, top);
    return out; //works
}

af::array LUTstack::centercrop(af::array X, af::dim4 S){ // crops at center of image X with size S
    //af_print(X);
    int xs = floor((X.dims(0) - S[0])/2);
    int ys = floor((X.dims(1) - S[1])/2);
    //std::cout << S[0] << std::endl;
    af::array Y = X.rows(xs+1, xs+S[0]);
    Y = Y.cols(ys+1, ys+S[1]);
    //std::cout << Y.dims(0) << "," << Y.dims(1) << std::endl;
    return Y;
}

bool LUTstack::ispow2(int x){ // checks whether a number is a power of 2
    return (x != 0) && ((x & (x - 1)) == 0);
}

const std::vector<af::array> &LUTstack::getC() const {
    return C;
}

void LUTstack::setC(const std::vector<af::array> &C) {
    LUTstack::C = C;
}

void LUTstack::setLUT(int idx, af::array newLUT) {
    C[idx] = newLUT;
    C[idx].eval();
}

