#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;
#include "pointcloud.h"
#include "LUTstack.h"
#include "hologram.h"
#include "angular_spectrum_kernel.h"
#include <arrayfire.h>

int main() {

    int nWRP = 6;
    int nPerWRP = ceil(267251/nWRP);
    int nPoints = 267251;
    int holoRes = 512;
    int Wcor = -75;
    //POINTCLOUD
    pointcloud PC = pointcloud(nWRP);
    //saving min and max depth of object
    float minZ = 100;
    float maxZ = 0;
    for (int i = 0; i < nPoints; i++) {
        if(PC.pointcloudmatrix[i][2] < minZ){
            minZ = PC.pointcloudmatrix[i][2];
        }
        else if(PC.pointcloudmatrix[i][2] > maxZ){
            maxZ = PC.pointcloudmatrix[i][2];
        }
    }
    float ZL[2] = {minZ, maxZ};

    //Make LUT
    int sigma = 3;
    int HRH = 256;
    int HRV = HRH; // size of point profiles
    af::array Gau = af::gaussianKernel(HRH, HRV, sigma, sigma);
    auto max = af::max(af::max(Gau));
    auto min = af::min(af::min(Gau));
    float maxval = max(0).scalar<float>();
    float minval = min(0).scalar<float>();
    Gau -= minval; //Normalisation
    Gau /= (maxval - minval); //Normalisation
    float pp = 8.00e-06;
    float wlen = 640e-9;
    float odepth = fabs(ZL[1] - ZL[0]);
    float dwrp = odepth/nWRP; // depth of each WRP zone OR distance between successive WRPs
    int QL = 63;

    vector<float> levels(QL); // af::array seq = af::seq(-1, 1, 2/(QL-1));
    float x = -1;
    for (int j = 0; j < QL; j++) {
        if (j == QL-1){
            levels[j] = 1 * (dwrp-(dwrp/QL))/2;
        }
        else {levels[j] = x * (dwrp-(dwrp/QL))/2;}
        x = x + ((float) 2/(float)((QL-1)));
    }

    double zl = fabs(fabs(levels[0])-fabs(levels[1]) );
    LUTstack LUT = LUTstack(Gau, wlen, pp, zl, QL);
    af::deviceGC();

    int corrfac = 6;
    vector<float> Wj = levels;
    vector<float> LI = levels;
    for (int k = 0; k < QL; ++k) { //fabs(Wj[k]) * tan(asin(wlen/(2*pp)))
        Wj[k] =  fabs(Wj[k]) * tan(asin(wlen/(2*pp)));
        //cout <<  Wj[k] << endl;
        LI[k] = ceil(Wj[k]/pp) + fmod(ceil(Wj[k]/pp), 2) + corrfac;
        //cout <<  LI[k] << endl;
    }

    //af::array WJmask;
    for (int jj = 0 ; jj < QL; jj++) { //63
        //cut circular apertures
        float cc = 2 * LI[jj]; float cx1 = LI[jj]; float cy1 = LI[jj]; float R1 = LI[jj];

        //Create Mask
        af::array WJmask = af::gaussianKernel(cc, cc, 0, 0);

        float cutoff = WJmask(10).scalar<float>(); // min = WJmask(10)
        WJmask(where(WJmask < cutoff)) = 0;
        WJmask(where(WJmask >= cutoff)) = 1;

        //SurSpe case
        float startcol = HRV/2-LI[jj]+1;
        float endcol = HRV/2+LI[jj];
        float startrow = HRH/2-LI[jj]+1;
        float endrow = HRH/2+LI[jj];
        af::array appliedMask = LUT.getC()[jj](af::seq(startrow,endrow), af::seq(startcol,endcol), 0, 0) * WJmask;
        appliedMask.eval();
        LUT.setLUT(jj, appliedMask); //.cols(startcol, endcol).rows(startrow, endrow) * WJmask);

    }
    af::deviceGC(); //66 af in mem

    af::array plane = af::complex(af::constant(0, holoRes, holoRes));
    hologram H = hologram(plane);
    af::printMemInfo(); //67 af in mem

    int occDim = 5; //size of occlusion mask
    af::array Mask = af::gaussianKernel(occDim, occDim, 100, 100);
    Mask = 1 - (Mask - af::min(af::min(Mask)).scalar<float>())/(af::max(af::max(Mask)).scalar<float>() - af::min(af::min(Mask)).scalar<float>()); //Inverse Gaussian
    Mask.eval();

    std::vector<double> WRPlevels(nWRP);  // z of WRPS
    for(int np = 0; np < nWRP; np++){
        WRPlevels[np] = ZL[1]-(2*np)*odepth/(2*nWRP);
    }
    float ocm = (occDim-1)/2;
    af::deviceGC(); //68 af in mem
    af::printMemInfo();
    for (int cWRP = 0; cWRP < nWRP; ++cWRP) { //for each WRP
        std::cout << "WRP: " << cWRP  << " of " << nWRP << std::endl;
        int lowerbound = PC.WRPindexes[cWRP][0];
        int upperbound = PC.WRPindexes[cWRP][1];

        for (int point = lowerbound; point < upperbound; ++point) { //for each point within WRP
            //occlusion per point
            if(point%1000 == 0){
                cout << point << endl;
                af::deviceGC();
                af::printMemInfo();

            }
            //apply occlusion mask on the WRP plane
            H.occludePoint(point,PC,holoRes,pp,ocm,Mask);
            af::deviceGC();
            float z = PC.pointcloudmatrix[point][2];
            //calculate which LUT should be used for the point
            int Q = QL/2+ (WRPlevels[cWRP] - z) / (dwrp/(QL-1));
            if(Q > QL){ Q = QL; } else if(Q < 1){ Q = 1; }
            //apply the LUT to the WRP plane
            H.applyPoint(point, PC, holoRes, pp, Q, LUT);
            af::deviceGC();
        }
        af_print(H.plane);
        //af_print(H.getPlane().row(251).col(213));
        //saveArray("plane", H.plane, "/Users/elsegroen/ClionProjects/Thesis/plane.txt");
        /*
        if(cWRP < nWRP-1){ // propagate to next WRP
            H.padPlane(holoRes/2); //add padding
            angular_spectrum_kernel angkernel = angular_spectrum_kernel(wlen, pp, zl, H.plane.dims());
            H.plane = af::ifft2(af::fft2(H.plane) * angkernel.getAngSpecKernel());
            H.removePadPlane(holoRes/2);
        }

        else if(cWRP == (nWRP-1)) {  // only for last WRP to holographic plane
            H.ang_spec_prop(wlen, zl, holoRes/2, pp);
            //zrec = HP.getWRPlevels()[HP.getWRPlevels().size()];
            //WRPHolo[chan] = ang_spec_prop(WRPHolo[chan], HS.getWlen()[chan], zrec, padsiz); //WRPHolo, pitch, zwrp, wlen

        }
        */
    }

    return 0;
}

