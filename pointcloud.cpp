//
// Created by Else Groen on 27/03/2017.
//

#include "pointcloud.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <tgmath.h>
using namespace std;

//compare function used for sorting in pointcloud
int compare(const void *aa, const void *bb) {
    int *a=(int *)aa;
    int *b=(int *)bb;
    if (a[2]<b[2])
        return -1;
    else if (a[2]==b[2])
        return 0;
    else
        return 1;

}

pointcloud::pointcloud(int nWRP) {
    //READ POINTCLOUD
    string line;
    ifstream myfile("/Users/elsegroen/ClionProjects/Thesis/Venus260K_normNew.txt");
    if (myfile.is_open())
    {
        int row = 0;
        while ( getline (myfile,line) )
        {
            for(int col = 0; col < 7; col++) {
                if(col == 6){
                    pointcloudmatrix[row][col] = 0;
                }
                else {
                    int found = line.find(",");
                    string value = line.substr(0, found);
                    pointcloudmatrix[row][col] = stof(value); //turn string value to float & add to matrix // for double: std::stod
                    line = line.substr(found + 1);
                }
            }
            row++;
        }
        myfile.close();
    }
    else cout << "Unable to open file";

    //SORT POINTS ON Z AXIS
    qsort(pointcloudmatrix, 267251, sizeof(*pointcloudmatrix), compare);

    //ASSIGN WRP TO POINTS
    int WRPcounter = 0;
    int nPerWRP = ceil(267251/nWRP);
    WRPindexes[WRPcounter][0] = 0;
    for(int i=0;i<267251;i++){
        pointcloudmatrix[i][6] = WRPcounter;
        if(i % nPerWRP == 0 && WRPcounter != nWRP && i != 0) {
            WRPindexes[WRPcounter][1] = i;
            WRPcounter++;
            WRPindexes[WRPcounter][0] = i+1;
        }
    }
    WRPindexes[WRPcounter-1][1] = 267251;
    //PRINT CLOUD
    /*
    for(int i=0;i<267251;i++)
    {
        for(int j=0;j<7;j++)
        {
            cout << pointcloudmatrix[i][j] <<"     ";
        }
        cout<<"\n";
    }
     */

}

