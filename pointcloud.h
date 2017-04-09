//
// Created by Else Groen on 27/03/2017.
//

#ifndef THESIS_POINTCLOUD_H
#define THESIS_POINTCLOUD_H


class pointcloud {


    public:
    float pointcloudmatrix[267251][7];
    int WRPindexes[6][2];
    pointcloud(int nWRP);

};


#endif //THESIS_POINTCLOUD_H
