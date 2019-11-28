#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#include "vtk.h"


#include "cfd_adaptor.h"

extern "C" void cfd_catalyst(  int *command, double *phic, double *grad, double *u, double *v, double *w, double *temp, double *px, double *py, double *pz, int *e2vx, int *npnts, int *ncvs,int *esec, int *etype, int *nsec,int *ne2vx_max, int *iout, int *max_it, double *delta_t)//local variables
{

    int np      = *npnts;
    int nc      = *ncvs;
    int ns      = *nsec;
    int ne_max  = *ne2vx_max;
    int itr     = *iout;
    int max_itr = *max_it;
    double  dt  = *delta_t;
    bool timeStep= false;
std::cout<<"iteration "<<itr<<std::endl;
//    if (strcmp(command,"i")==0) {
    if (*command==1) {
        std::cout<<"init starting"<<std::endl;
        FEAdaptor::Initialize("input.py");
        std::cout<<"init completed"<<std::endl;
    }

    //if (strcmp(command,"c")==0) {
    if (*command==2) {
        if ( itr == max_itr - 1) timeStep= true;
        std::cout<<"run starting"<<std::endl;
        FEAdaptor::CoProcess( phic, grad, u, v, w, temp, px, py, pz, e2vx, np, nc, esec, etype, ns, ne_max, dt, itr, timeStep);
        std::cout<<"run completed"<<std::endl;
    }

    //if (strcmp(command,"f")==0) {
    if (*command==3) {
        std::cout<<"finilizing starting"<<std::endl;
        FEAdaptor::Finalize();
        std::cout<<"finilizing completed"<<std::endl;
    }
}

