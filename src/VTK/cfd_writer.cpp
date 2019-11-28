
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <string>

///////////// VTK include files ///////////////////
#include "vtk.h"


extern "C" void c_vtk_writer(  double *phic, double *gradian,double *u, double *v,double *w, double *temperature, double *px, double *py, double *pz, int *e2vx, int *npnts, int *ncvs,int *esec, int *etype, int *nsec,int *ne2vx_max, int *it)//local variables
{


	int x, y, z, i, j,np,ll,ne, ess,eee,type_no;
	bool obsval;
	double u_x, u_y, u_z, d_loc, press;
	double xx, yy, zz;

    std::ostringstream fileNameStream("output/output_no_"+ std::to_string(*it) + ".vtu" );
    std::string fileName = fileNameStream.str();

                         //1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20/
    int element_nvx[20] = {1,1,2,3,3,6,4,8,9, 4,10, 5,14, 6,15,18, 8,20,27,0 };
    int ElementTypeDim[24] = {0,3,0,1,1,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,2,3};


	x = *ncvs/2;
	ne = *ncvs;
	np = *npnts;
	int elm_max= *ne2vx_max;
 	int ns= *nsec ;

    vtkSmartPointer<vtkUnstructuredGrid> VTKGrid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();

    vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

    points->SetDataTypeToDouble();
  // create the points -- slowest in the x and fastest in the z directions
  for ( i = 0; i < np; i++)
  {
        points->InsertNextPoint(px[i] , py[i] , pz[i] );
  }

  VTKGrid->SetPoints(points);

  // create the cells
  for ( int s=0;s<ns; s++){
    int cgns_type=etype[s];
    int vtk_type = cgns_type;
//    std::cout << "Etype is: " << cgns_type << std::endl;
//    std::cout << "ElementTypeDim is: " << ElementTypeDim[cgns_type] << std::endl;
    if (ElementTypeDim[cgns_type]!=2){
    switch(cgns_type) {
        case 10 : type_no=10; break;//tetrahedron
        case 17 : type_no=12; break;//hexahedron
        case 12 : type_no=14; break;//pyramid
        case 14 : type_no=13; break;//prism or wedge
        case 23 : type_no=42; break;//polyhedron
    }
    int elm=element_nvx[cgns_type-1];
//   std::cout << "type_no is: " << type_no << std::endl;
//    std::cout << "Element node number is: " << elm << std::endl;
//    int cellPoints[elm];
    vtkIdType tmp[elm];

    for ( i = esec[s*2]-1; i <esec[s*2+1] ; i++){
            for (j=0; j<elm_max; j++){
                if (j<elm) {
                    tmp[j] = e2vx[i*elm_max + j ]-1;
//                    tmp[j] = cellPoints[j]-1;
                }
            }
        switch(cgns_type) {
        case 10 : VTKGrid->InsertNextCell(VTK_TETRA, elm, tmp);         break;//tetrahedron
        case 17 : VTKGrid->InsertNextCell(VTK_HEXAHEDRON, elm, tmp);    break;//hexahedron
        case 12 : VTKGrid->InsertNextCell(VTK_PYRAMID, elm, tmp);       break;//pyramid
        case 14 : VTKGrid->InsertNextCell(VTK_WEDGE, elm, tmp);         break;//prism or wedge
        case 23 : VTKGrid->InsertNextCell(VTK_POLYHEDRON, elm, tmp);    break;//polyhedron
        }
    }
    }
  }


      // Pressure gradian array
      vtkNew<vtkDoubleArray> grad;
      grad->SetName("GradP");
      grad->SetNumberOfComponents(3);
      grad->SetNumberOfTuples(static_cast<vtkIdType>(ne));
      VTKGrid->GetCellData()->AddArray(grad.GetPointer());
      grad->SetArray(gradian, static_cast<vtkIdType>(ne*3), 1);             ////test

      // velocity Vectors
      vtkNew<vtkDoubleArray> vel;
      vel->SetName("Velocity");
      vel->SetNumberOfComponents(3);
      vel->SetNumberOfTuples(static_cast<vtkIdType>(ne));
      VTKGrid->GetCellData()->AddArray(vel.GetPointer());

    double* velocityData1 = u;
    double* velocityData2 = v;
    double* velocityData3 = w;
    vtkIdType numTuples = vel->GetNumberOfTuples();
    for (vtkIdType i = 0; i < numTuples; i++)
    {
      double values[3] = { velocityData1[i], velocityData2[i], velocityData3[i] };
      vel->SetTuple(i, values);
    }


      // pressure array
      vtkNew<vtkDoubleArray> pressure;
      pressure->SetName("pressure");
      pressure->SetNumberOfComponents(1);
      pressure->SetNumberOfTuples(static_cast<vtkIdType>(ne));
      VTKGrid->GetCellData()->AddArray(pressure.GetPointer());

    double* pressureData = phic;
    pressure->SetArray(pressureData, static_cast<vtkIdType>(ne), 1);

    // temperature array
    vtkNew<vtkDoubleArray> temperatur;
    temperatur->SetName("temperature");
    temperatur->SetNumberOfComponents(1);
    temperatur->SetNumberOfTuples(static_cast<vtkIdType>(ne));
    VTKGrid->GetCellData()->AddArray(temperatur.GetPointer());

    double* temperaturData = temperature;
    temperatur->SetArray(temperaturData, static_cast<vtkIdType>(ne), 1);

  vtkSmartPointer<vtkDirectory> directory = vtkSmartPointer<vtkDirectory>::New();
  int opened = directory->Open("output");

  if(!opened)
  {
    int created = directory->MakeDirectory("output");
  }

  // Write file
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
//  writer->SetFileName("output/test.vtu");
  writer->SetFileName(fileName.c_str());
  writer->SetInputData(VTKGrid);
  writer->SetDataModeToAppended();
  writer->SetEncodeAppendedData(false);
  writer->SetHeaderTypeToUInt64();
  writer->SetIdTypeToInt32();
  writer->SetCompressorTypeToNone();
//  writer->SetDataModeToBinary();
  writer->Write();


}

