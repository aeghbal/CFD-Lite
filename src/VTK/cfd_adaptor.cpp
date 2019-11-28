#include <iostream>
#include "cfd_adaptor.h"
#include "vtk.h"
#include "catalyst.h"

namespace
{
vtkCPProcessor* Processor = NULL;
vtkUnstructuredGrid* VTKGrid;
//vtkUnstructuredGrid* VTKGrid;

void BuildVTKGrid( double *px, double *py, double *pz, int *e2vx, int npnts, int ncvs,int *esec, int *etype, int nsec, int ne2vx_max)
{
	int x, y, z, i, j,np,ll,ne, ess,eee,type_no;
	bool obsval;
	double u_x, u_y, u_z, d_loc, press;
	double xx, yy, zz;

                         //1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20/
    int element_nvx[20] = {1,1,2,3,3,6,4,8,9, 4,10, 5,14, 6,15,18, 8,20,27,0 };
    int ElementTypeDim[24] = {0,3,0,1,1,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,2,3};

    x = ncvs/2;
	ne = ncvs;
	np = npnts;
	int elm_max= ne2vx_max;
 	int ns= nsec ;


//    vtkSmartPointer<vtkUnstructuredGrid> VTKGrid =
//    vtkSmartPointer<vtkUnstructuredGrid>::New();

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
  for ( int s=0;s<nsec; s++){
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
//    std::cout << "type_no is: " << type_no << std::endl;
//    std::cout << "Element node number is: " << elm << std::endl;

    vtkIdType tmp[elm];

    for ( i = esec[s*2]-1; i <esec[s*2+1] ; i++){
            for (j=0; j<elm_max; j++){
                if (j<elm) {
                    tmp[j] = e2vx[i*elm_max + j ]-1;
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
  std::cout<<"VTKGrid built:"<<std::endl;

}

void UpdateVTKAttributes(double *phic, double *gradian, double *u, double *v, double *w, double *temp, int npnts, int ncvs, vtkCPInputDataDescription* idd)
{


  if (idd->IsFieldNeeded("pressure", vtkDataObject::CELL) == true)
  {
//    std::cout<<"pressure values building:"<<std::endl;
    if (VTKGrid->GetCellData()->GetNumberOfArrays() == 0)
    {
      // pressure array
      vtkNew<vtkDoubleArray> pressure;
      pressure->SetName("pressure");
      pressure->SetNumberOfComponents(1);

      pressure->SetNumberOfTuples(static_cast<vtkIdType>(ncvs));      /////////////
      VTKGrid->GetCellData()->AddArray(pressure.GetPointer());
      std::cout<<"pressure array built:"<<std::endl;
    }
    vtkDoubleArray* pressure =
      vtkDoubleArray::SafeDownCast(VTKGrid->GetCellData()->GetArray("pressure"));
//      std::cout<<"pressure safe down cast built:"<<std::endl;
    // The pressure array is a scalar array so we can reuse
    // memory as long as we ordered the points properly.
    double* pressureData = phic;
//    std::cout<<"pressure array addressed:"<<std::endl;
    pressure->SetArray(pressureData, static_cast<vtkIdType>(ncvs), 1);
//    std::cout<<"pressure values Built:"<<std::endl;
   }


  if (idd->IsFieldNeeded("GradP", vtkDataObject::CELL) == true)
  {
//    std::cout<<"Grad values building:"<<std::endl;
    if (VTKGrid->GetCellData()->GetNumberOfArrays() == 1)
    {
      // velocity array
      vtkNew<vtkDoubleArray> grad;
      grad->SetName("GradP");
      grad->SetNumberOfComponents(3);
      grad->SetNumberOfTuples(static_cast<vtkIdType>(ncvs));
      VTKGrid->GetCellData()->AddArray(grad.GetPointer());
      std::cout<<"GradP array built:"<<std::endl;
    }
    vtkDoubleArray* grad =
      vtkDoubleArray::SafeDownCast(VTKGrid->GetCellData()->GetArray("GradP"));
    grad->SetArray(gradian, static_cast<vtkIdType>(ncvs*3), 1);
//    std::cout<<"Grad values built:"<<std::endl;
  }
  if (idd->IsFieldNeeded("Velocity", vtkDataObject::CELL) == true)
  {
//    std::cout<<"Grad values building:"<<std::endl;
    if (VTKGrid->GetCellData()->GetNumberOfArrays() == 2)
    {
      // velocity array
      vtkNew<vtkDoubleArray> vel;
      vel->SetName("Velocity");
      vel->SetNumberOfComponents(3);
      vel->SetNumberOfTuples(static_cast<vtkIdType>(ncvs));
      VTKGrid->GetCellData()->AddArray(vel.GetPointer());
      std::cout<<"Velocity array built:"<<std::endl;
    }
    vtkDoubleArray* vel =
      vtkDoubleArray::SafeDownCast(VTKGrid->GetCellData()->GetArray("Velocity"));
    double* velocityData1 = u;
    double* velocityData2 = v;
    double* velocityData3 = w;
    vtkIdType numTuples = vel->GetNumberOfTuples();
    for (vtkIdType i = 0; i < numTuples; i++)
    {
      double values[3] = { velocityData1[i], velocityData2[i], velocityData3[i] };
      vel->SetTuple(i, values);
    }

//    std::cout<<"Grad values built:"<<std::endl;
  }

  if (idd->IsFieldNeeded("temperature", vtkDataObject::CELL) == true)
  {
//    std::cout<<"temperature values building:"<<std::endl;
    if (VTKGrid->GetCellData()->GetNumberOfArrays() == 3)
    {
      // pressure array
      vtkNew<vtkDoubleArray> temperature;
      temperature->SetName("temperature");
      temperature->SetNumberOfComponents(1);

      temperature->SetNumberOfTuples(static_cast<vtkIdType>(ncvs));      /////////////
      VTKGrid->GetCellData()->AddArray(temperature.GetPointer());
      std::cout<<"pressure array built:"<<std::endl;
    }
    vtkDoubleArray* temperature =
      vtkDoubleArray::SafeDownCast(VTKGrid->GetCellData()->GetArray("temperature"));
//      std::cout<<"temperature safe down cast built:"<<std::endl;
    // The temperature array is a scalar array so we can reuse
    // memory as long as we ordered the points properly.
    double* temperatureData = temp;
//    std::cout<<"pressure array addressed:"<<std::endl;
    temperature->SetArray(temperatureData, static_cast<vtkIdType>(ncvs), 1);
//    std::cout<<"temperature values Built:"<<std::endl;
   }

}

void BuildVTKDataStructures(double *phic, double *grad, double *u, double *v, double *w, double *temp, double *px, double *py, double *pz, int *e2vx, int npnts, int ncvs,int *esec, int *etype, int nsec,int ne2vx_max, vtkCPInputDataDescription* idd)
{
  if (VTKGrid == NULL)
  {
    // The grid structure isn't changing so we only build it
    // the first time it's needed. If we needed the memory
    // we could delete it and rebuild as necessary.
    VTKGrid = vtkUnstructuredGrid::New();
    BuildVTKGrid(px, py, pz, e2vx, npnts, ncvs, esec, etype, nsec, ne2vx_max);
  }
  UpdateVTKAttributes(phic, grad, u, v, w, temp, npnts, ncvs, idd);
}
}

namespace FEAdaptor
{

void Initialize(const char* scripts)
{
  if (Processor == NULL)
  {
    Processor = vtkCPProcessor::New();
    Processor->Initialize();
  }
  else
  {
    Processor->RemoveAllPipelines();
  }
    vtkNew<vtkCPPythonScriptPipeline> pipeline;
    pipeline->Initialize(scripts);
    Processor->AddPipeline(pipeline.GetPointer());
}

void Finalize()
{
  if (Processor)
  {
    Processor->Delete();
    Processor = NULL;
  }
  if (VTKGrid)
  {
    VTKGrid->Delete();
    VTKGrid = NULL;
  }
}

void CoProcess(
  double *phic, double *grad, double *u, double *v, double *w, double *temp, double *px, double *py, double *pz, int *e2vx, int npnts, int ncvs,int *esec, int *etype, int nsec,int ne2vx_max, double time, unsigned int timeStep, bool lastTimeStep)
{
  vtkNew<vtkCPDataDescription> dataDescription;
  dataDescription->AddInput("input");
  dataDescription->SetTimeData(time*timeStep, timeStep);
  if (lastTimeStep == true)
  {
    std::cout<<"last time step"<<std::endl;
    // assume that we want to all the pipelines to execute if it
    // is the last time step.
    dataDescription->ForceOutputOn();
  }
  if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0)
  {
    vtkCPInputDataDescription* idd = dataDescription->GetInputDescriptionByName("input");
    BuildVTKDataStructures(phic, grad, u, v, w, temp, px, py, pz, e2vx, npnts, ncvs, esec, etype, nsec, ne2vx_max, idd);
    idd->SetGrid(VTKGrid);
    Processor->CoProcess(dataDescription.GetPointer());
  }
}
} // end of Catalyst namespace

