#ifndef CatalystADAPTOR_HEADER
#define CatalystADAPTOR_HEADER

namespace FEAdaptor
{
void Initialize( const char* scripts);

void Finalize();

void CoProcess(
  double *phic, double *grad, double *u, double *v, double *w, double *temp, double *px, double *py, double *pz, int *e2vx, int npnts, int ncvs,int *esec, int *etype, int nsec,int ne2vx_max, double time, unsigned int timeStep, bool lastTimeStep);
}

#endif
