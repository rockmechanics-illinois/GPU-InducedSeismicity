#include "stdio.h"
#include "cuda.h"

#define load(A,nx,ny,nz,Aname)  double *A##_h,*A##_d;\
A##_h = (double*)malloc((nx)*(ny)*(nz)*sizeof(double));\
FILE* A##fid=fopen(Aname, "rb"); fread(A##_h, sizeof(double), (nx)*(ny)*(nz), A##fid); fclose(A##fid);\
cudaMalloc(&A##_d,((nx)*(ny)*(nz))*sizeof(double));\
cudaMemcpy(A##_d,A##_h,((nx)*(ny)*(nz))*sizeof(double),cudaMemcpyHostToDevice)

#define save(A,nx,ny,nz,Aname)  cudaMemcpy(A##_h,A##_d,((nx)*(ny)*(nz))*sizeof(double),cudaMemcpyDeviceToHost);\
FILE* A##fidw=fopen(Aname, "wb"); fwrite(A##_h, sizeof(double), ((nx)*(ny)*(nz)), A##fidw); fclose(A##fidw)




#define Pr(ix,iy,iz)         Mnxny_d[1  - 1 + nMnxny*(ix-1) + nMnxny*(iy-1)*nx + nMnxny*(iz-1)*ny*nx  ]
#define Pf(ix,iy,iz)         Mnxny_d[2  - 1 + nMnxny*(ix-1) + nMnxny*(iy-1)*nx + nMnxny*(iz-1)*ny*nx  ]
#define Txx(ix,iy,iz)        Mnxny_d[3  - 1 + nMnxny*(ix-1) + nMnxny*(iy-1)*nx + nMnxny*(iz-1)*ny*nx  ]
#define Tyy(ix,iy,iz)        Mnxny_d[4  - 1 + nMnxny*(ix-1) + nMnxny*(iy-1)*nx + nMnxny*(iz-1)*ny*nx  ]
#define Tzz(ix,iy,iz)        Mnxny_d[5  - 1 + nMnxny*(ix-1) + nMnxny*(iy-1)*nx + nMnxny*(iz-1)*ny*nx  ]
#define MatIndex(ix,iy,iz)   Mnxny_d[6  - 1 + nMnxny*(ix-1) + nMnxny*(iy-1)*nx + nMnxny*(iz-1)*ny*nx  ]
#define Qs(ix,iy,iz)         Mnxny_d[7  - 1 + nMnxny*(ix-1) + nMnxny*(iy-1)*nx + nMnxny*(iz-1)*ny*nx  ]
#define phi(ix,iy,iz)        Mnxny_d[8  - 1 + nMnxny*(ix-1) + nMnxny*(iy-1)*nx + nMnxny*(iz-1)*ny*nx  ]
#define k_etaf(ix,iy,iz)     Mnxny_d[9  - 1 + nMnxny*(ix-1) + nMnxny*(iy-1)*nx + nMnxny*(iz-1)*ny*nx  ]
#define Txy(ix,iy,iz)      Mintern_d[1  - 1 + nMintern*(ix-1) + nMintern*(iy-1)*(nx-1) + nMintern*(iz-1)*(ny-1)*(nx-1)]
#define Txz(ix,iy,iz)      Mintern_d[2  - 1 + nMintern*(ix-1) + nMintern*(iy-1)*(nx-1) + nMintern*(iz-1)*(ny-1)*(nx-1)]
#define Tyz(ix,iy,iz)      Mintern_d[3  - 1 + nMintern*(ix-1) + nMintern*(iy-1)*(nx-1) + nMintern*(iz-1)*(ny-1)*(nx-1)]
#define  Vx(ix,iy,iz)           Fx_d[1  - 1 + nFx*(ix-1) + nFx*(iy-1)*(nx+1) + nFx*(iz-1)*ny*(nx+1)  ]
#define  qx(ix,iy,iz)           Fx_d[2  - 1 + nFx*(ix-1) + nFx*(iy-1)*(nx+1) + nFx*(iz-1)*ny*(nx+1)  ]
#define  Vy(ix,iy,iz)           Fy_d[1  - 1 + nFy*(ix-1) + nFy*(iy-1)*(nx) + nFy*(iz-1)*(ny+1)*(nx)  ]
#define  qy(ix,iy,iz)           Fy_d[2  - 1 + nFy*(ix-1) + nFy*(iy-1)*(nx) + nFy*(iz-1)*(ny+1)*(nx)  ]
#define  Vz(ix,iy,iz)           Fz_d[1  - 1 + nFz*(ix-1) + nFz*(iy-1)*(nx) + nFz*(iz-1)*(ny)*(nx)  ]
#define  qz(ix,iy,iz)           Fz_d[2  - 1 + nFz*(ix-1) + nFz*(iy-1)*(nx) + nFz*(iz-1)*(ny)*(nx)  ]
#define  MatProp(ix,iy)    MatProp_d[             (ix-1) + (iy-1)*nMatProp1] 

__global__ void compute_V(double* Fx_d, double* Fy_d, double* Fz_d, double* Mnxny_d, double* Mintern_d, double* pa, double* MatProp_d, const  int nx, const  int ny, const int nz, const int nMnxny, const int nMintern, const int nFx, const int nFy, const int nFz, const int nMatProp1)
{
    int ix = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int iy = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int iz = blockIdx.z * blockDim.z + threadIdx.z + 1;
    double dx = pa[0], dy = pa[1], dz = pa[2], dt = pa[3], rho = pa[4], dampX = pa[5], dampY=pa[6],dampZ=pa[7];//, rhofg = pa[8], rhotg = pa[9];
    if (ix >= 2 && ix <= nx && iy >= 2 && iy <= (ny - 1) && iz >= 2 && iz <= (nz-1))
    {
       Vx(ix,iy,iz) = Vx(ix, iy,iz)*dampX + 
                                       dt/rho*(
                                              -( Pr(ix  , iy  , iz  ) -  Pr(ix-1, iy  , iz  ))/dx +
                                               (Txx(ix  , iy  , iz  ) - Txx(ix-1, iy  , iz  ))/dx +
                                               (Txy(ix-1, iy  , iz  ) - Txy(ix-1, iy-1, iz  ))/dy +
                                               (Txz(ix-1, iy  , iz  ) - Txz(ix-1, iy  , iz-1))/dz
                                             );
//         qx(ix,iy,iz) = -MatProp((int)MatIndex(ix,iy,iz),6)*(Pf(ix, iy,iz) - Pf(ix-1,iy,iz))/dx;
        qx(ix,iy,iz) = -0.5*(k_etaf(ix,iy,iz)+k_etaf(ix-1,iy,iz))*(Pf(ix, iy,iz) - Pf(ix-1,iy,iz))/dx;
    }
    if (ix >= 2 && ix <= (nx - 1) && iy >= 2 && iy <= ny && iz >= 2 && iz <= (nz - 1))
    {
        Vy(ix,iy,iz) = Vy(ix, iy,iz)*dampY + 
                                       dt/rho*(
                                              -( Pr(ix  , iy  , iz  ) -  Pr(ix  , iy-1, iz  ))/dy +
                                               (Tyy(ix  , iy  , iz  ) - Tyy(ix  , iy-1, iz  ))/dy +
                                               (Txy(ix  , iy-1, iz  ) - Txy(ix-1, iy-1, iz  ))/dx +
                                               (Tyz(ix  , iy-1, iz  ) - Tyz(ix  , iy-1, iz-1))/dz 
                                              );
//         qy(ix,iy,iz) = -MatProp((int)MatIndex(ix,iy,iz),6)*(Pf(ix,iy,iz) - Pf(ix,iy-1,iz))/dy;
           qy(ix,iy,iz) = -0.5*(k_etaf(ix,iy,iz)+k_etaf(ix,iy-1,iz))*(Pf(ix,iy,iz) - Pf(ix,iy-1,iz))/dy;
    }
    if (ix >= 2 && ix <= (nx - 1) && iy >= 2 && iy <= (ny-1) && iz >= 2 && iz <= (nz))
    {
       Vz(ix,iy,iz) = Vz(ix, iy,iz)*dampZ + 
                                       dt/rho*(
                                              -( Pr(ix  , iy  , iz  ) -  Pr(ix  , iy  , iz-1))/dz +
                                               (Tzz(ix  , iy  , iz  ) - Tzz(ix  , iy  , iz-1))/dz +
                                               (Tyz(ix  , iy  , iz-1) - Tyz(ix  , iy-1, iz-1))/dy +
                                               (Txz(ix  , iy  , iz-1) - Txz(ix-1, iy  , iz-1))/dx +
                                                MatProp(5,7)*2.7*(1-phi(ix,iy,iz))+MatProp(5,7)*phi(ix,iy,iz)   //Heterogeneous 3d porosity from seismic
                                              );
           qz(ix,iy,iz) = -0.5*(k_etaf(ix,iy,iz)+k_etaf(ix,iy,iz-1))*((Pf(ix,iy,iz) - Pf(ix,iy,iz - 1))/dz+MatProp(5,7));
// Equation for case without gravity
//        Vz(ix,iy,iz) = Vz(ix, iy,iz)*dampZ + 
//                                        dt/rho*(
//                                               -( Pr(ix  , iy  , iz  ) -  Pr(ix  , iy  , iz-1))/dz +
//                                                (Tzz(ix  , iy  , iz  ) - Tzz(ix  , iy  , iz-1))/dz +
//                                                (Tyz(ix  , iy  , iz-1) - Tyz(ix  , iy-1, iz-1))/dy +
//                                                (Txz(ix  , iy  , iz-1) - Txz(ix-1, iy  , iz-1))/dx
//                                               );
//            qz(ix,iy,iz) = -0.5*(k_etaf(ix,iy,iz)+k_etaf(ix,iy,iz-1))*((Pf(ix,iy,iz) - Pf(ix,iy,iz - 1))/dz);
    }  
}

__global__ void compute_P(double* Fx_d, double* Fy_d,  double* Fz_d, double* Mnxny_d, double* Mintern_d, double* pa, double* MatProp_d, const  int nx, const  int ny, const int nz, const int nMnxny, const int nMintern, const int nFx, const int nFy, const int nFz, const int nMatProp1, const int indWX, const int indWY, const int indWZ, const int nWX, const int nWY,const int nWZ)
{
    int ix = blockIdx.x * blockDim.x + threadIdx.x + 1; if (ix >= (nx + 1)) return;
    int iy = blockIdx.y * blockDim.y + threadIdx.y + 1; if (iy >= (ny + 1)) return;
    int iz = blockIdx.z * blockDim.z + threadIdx.z + 1; if (iz >= (nz + 1)) return;
    double dx = pa[0], dy = pa[1], dz = pa[2], dt = pa[3], Ks_phys = pa[8];
    double divV   = (Vx(ix+1,iy  ,iz  )-Vx(ix,iy,iz))/dx
                  + (Vy(ix  ,iy+1,iz  )-Vy(ix,iy,iz))/dy
                  + (Vz(ix  ,iy  ,iz+1)-Vz(ix,iy,iz))/dz; 
    double divq   = (qx(ix+1,iy  ,iz  )-qx(ix,iy,iz))/dx + 
                    (qy(ix  ,iy+1,iz  )-qy(ix,iy,iz))/dy +  
                    (qz(ix  ,iy  ,iz+1)-qz(ix,iy,iz))/dz;
    Pr (ix,iy,iz) = Pr(ix,iy,iz) - dt * ( MatProp((int)MatIndex(ix,iy,iz),1)*(                                   divV + MatProp((int)MatIndex(ix,iy,iz),3)                                   *divq)) + (MatProp((int)MatIndex(ix,iy,iz),1)*MatProp((int)MatIndex(ix,iy,iz),3)                                   *Qs(ix,iy,iz)); 
    Pf (ix,iy,iz) = Pf(ix,iy,iz) - dt * ( MatProp((int)MatIndex(ix,iy,iz),1)*(MatProp((int)MatIndex(ix,iy,iz),3)*divV + MatProp((int)MatIndex(ix,iy,iz),3)/MatProp((int)MatIndex(ix,iy,iz),4)*divq)) + (MatProp((int)MatIndex(ix,iy,iz),1)*MatProp((int)MatIndex(ix,iy,iz),3)/MatProp((int)MatIndex(ix,iy,iz),4)*Qs(ix,iy,iz));
//     Pf (ix,iy,iz) = Pf(ix,iy,iz) - dt * ((Pf(ix,iy,iz)-Pf_old(ix,iy,iz))/dt_phys + MatProp((int)MatIndex(ix,iy,iz),1)*(MatProp((int)MatIndex(ix,iy,iz),3)*divV + MatProp((int)MatIndex(ix,iy,iz),3)/MatProp((int)MatIndex(ix,iy,iz),4)*divq)) + (MatProp((int)MatIndex(ix,iy,iz),1)*MatProp((int)MatIndex(ix,iy,iz),3)/MatProp((int)MatIndex(ix,iy,iz),4)*Qs(ix,iy,iz));
    Txx(ix,iy,iz) = Txx(ix,iy,iz) + 2*MatProp((int)MatIndex(ix,iy,iz),2)*dt*((Vx(ix+1, iy  , iz  ) - Vx(ix, iy, iz))/dx-1.0/3.0*divV);
    Tyy(ix,iy,iz) = Tyy(ix,iy,iz) + 2*MatProp((int)MatIndex(ix,iy,iz),2)*dt*((Vy(ix  , iy+1, iz  ) - Vy(ix, iy, iz))/dy-1.0/3.0*divV);
    Tzz(ix,iy,iz) = Tzz(ix,iy,iz) + 2*MatProp((int)MatIndex(ix,iy,iz),2)*dt*((Vz(ix  , iy  , iz+1) - Vz(ix, iy, iz))/dz-1.0/3.0*divV);
    if (ix >= 1 && ix <= (nx - 1) && iy >= 1 && iy <= (ny-1) && iz >= 2 && iz <= (nz-1) )
    {
        Txy(ix,iy,iz)    = Txy(ix,iy,iz) + dt*
            (MatProp((int)MatIndex(ix,iy,iz),2)+MatProp((int)MatIndex(ix+1,iy,iz),2)+MatProp((int)MatIndex(ix,iy+1,iz),2)+MatProp((int)MatIndex(ix+1,iy+1,iz),2))/4.0*
                                                (
                                                  (Vx(ix+1,iy+1,iz  )-Vx(ix+1,iy  ,iz  ))/dy+
                                                  (Vy(ix+1,iy+1,iz  )-Vy(ix  ,iy+1,iz  ))/dx 
                                                );
    }
        if (ix >= 1 && ix <= (nx - 1) && iy >= 2 && iy <= (ny-1) && iz >= 1 && iz <= (nz-1) )
    {
        Txz(ix,iy,iz)    = Txz(ix,iy,iz) + dt*(MatProp((int)MatIndex(ix,iy,iz),2)+MatProp((int)MatIndex(ix+1,iy,iz),2)+MatProp((int)MatIndex(ix,iy,iz+1),2)+MatProp((int)MatIndex(ix+1,iy,iz+1),2))/4.0*
                                                (
                                                 (Vz(ix+1,iy,iz+1)-Vz(ix,iy,iz+1))/dx+
                                                 (Vx(ix+1,iy,iz+1)-Vx(ix+1,iy,iz))/dz   
                                                );
    }
    if (ix >= 2 && ix <= (nx-1) && iy >= 1 && iy <= (ny-1) && iz >= 1 && iz <= (nz-1))
    {
        Tyz(ix,iy,iz)    = Tyz(ix,iy,iz) + dt*(MatProp((int)MatIndex(ix,iy,iz),2)+MatProp((int)MatIndex(ix,iy+1,iz),2)+MatProp((int)MatIndex(ix,iy,iz+1),2)+MatProp((int)MatIndex(ix,iy+1,iz+1),2))/4.0*
                                                (
                                                    (Vy(ix,iy+1,iz+1)-Vy(ix,iy+1,iz  ))/dz+
                                                    (Vz(ix,iy+1,iz+1)-Vz(ix,iy  ,iz+1))/dy
                                                );
    }
    if (ix ==  1) {Pr(ix, iy, iz)  = Pr(ix+1,iy,iz); Pf(ix,iy,iz) = Pf(ix+1,iy,iz); Txx(ix,iy,iz) = Txx(ix+1,iy,iz); Tyy(ix,iy,iz) = Tyy(ix+1,iy,iz); Tzz(ix,iy,iz) = Tzz(ix+1,iy,iz);}
    if (iy ==  1) {Pr(ix, iy, iz)  = Pr(ix,iy+1,iz); Pf(ix,iy,iz) = Pf(ix,iy+1,iz); Txx(ix,iy,iz) = Txx(ix,iy+1,iz); Tyy(ix,iy,iz) = Tyy(ix,iy+1,iz); Tzz(ix,iy,iz) = Tzz(ix,iy+1,iz);}
    if (ix == nx) {Pr(ix, iy, iz)  = Pr(ix-1,iy,iz); Pf(ix,iy,iz) = Pf(ix-1,iy,iz); Txx(ix,iy,iz) = Txx(ix-1,iy,iz); Tyy(ix,iy,iz) = Tyy(ix-1,iy,iz); Tzz(ix,iy,iz) = Tzz(ix-1,iy,iz);}
    if (iy == ny) {Pr(ix, iy, iz)  = Pr(ix,iy-1,iz); Pf(ix,iy,iz) = Pf(ix,iy-1,iz); Txx(ix,iy,iz) = Txx(ix,iy-1,iz); Tyy(ix,iy,iz) = Tyy(ix,iy-1,iz); Tzz(ix,iy,iz) = Tzz(ix,iy-1,iz);}


    if(ix == 1)    {Txz(ix,iy,iz)=Txz(ix+1,iy,iz);}
    if(iy == 1)    {Tyz(ix,iy,iz)=Tyz(ix,iy+1,iz);}
    if(ix == nx-1) {Txz(ix,iy,iz)=Txz(ix-1,iy,iz);}
    if(iy == ny-1) {Tyz(ix,iy,iz)=Tyz(ix,iy-1,iz);}
//     if((ix >= indWX - nWX)&(ix <= indWX + nWX)&(iy >= indWY - nWY)&(iy <= indWY + nWY)&(iz >= indWZ - nWZ)&(iz <= indWZ + nWZ))  {Pf(ix,iy,iz) = 2/Ks_phys;}

    
//     if (ix ==  1) {Txy(ix,iy,iz) = Txy(ix+1,iy,iz); Tyz(ix,iy,iz) = Tyz(ix+1,iy,iz); Txz(ix,iy,iz) = Txz(ix+1,iy,iz);}
//     if (iy ==  1) {Txy(ix,iy,iz) = Txy(ix,iy+1,iz); Tyz(ix,iy,iz) = Tyz(ix,iy+1,iz); Txz(ix,iy,iz) = Txz(ix,iy+1,iz);}
//     if (iz ==  1) {Txy(ix,iy,iz) = Txy(ix,iy,iz+1); Tyz(ix,iy,iz) = Tyz(ix,iy,iz+1); Txz(ix,iy,iz) = Txz(ix,iy,iz+1);}
//     if (ix == nx-1) {Txy(ix,iy,iz) = Txy(ix-1,iy,iz); Tyz(ix,iy,iz) = Tyz(ix-1,iy,iz); Txz(ix,iy,iz) = Txz(ix-1,iy,iz);}
//     if (iy == ny-1) {Txy(ix,iy,iz) = Txy(ix,iy-1,iz); Tyz(ix,iy,iz) = Tyz(ix,iy-1,iz); Txz(ix,iy,iz) = Txz(ix,iy-1,iz);}
//     if (iz ==  1) {Txy(ix,iy,iz) = Txy(ix,iy,iz+1); Tyz(ix,iy,iz) = Tyz(ix,iy,iz+1); Txz(ix,iy,iz) = Txz(ix,iy,iz+1);}   
//     if ((ix==1)&&(iy==1)) {Txy(ix,iy,iz) = Txy(ix+1,iy+1,iz); Tyz(ix,iy,iz) = Tyz(ix+1,iy+1,iz); Txz(ix,iy,iz) = Txz(ix+1,iy+1,iz);}

}


int main()
{
    const int ip=18;
    int ipa[ip],ngridX,ngridY,ngridZ,nt,npars, nMatProp1,nMatProp2, nMnxny, nMintern,nFx,nFy,nFz, indWX, indWY, indWZ, nWX, nWY, nWZ; FILE *fid = fopen("ipa.dat","rb");
    fread(ipa,sizeof(int),ip,fid); fclose(fid);
    ngridX   = ipa[0];
    ngridY   = ipa[1];
    ngridZ   = ipa[2];   
    nt       = ipa[3];
    npars    = ipa[4];
    nMatProp1= ipa[5];
    nMatProp2= ipa[6];
    nMnxny   = ipa[7];
    nMintern = ipa[8];
    nFx      = ipa[9];
    nFy      = ipa[10];
    nFz      = ipa[11];
    indWX    = ipa[12];
    indWY    = ipa[13];
    indWZ    = ipa[14];
    nWX      = ipa[15];
    nWY      = ipa[16];
    nWZ      = ipa[17];
    dim3 grid, block;
    block.x = 8; block.y = 8; block.z = 8;
    grid.x = ngridX; grid.y = ngridY; grid.z = ngridZ;
    const long int nx = block.x * grid.x;
    const long int ny = block.y * grid.y;
    const long int nz = block.z * grid.z;
    cudaSetDevice(0); cudaDeviceReset();
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    load(Mnxny, nMnxny*nx,ny,nz,"Mnxny.dat");
    load(Mintern, nMintern*(nx-1),ny-1,nz-1, "Mintern.dat");
    load(Fx, nFx*(nx+1),     ny   ,     nz   , "Fx.dat");
    load(Fy,      nx   ,nFy*(ny+1),     nz   , "Fy.dat");
    load(Fz,      nx   ,     ny   ,nFz*(nz+1), "Fz.dat");
    load(MatProp, nMatProp1,nMatProp2, 1, "MatProp.dat");
    load(pa,   npars, 1, 1, "pa.dat");
    for (int it = 1; it <= nt; it++)
    {
       compute_V << <grid, block >> > (Fx_d, Fy_d, Fz_d, Mnxny_d, Mintern_d, pa_d, MatProp_d, nx, ny, nz, nMnxny, nMintern, nFx, nFy, nFz, nMatProp1); cudaDeviceSynchronize();
       compute_P << <grid, block >> > (Fx_d, Fy_d, Fz_d, Mnxny_d, Mintern_d, pa_d, MatProp_d, nx, ny, nz, nMnxny, nMintern, nFx, nFy, nFz, nMatProp1, indWX, indWY, indWZ, nWX, nWY, nWZ); cudaDeviceSynchronize();
    }
    save(Mnxny, nMnxny*nx,ny,nz,"Mnxny.dat");
    save(Mintern, nMintern*(nx-1),ny-1,nz-1, "Mintern.dat");
    save(Fx, nFx*(nx+1),     ny   ,     nz   , "Fx.dat");
    save(Fy,      nx   ,nFy*(ny+1),     nz   , "Fy.dat");
    save(Fz,      nx   ,     ny   ,nFz*(nz+1), "Fz.dat");
    cudaDeviceReset();
}