clear, figure(8), clf, colormap jet, drawnow
delete *.dat *.exe *.exp *.lib
%% Define modules to run
matlab_run     = 0;            % Run CPU code to compare with the GPU
formations_run = 1;            % Initialize MTS, ARG, and PRE layers from seismic data
faults_run     = 1;            % Add faults in the model based on seismic 
injection_run  = 0;            % Initialize non-zero source term at the well area
phik_run       = 1;            % Use seismic to get heterogeneous porosity and permeability
%% Initialize geometry (real)
% Number of grid points nx = 8*ngridX, ny = 8*ngridY, nz = 8*ngridZ
ngridX       = 18;
ngridY       = 28;
ngridZ       = 13;
run GeometryInitFaultsAndPorosity.m
load Geometry.mat
% dx*Lz_phys;
% dy*Lz_phys;
% dz*Lz_phys;
%% Initialize simple geometry for validation
% ng           = 10;
% Lx           = 10;
% Ly           = 10;
% Lz           = 10;
% Lz_phys      = 100;
% ngridX       = ng;
% ngridY       = ng;
% ngridZ       = ng;
% nx           = 8*ngridX;
% ny           = 8*ngridY;
% nz           = 8*ngridZ;
% dx           = Lx/(nx-1);
% dy           = Ly/(ny-1);
% dz           = Lz/(nz-1);
% x            = -Lx/2:dx:Lx/2;
% y            = -Ly/2:dy:Ly/2;
% z            = -Lz/2:dz:Lz/2;
% [x,y,z]      = ndgrid(x,y,z);
% indWX        = fix(nx/2); indWY = fix(ny/2); indWZ = fix(nz/2);
% nWX          = 0; nWY = 0; nWZ = 0;
% MatIndex     = ones(nx,ny,nz);
%% Numerics
CFL          = 1/2;
dampX        = (1.0-5.0/nx);
dampY        = (1.0-5.0/ny);
dampZ        = (1.0-5.0/nz);
nWX = 1; nWY = 1; nWZ = 2; % Number of grid points associated with injection zone
%% Initialization of material properties matrix
run MaterialProperties.m
if phik_run == 1                 % Heterogeneous porosity
    dt          = 1/8*(        min([dx dy dz])).^2./max(           k_etaf(:))./max(            MatProp(:,1));
    dt_real     = 1/8*(Lz_phys*min([dx dy dz])).^2./max(ketaf_phys*k_etaf(:))./max(Ks_phys*1e6*MatProp(:,1));
else                             % Homogeneous porosity
    dt          = 1/8*(        min([dx dy dz])).^2./max(           MatProp(:,6))./max(            MatProp(:,1));
    dt_real     = 1/8*(Lz_phys*min([dx dy dz])).^2./max(ketaf_phys*MatProp(:,6))./max(Ks_phys*1e6*MatProp(:,1));
end
time        = 0;           % Time at the begining of calculations
nt          = 20;
nt2         = 10000;
niter       = 10;
Vsound       = min([dx dy dz])/dt*CFL;
rho0         = 1/2*(max(MatProp(:,1))+4/3*max(MatProp(:,2)))/Vsound^2;
nMnxny       = 9;
nMintern     = 3; 
nFx          = 2;
nFy          = 2;
nFz          = 2;
%               0  1  2   3  4      5   6     7       8       9
pa           = [dx dy dz dt rho0 dampX dampY dampZ Ks_phys];
%                 0     1        2    3      4         5     6      7        8    9    10    11   12      13     14   15    16  17  
ipa          = [ngridX,ngridY,ngridZ, nt,length(pa),nMat, nProp, nMnxny, nMintern nFx, nFy, nFz, indWX, indWY, indWZ, nWX, nWY, nWZ];
fid          = fopen('pa.dat','wb');fwrite(fid,pa(:),'double');fclose(fid);
fid          = fopen('ipa.dat','wb');fwrite(fid,ipa,'int32'); fclose(fid);
%% Initial conditions
% Pf           = 0*exp(-((x-max(x(:))/2)).^2-((y-max(y(:))/2)).^2-((z-min(z(:))/2)).^2);
% Pr           = 0*exp(-((x-max(x(:))/2)).^2-((y-max(y(:))/2)).^2-((z-min(z(:))/2)).^2);
% Pf           = 0*exp(-(x).^2-(y).^2-(z).^2);
% subplot(3,1,1), plot(squeeze(x(:,ny/2,nz/2)),squeeze(Pf(:,ny/2,nz/2)))
% subplot(3,1,2), plot(squeeze(y(nx/2,:,nz/2)),squeeze(Pf(nx/2,:,nz/2)))
% subplot(3,1,3), plot(squeeze(z(nx/2,ny/2,:)),squeeze(Pf(nx/2,ny/2,:)))
% plot(squeeze(x(:,ny/2,nz/2)),squeeze(Pf(:,ny/2,nz/2)))
% drawnow
Pr           = 1*reshape(MatProp(MatIndex,7),nx,ny,nz).*z-10/Ks_phys;
Pf           = 1*-MatProp(5,7)*z+17.5/Ks_phys;


Qs           = zeros(nx,ny,nz);         % Source term for the injection
if injection_run == 1
    Qs(indWX-nWX:indWX+nWX,indWY-nWY:indWY+nWY,indWZ:indWZ+nWZ*2) = 1e-7;
end
Txx          = zeros(nx,ny,nz);
Tyy          = zeros(nx,ny,nz);
Tzz          = zeros(nx,ny,nz);
Txy          = zeros(nx-1,ny-1,nz-1);
Txz          = zeros(nx-1,ny-1,nz-1);
Tyz          = zeros(nx-1,ny-1,nz-1);
Vx           = zeros(nx+1,ny  ,nz  );
Vy           = zeros(nx  ,ny+1,nz  );
Vz           = zeros(nx  ,ny  ,nz+1);
qx           = zeros(nx+1,ny  ,nz  );
qy           = zeros(nx  ,ny+1,nz  );
qz           = zeros(nx  ,ny  ,nz+1);
%% Array allocation
Mnxny        = zeros(nMnxny,nx,ny,nz);
Mintern      = zeros(nMintern,nx-1,ny-1,nz-1);
Fx           = zeros(2,nx+1,ny  ,nz  );
Fy           = zeros(2,nx  ,ny+1,nz  );
Fz           = zeros(2,nx  ,ny  ,nz+1);
%%
Lx = max(x(:))- min(x(:)); Ly = max(y(:))- min(y(:)); Lz = max(z(:))- min(z(:));
[xVx, yVx, zVx] = ndgrid(-(Lx+dx)/2:dx:(Lx+dx)/2, -(Ly)/2   :dy:Ly/2     , -Lz/2     :dz:Lz/2     );
[xVy, yVy, zVy] = ndgrid(-(Lx)/2   :dx:(Lx)/2   , -(Ly+dy)/2:dy:(Ly+dy)/2, -Lz/2     :dz:Lz/2     );
[xVz, yVz, zVz] = ndgrid(-(Lx)/2   :dx:(Lx)/2   , -(Ly)/2   :dy:(Ly)/2   , -(Lz+dz)/2:dz:(Lz+dz)/2);

InitialCond      = -1*1e5/nt*[8 0 0; 0 1.32 0; 0 0 4.5];
Fx(1,:,:,:)      = 1*(-InitialCond(1,1)*xVx - InitialCond(1,2)*yVx - InitialCond(1,3)*zVx);
Fy(1,:,:,:)      = 1*(-InitialCond(2,1)*xVy - InitialCond(2,2)*yVy - InitialCond(2,3)*zVy);
Fz(1,:,:,:)      = 1*(-InitialCond(3,1)*xVz - InitialCond(3,2)*yVz - InitialCond(3,3)*zVz);


%% Mapping
Mnxny(1,:,:,:)   = Pr;
Mnxny(2,:,:,:)   = Pf;
Mnxny(3,:,:,:)   = Txx;
Mnxny(4,:,:,:)   = Tyy;
Mnxny(5,:,:,:)   = Tzz;
Mnxny(6,:,:,:)   = MatIndex;
Mnxny(7,:,:,:)   = Qs;
if formations_run==1
    Mnxny(8,:,:,:)   = phi;
    Mnxny(9,:,:,:)   = k_etaf;
end
Mintern(1,:,:,:) = Txy;
Mintern(2,:,:,:) = Txz;
Mintern(3,:,:,:) = Tyz;
% Fx(1,:,:,:)      = Vx;
Fx(2,:,:,:)      = qx;
% Fy(1,:,:,:)      = Vy;
Fy(2,:,:,:)      = qy;
% Fz(1,:,:,:)      = Vz;
Fz(2,:,:,:)      = qz;
%% Writting dat files for CUDA input
fid = fopen('Mnxny.dat','wb')  ; fwrite(fid, Mnxny(:),'double');  fclose(fid);
fid = fopen('Mintern.dat','wb'); fwrite(fid, Mintern(:),'double');fclose(fid);
fid = fopen('MatProp.dat','wb'); fwrite(fid,MatProp(:),'double'); fclose(fid);
fid = fopen('Fx.dat','wb');      fwrite(fid, Fx(:),'double');     fclose(fid);
fid = fopen('Fy.dat','wb');      fwrite(fid, Fy(:),'double');     fclose(fid);
fid = fopen('Fz.dat','wb');      fwrite(fid, Fz(:),'double');     fclose(fid);
%% Action on GPU
system(['nvcc HM_coupled_05_10_24.cu']);
if ispc
    tic,system('a.exe');GPU_time = toc
else
    tic,system('./a.out');GPU_time = toc
end

fid      = fopen( 'Mnxny.dat');Mnxny      = fread(fid,'double');fclose(fid);
Mnxny    = reshape(Mnxny, [nMnxny,nx,ny,nz]);
Fx(1,:,:,:)      = 0;
Fy(1,:,:,:)      = 0;
Fz(1,:,:,:)      = 0;
Mnxny(2,:,:,:)   = 1*-MatProp(5,7)*z+17.5/Ks_phys;
ipa          = [ngridX,ngridY,ngridZ, nt2,length(pa),nMat, nProp, nMnxny, nMintern nFx, nFy, nFz, indWX, indWY, indWZ, nWX, nWY, nWZ];
fid          = fopen('ipa.dat','wb');fwrite(fid,ipa,'int32'); fclose(fid);
fid = fopen( 'Fx.dat','wb');      fwrite(fid, Fx(:),'double');      fclose(fid);
fid = fopen( 'Fy.dat','wb');      fwrite(fid, Fy(:),'double');      fclose(fid);
fid = fopen( 'Fz.dat','wb');      fwrite(fid, Fz(:),'double');      fclose(fid);
fid = fopen('Mnxny.dat','wb')  ; fwrite(fid, Mnxny(:),'double');  fclose(fid);

% GBs      = nx*ny*nz*nt*16*2*8/GPU_time/1e9;

time_flag = 75;

for iter = 1:niter
time = time + dt_real*nt2;
fid      = fopen( 'Mnxny.dat');Mnxny      = fread(fid,'double');fclose(fid);
Mnxny    = reshape(Mnxny, [nMnxny,nx,ny,nz]);
Mnxny(2,:,:,:)   = 1*-MatProp(5,7)*z+17.5/Ks_phys;
fid = fopen('Mnxny.dat','wb')  ; fwrite(fid, Mnxny(:),'double');  fclose(fid);
    system('a.exe');
                fid      = fopen( 'Mnxny.dat');Mnxny      = fread(fid,'double');fclose(fid);
                fid      = fopen( 'Mintern.dat');Mintern  = fread(fid,'double');fclose(fid);
                fid      = fopen( 'Fx.dat');     Fx       = fread(fid,'double');fclose(fid);
                fid      = fopen( 'Fy.dat');     Fy       = fread(fid,'double');fclose(fid);
                fid      = fopen( 'Fz.dat');     Fz       = fread(fid,'double');fclose(fid);
                Mnxny    = reshape(Mnxny, [nMnxny,nx,ny,nz]);
                Mintern  = reshape(Mintern, [nMintern,nx-1,ny-1,nz-1]);
                Fx       = reshape(Fx, [nFx,nx+1,ny  ,nz  ]);
                Fy       = reshape(Fy, [nFy,nx  ,ny+1,nz  ]);
                Fz       = reshape(Fz, [nFz,nx  ,ny  ,nz+1]);
                PrC = squeeze(Mnxny(1,:,:,:));
                PfC = squeeze(Mnxny(2,:,:,:));
                TxxC = squeeze(Mnxny(3,:,:,:));
                TyyC = squeeze(Mnxny(4,:,:,:));
                TzzC = squeeze(Mnxny(5,:,:,:));
                TxyC = squeeze(Mintern(1,:,:,:));
                TxzC = squeeze(Mintern(2,:,:,:));
                TyzC = squeeze(Mintern(3,:,:,:));
                VxC = squeeze(Fx(1,:,:,:));
                qxC = squeeze(Fx(2,:,:,:));
                VyC = squeeze(Fy(1,:,:,:));
                qyC = squeeze(Fy(2,:,:,:));
                VzC = squeeze(Fz(1,:,:,:));
                qzC = squeeze(Fz(2,:,:,:));
                SxxC = -PrC + TxxC;
                SyyC = -PrC + TyyC;
                SzzC = -PrC + TzzC;
                tiledlayout(3,5)
                nexttile, contourf(Lz_phys*squeeze(x(:,:,indWZ)),Lz_phys*squeeze(y(:,:,indWZ)),squeeze(Ks_phys*PrC(:,:,indWZ))), colorbar, axis image, title('PrC_{xy}');
                nexttile, contourf(Lz_phys*squeeze(x(:,:,indWZ)),Lz_phys*squeeze(y(:,:,indWZ)),squeeze(Ks_phys*(PfC(:,:,indWZ)))), colorbar, axis image, title('PfC_{xy}'); %xlim([900 1200]), ylim([700 1000]);
                nexttile, contourf(Lz_phys*squeeze(x(:,:,indWZ)),Lz_phys*squeeze(y(:,:,indWZ)),squeeze(Ks_phys*SxxC(:,:,indWZ))), colorbar, axis image, title('SxxC_{xy}');                
                nexttile, contourf(Lz_phys*squeeze(x(:,:,indWZ)),Lz_phys*squeeze(y(:,:,indWZ)),squeeze(Ks_phys*SyyC(:,:,indWZ))), colorbar, axis image, title('SyyC_{xy}');  
                nexttile, contourf(Lz_phys*squeeze(x(:,:,indWZ)),Lz_phys*squeeze(y(:,:,indWZ)),squeeze(Ks_phys*SzzC(:,:,indWZ))), colorbar, axis image, title('SzzC_{xy}');       
                nexttile, contourf(Lz_phys*squeeze(x(:,indWY,:)),Lz_phys*squeeze(z(:,indWY,:)),squeeze(Ks_phys*PrC(:,indWY,:))), colorbar, axis image, title('PrC_{xz}');
                nexttile, contourf(Lz_phys*squeeze(x(:,indWY,:)),Lz_phys*squeeze(z(:,indWY,:)),squeeze(Ks_phys*(PfC(:,indWY,:)))), colorbar, axis image, title('PfC_{xz}'); %xlim([900 1200]), ylim([-300 0]);
                nexttile, contourf(Lz_phys*squeeze(x(:,indWY,:)),Lz_phys*squeeze(z(:,indWY,:)),squeeze(Ks_phys*SxxC(:,indWY,:))), colorbar, axis image, title('SxxC_{xz}');
                nexttile, contourf(Lz_phys*squeeze(x(:,indWY,:)),Lz_phys*squeeze(z(:,indWY,:)),squeeze(Ks_phys*SyyC(:,indWY,:))), colorbar, axis image, title('SyyC_{xz}');
                nexttile, contourf(Lz_phys*squeeze(x(:,indWY,:)),Lz_phys*squeeze(z(:,indWY,:)),squeeze(Ks_phys*SzzC(:,indWY,:))), colorbar, axis image, title('SzzC_{xz}');
                nexttile, contourf(Lz_phys*squeeze(y(indWX,:,:)),Lz_phys*squeeze(z(indWX,:,:)),squeeze(Ks_phys*PrC(indWX,:,:))), colorbar, axis image, title('PrC_{yz}');                   
                nexttile, contourf(Lz_phys*squeeze(y(indWX,:,:)),Lz_phys*squeeze(z(indWX,:,:)),squeeze(Ks_phys*PfC(indWX,:,:))), colorbar, axis image, title('PfC_{yz}'); %xlim([700 1000]), ylim([-300 0]);                       
                nexttile, contourf(Lz_phys*squeeze(y(indWX,:,:)),Lz_phys*squeeze(z(indWX,:,:)),squeeze(Ks_phys*SxxC(indWX,:,:))), colorbar, axis image, title('SxxC_{yz}');
                nexttile, contourf(Lz_phys*squeeze(y(indWX,:,:)),Lz_phys*squeeze(z(indWX,:,:)),squeeze(Ks_phys*SyyC(indWX,:,:))), colorbar, axis image, title('SyyC_{yz}');
                nexttile, contourf(Lz_phys*squeeze(y(indWX,:,:)),Lz_phys*squeeze(z(indWX,:,:)),squeeze(Ks_phys*SzzC(indWX,:,:))), colorbar, axis image, title('SzzC_{yz}');           
                drawnow
end


%% Action on CPU
if matlab_run == 1, tic
Vx      = 1*(-InitialCond(1,1)*xVx - InitialCond(1,2)*yVx - InitialCond(1,3)*zVx);
Vy      = 1*(-InitialCond(2,1)*xVy - InitialCond(2,2)*yVy - InitialCond(2,3)*zVy);
Vz      = 1*(-InitialCond(3,1)*xVz - InitialCond(3,2)*yVz - InitialCond(3,3)*zVz);
for iter = 1:niter+1 
    Pf_old = Pf;
    for it = 1:nt
        ix = 2:nx;
        iy = 2:ny-1;
        iz = 2:nz-1;
            Vx(ix,iy,iz) = Vx(ix,iy,iz)*dampX +...
                                          dt/rho0*(...
                                                   -( Pr(ix  ,iy,iz) -  Pr(ix-1,iy  ,iz  ))/dx+...
                                                    (Txx(ix  ,iy,iz) - Txx(ix-1,iy  ,iz  ))/dx+...
                                                    (Txy(ix-1,iy,iz) - Txy(ix-1,iy-1,iz  ))/dy+...
                                                    (Txz(ix-1,iy,iz) - Txz(ix-1,iy  ,iz-1))/dz ...
                                                  );
            qx(ix,iy,iz) = -reshape(MatProp(MatIndex(ix,iy,iz),6),nx-1,ny-2,nz-2).*(Pf(ix,iy,iz) - Pf(ix - 1,iy,iz))/dx;
        ix = 2:nx-1;
        iy = 2:ny;
        iz = 2:nz-1;       
            Vy(ix,iy,iz) = Vy(ix,iy,iz)*dampY +...
                                          dt/rho0*(...
                                                   -( Pr(ix,iy  ,iz) -  Pr(ix  ,iy-1,iz  ))/dy+...
                                                    (Tyy(ix,iy  ,iz) - Tyy(ix  ,iy-1,iz  ))/dy+...
                                                    (Txy(ix,iy-1,iz) - Txy(ix-1,iy-1,iz  ))/dx+...
                                                    (Tyz(ix,iy-1,iz) - Tyz(ix  ,iy-1,iz-1))/dz...
                                                  );
            qy(ix,iy,iz) = -reshape(MatProp(MatIndex(ix,iy,iz),6),nx-2,ny-1,nz-2).*(Pf(ix,iy,iz) - Pf(ix,iy - 1,iz))/dy;
        ix = 2:nx-1;
        iy = 2:ny-1;
        iz = 2:nz;    
            Vz(ix,iy,iz) = Vz(ix,iy,iz)*dampZ +...
                                          dt/rho0*(...
                                                   -( Pr(ix,iy,iz)   -  Pr(ix  ,iy  ,iz-1))/dz+...
                                                    (Tzz(ix,iy,iz)   - Tzz(ix  ,iy  ,iz-1))/dz+...
                                                    (Tyz(ix,iy,iz-1) - Tyz(ix  ,iy-1,iz-1))/dy+...
                                                    (Txz(ix,iy,iz-1) - Txz(ix-1,iy  ,iz-1))/dx ...
                                                    +reshape(MatProp(MatIndex(ix,iy,iz),7), nx-2,ny-2,nz-1)...
                                                  );
            qz(ix,iy,iz) = -reshape(MatProp(MatIndex(ix,iy,iz),6),nx-2,ny-2,nz-1).*((Pf(ix,iy,iz) - Pf(ix,iy,iz - 1))/dz+MatProp(5,7));
%             Vz(ix,iy,iz) = Vz(ix,iy,iz)*dampZ +...
%                                           dt/rho0*(...
%                                                    -( Pr(ix,iy,iz)   -  Pr(ix  ,iy  ,iz-1))/dz+...
%                                                     (Tzz(ix,iy,iz)   - Tzz(ix  ,iy  ,iz-1))/dz+...
%                                                     (Tyz(ix,iy,iz-1) - Tyz(ix  ,iy-1,iz-1))/dy+...
%                                                     (Txz(ix,iy,iz-1) - Txz(ix-1,iy  ,iz-1))/dx ...
%                                                   );
%             qz(ix,iy,iz) = -reshape(MatProp(MatIndex(ix,iy,iz),6),nx-2,ny-2,nz-1).*((Pf(ix,iy,iz) - Pf(ix,iy,iz - 1))/dz);

        ix        = 1:nx;
        iy        = 1:ny;
        iz        = 1:nz; 
        divV          = (Vx(ix+1,iy  ,iz  )-Vx(ix,iy,iz))/dx ...
                      + (Vy(ix  ,iy+1,iz  )-Vy(ix,iy,iz))/dy ...
                      + (Vz(ix  ,iy  ,iz+1)-Vz(ix,iy,iz))/dz;
        divq          = (qx(ix+1,iy  ,iz  )-qx(ix,iy,iz))/dx ...
                      + (qy(ix  ,iy+1,iz  )-qy(ix,iy,iz))/dy ...
                      + (qz(ix  ,iy  ,iz+1)-qz(ix,iy,iz))/dz;          
        Pr(ix,iy,iz)  = Pr(ix,iy,iz) - dt*reshape(MatProp(MatIndex(ix,iy,iz),1),nx,ny,nz).*(1.0*divV + reshape(MatProp(MatIndex(ix,iy,iz),3),nx,ny,nz).*divq);
        Pf(ix,iy,iz)  = Pf(ix,iy,iz) - dt*reshape(MatProp(MatIndex(ix,iy,iz),1),nx,ny,nz).*(reshape(MatProp(MatIndex(ix,iy,iz),3),nx,ny,nz).*divV + reshape(MatProp(MatIndex(ix,iy,iz),3),nx,ny,nz)./reshape(MatProp(MatIndex(ix,iy,iz),4),nx,ny,nz).*divq);
%         Pf(ix,iy,iz)  = Pf(ix,iy,iz) - dt*((Pf(ix,iy,iz)-Pf_old(ix,iy,iz))/dt_phys+reshape(MatProp(MatIndex(ix,iy,iz),1),nx,ny,nz).*(reshape(MatProp(MatIndex(ix,iy,iz),3),nx,ny,nz).*divV + reshape(MatProp(MatIndex(ix,iy,iz),3),nx,ny,nz)./reshape(MatProp(MatIndex(ix,iy,iz),4),nx,ny,nz).*divq));
        Txx(ix,iy,iz) =Txx(ix,iy,iz) + 2*reshape(MatProp(MatIndex(ix,iy,iz),2),nx,ny,nz)*dt.*((Vx(ix+1,iy  ,iz  )-Vx(ix,iy,iz))/dx - divV/3);
        Tyy(ix,iy,iz) =Tyy(ix,iy,iz) + 2*reshape(MatProp(MatIndex(ix,iy,iz),2),nx,ny,nz)*dt.*((Vy(ix  ,iy+1,iz  )-Vy(ix,iy,iz))/dy - divV/3);
        Tzz(ix,iy,iz) =Tzz(ix,iy,iz) + 2*reshape(MatProp(MatIndex(ix,iy,iz),2),nx,ny,nz)*dt.*((Vz(ix  ,iy  ,iz+1)-Vz(ix,iy,iz))/dz - divV/3);
        ix        = 1:nx-1;
        iy        = 1:ny-1;
        iz        = 2:nz-1;     
            Txy(ix,iy,iz)    = Txy(ix,iy,iz) + dt*(reshape(MatProp(MatIndex(ix,iy,iz),2),nx-1,ny-1,nz-2)+reshape(MatProp(MatIndex(ix+1,iy,iz),2),nx-1,ny-1,nz-2)+reshape(MatProp(MatIndex(ix,iy+1,iz),2),nx-1,ny-1,nz-2)+reshape(MatProp(MatIndex(ix+1,iy+1,iz),2),nx-1,ny-1,nz-2))/4.0.*...
                                                    ( ...
                                                      (Vx(ix+1,iy+1,iz)-Vx(ix+1,iy  ,iz))/dy+...
                                                      (Vy(ix+1,iy+1,iz)-Vy(ix  ,iy+1,iz))/dx ...
                                                    );
        ix        = 1:nx-1;
        iy        = 2:ny-1;
        iz        = 1:nz-1;
            Txz(ix,iy,iz)    = Txz(ix,iy,iz) +dt*(reshape(MatProp(MatIndex(ix,iy,iz),2),nx-1,ny-2,nz-1)+reshape(MatProp(MatIndex(ix+1,iy,iz),2),nx-1,ny-2,nz-1)+reshape(MatProp(MatIndex(ix,iy,iz+1),2),nx-1,ny-2,nz-1)+reshape(MatProp(MatIndex(ix+1,iy,iz+1),2),nx-1,ny-2,nz-1))/4.0.*...
                                                    ( ...
                                                     (Vz(ix+1,iy,iz+1)-Vz(ix,iy,iz+1))/dx+...
                                                     (Vx(ix+1,iy,iz+1)-Vx(ix+1,iy,iz))/dz ...
                                                    );
        ix        = 2:nx-1;
        iy        = 1:ny-1;
        iz        = 1:nz-1;
            Tyz(ix,iy,iz)    = Tyz(ix,iy,iz) +dt*(reshape(MatProp(MatIndex(ix,iy,iz),2),nx-2,ny-1,nz-1)+reshape(MatProp(MatIndex(ix,iy+1,iz),2),nx-2,ny-1,nz-1)+reshape(MatProp(MatIndex(ix,iy,iz+1),2),nx-2,ny-1,nz-1)+reshape(MatProp(MatIndex(ix,iy+1,iz+1),2),nx-2,ny-1,nz-1))/4.0.*...
                                                   ( ...
                                                    (Vy(ix,iy+1,iz+1)-Vy(ix,iy+1,iz  ))/dz+...
                                                    (Vz(ix,iy+1,iz+1)-Vz(ix,iy  ,iz+1))/dy ...
                                                   );
        Pr(1,:,:) = Pr(2,:,:); Pr(nx,:,:) = Pr(nx-1,:,:);
        Pr(:,1,:) = Pr(:,2,:); Pr(:,ny,:) = Pr(:,ny-1,:);
        Txx(1,:,:) = Txx(2,:,:); Txx(nx,:,:) = Txx(nx-1,:,:);
        Txx(:,1,:) = Txx(:,2,:); Txx(:,ny,:) = Txx(:,ny-1,:);
        Tyy(1,:,:) = Tyy(2,:,:); Tyy(nx,:,:) = Tyy(nx-1,:,:);
        Tyy(:,1,:) = Tyy(:,2,:); Tyy(:,ny,:) = Tyy(:,ny-1,:);      
        Tzz(1,:,:) = Tzz(2,:,:); Tzz(nx,:,:) = Tzz(nx-1,:,:);
        Tzz(:,1,:) = Tzz(:,2,:); Tzz(:,ny,:) = Tzz(:,ny-1,:);
        Pf(1,:,:) = Pf(2,:,:); Pf(nx,:,:) = Pf(nx-1,:,:);
        Pf(:,1,:) = Pf(:,2,:); Pf(:,ny,:) = Pf(:,ny-1,:);
        Txz(1,:,:) = Txz(2,:,:); Txz(nx-1,:,:) = Txz(nx-2,:,:);
        Tyz(:,1,:) = Tyz(:,2,:); Tyz(:,ny-1,:) = Tyz(:,ny-2,:);
        
    end
    if iter == 1
        Vx(:,:,:) = 0; Vy(:,:,:) = 0; Vz(:,:,:) = 0;
    end
%     CPU_time = toc
%     speedup  = CPU_time/GPU_time
end
    % posprocessing
    fid      = fopen( 'Mnxny.dat');Mnxny      = fread(fid,'double');fclose(fid);
    fid      = fopen( 'Mintern.dat');Mintern  = fread(fid,'double');fclose(fid);
    fid      = fopen( 'Fx.dat');     Fx       = fread(fid,'double');fclose(fid);
    fid      = fopen( 'Fy.dat');     Fy       = fread(fid,'double');fclose(fid);
    fid      = fopen( 'Fz.dat');     Fz       = fread(fid,'double');fclose(fid);
    Mnxny    = reshape(Mnxny, [nMnxny,nx,ny,nz]);
    Mintern  = reshape(Mintern, [nMintern,nx-1,ny-1,nz-1]);
    Fx       = reshape(Fx, [nFx,nx+1,ny  ,nz  ]);
    Fy       = reshape(Fy, [nFy,nx  ,ny+1,nz  ]);
    Fz       = reshape(Fz, [nFz,nx  ,ny  ,nz+1]);
    PrC = squeeze(Mnxny(1,:,:,:));
    PfC = squeeze(Mnxny(2,:,:,:));
    TxxC = squeeze(Mnxny(3,:,:,:));
    TyyC = squeeze(Mnxny(4,:,:,:));
    TzzC = squeeze(Mnxny(5,:,:,:));
    TxyC = squeeze(Mintern(1,:,:,:));
    TxzC = squeeze(Mintern(2,:,:,:));
    TyzC = squeeze(Mintern(3,:,:,:));
figure(1),clf
    subplot(3,3,1), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(Pr(:,:,indWZ))), colorbar, axis image, title('Pr_{xy}');
    subplot(3,3,4), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(Pr(:,indWY,:))), colorbar, axis image, title('Pr_{xz}');
    subplot(3,3,7), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(Pr(indWX,:,:))), colorbar, axis image, title('Pr_{yz}');
    subplot(3,3,2), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(PrC(:,:,indWZ))), colorbar, axis image, title('Pr^{Cuda}_{xy}');
    subplot(3,3,5), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(PrC(:,indWY,:))), colorbar, axis image, title('Pr^{Cuda}_{xz}');
    subplot(3,3,8), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(PrC(indWX,:,:))), colorbar, axis image, title('Pr^{Cuda}_{yz}');
    subplot(3,3,3), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(Pr(:,:,indWZ)-PrC(:,:,indWZ))), colorbar, axis image, title('Matlab-Cuda');
    subplot(3,3,6), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(Pr(:,indWY,:)-PrC(:,indWY,:))), colorbar, axis image, title('Matlab-Cuda');
    subplot(3,3,9), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(Pr(indWX,:,:)-PrC(indWX,:,:))), colorbar, axis image, title('Matlab-Cuda');
figure(2),clf
    subplot(3,3,1), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(Txx(:,:,indWZ))), colorbar, axis image, title('Txx_{xy}');
    subplot(3,3,4), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(Txx(:,indWY,:))), colorbar, axis image, title('Txx_{xz}');
    subplot(3,3,7), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(Txx(indWX,:,:))), colorbar, axis image, title('Txx_{yz}');
    subplot(3,3,2), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(TxxC(:,:,indWZ))), colorbar, axis image, title('Txx^{Cuda}_{xy}');
    subplot(3,3,5), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(TxxC(:,indWY,:))), colorbar, axis image, title('Txx^{Cuda}_{xz}');
    subplot(3,3,8), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(TxxC(indWX,:,:))), colorbar, axis image, title('Txx^{Cuda}_{yz}');
    subplot(3,3,3), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(Txx(:,:,indWZ)-TxxC(:,:,indWZ))), colorbar, axis image, title('Matlab-Cuda');
    subplot(3,3,6), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(Txx(:,indWY,:)-TxxC(:,indWY,:))), colorbar, axis image, title('Matlab-Cuda');
    subplot(3,3,9), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(Txx(indWX,:,:)-TxxC(indWX,:,:))), colorbar, axis image, title('Matlab-Cuda');
figure(3),clf
    subplot(3,3,1), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(Tyy(:,:,indWZ))), colorbar, axis image, title('Tyy_{xy}');
    subplot(3,3,4), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(Tyy(:,indWY,:))), colorbar, axis image, title('Tyy_{xz}');
    subplot(3,3,7), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(Tyy(indWX,:,:))), colorbar, axis image, title('Tyy_{yz}');
    subplot(3,3,2), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(TyyC(:,:,indWZ))), colorbar, axis image, title('Tyy^{Cuda}_{xy}');
    subplot(3,3,5), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(TyyC(:,indWY,:))), colorbar, axis image, title('Tyy^{Cuda}_{xz}');
    subplot(3,3,8), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(TyyC(indWX,:,:))), colorbar, axis image, title('Tyy^{Cuda}_{yz}');
    subplot(3,3,3), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(Tyy(:,:,indWZ)-TyyC(:,:,indWZ))), colorbar, axis image, title('Matlab-Cuda');
    subplot(3,3,6), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(Tyy(:,indWY,:)-TyyC(:,indWY,:))), colorbar, axis image, title('Matlab-Cuda');
    subplot(3,3,9), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(Tyy(indWX,:,:)-TyyC(indWX,:,:))), colorbar, axis image, title('Matlab-Cuda');
figure(4),clf
    subplot(3,3,1), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(Tzz(:,:,indWZ))), colorbar, axis image, title('Tzz_{xy}');
    subplot(3,3,4), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(Tzz(:,indWY,:))), colorbar, axis image, title('Tzz_{xz}');
    subplot(3,3,7), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(Tzz(indWX,:,:))), colorbar, axis image, title('Tzz_{yz}');
    subplot(3,3,2), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(TzzC(:,:,indWZ))), colorbar, axis image, title('Tzz^{Cuda}_{xy}');
    subplot(3,3,5), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(TzzC(:,indWY,:))), colorbar, axis image, title('Tzz^{Cuda}_{xz}');
    subplot(3,3,8), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(TzzC(indWX,:,:))), colorbar, axis image, title('Tzz^{Cuda}_{yz}');
    subplot(3,3,3), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(Tzz(:,:,indWZ)-TzzC(:,:,indWZ))), colorbar, axis image, title('Matlab-Cuda');
    subplot(3,3,6), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(Tzz(:,indWY,:)-TzzC(:,indWY,:))), colorbar, axis image, title('Matlab-Cuda');
    subplot(3,3,9), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(Tzz(indWX,:,:)-TzzC(indWX,:,:))), colorbar, axis image, title('Matlab-Cuda');

%     figure(1), clf
%     subplot(3,3,1), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),Ks_phys*squeeze(Pr(:,:,indWZ))), colorbar, axis image, title('xy');
%     subplot(3,3,4), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),Ks_phys*squeeze(Pr(:,indWY,:))), colorbar, axis image, title('xz');
%     subplot(3,3,7), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),Ks_phys*squeeze(Pr(indWX,:,:))), colorbar, axis image, title('yz');
%     subplot(3,3,2), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),Ks_phys*squeeze(PrC(:,:,indWZ))), colorbar, axis image, title('Cuda:xy');
%     subplot(3,3,5), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),Ks_phys*squeeze(PrC(:,indWY,:))), colorbar, axis image, title('Cuda:xz');
%     subplot(3,3,8), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),Ks_phys*squeeze(PrC(indWX,:,:))), colorbar, axis image, title('Cuda:yz');
%     subplot(3,3,3), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),Ks_phys*squeeze(Pr(:,:,indWZ)-PrC(:,:,indWZ))), colorbar, axis image, title('Cuda:xy');
%     subplot(3,3,6), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),Ks_phys*squeeze(Pr(:,indWY,:)-PrC(:,indWY,:))), colorbar, axis image, title('Cuda:xz');
%     subplot(3,3,9), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),Ks_phys*squeeze(Pr(indWX,:,:)-PrC(indWX,:,:))), colorbar, axis image, title('Cuda:yz');
    figure(5),clf
    subplot(3,3,1), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(Pf(:,:,indWZ))), colorbar, axis image, title('Pf_{xy}');
    subplot(3,3,4), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(Pf(:,indWY,:))), colorbar, axis image, title('Pf_{xz}');
    subplot(3,3,7), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(Pf(indWX,:,:))), colorbar, axis image, title('Pf_{yz}');
    subplot(3,3,2), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(PfC(:,:,indWZ))), colorbar, axis image, title('Pf^{Cuda}_{xy}');
    subplot(3,3,5), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(PfC(:,indWY,:))), colorbar, axis image, title('Pf^{Cuda}_{xz}');
    subplot(3,3,8), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(PfC(indWX,:,:))), colorbar, axis image, title('Pf^{Cuda}_{yz}');
    subplot(3,3,3), contourf(squeeze(x(:,:,indWZ)),squeeze(y(:,:,indWZ)),squeeze(Pf(:,:,indWZ)-PfC(:,:,indWZ))), colorbar, axis image, title('Matlab-Cuda');
    subplot(3,3,6), contourf(squeeze(x(:,indWY,:)),squeeze(z(:,indWY,:)),squeeze(Pf(:,indWY,:)-PfC(:,indWY,:))), colorbar, axis image, title('Matlab-Cuda');
    subplot(3,3,9), contourf(squeeze(y(indWX,:,:)),squeeze(z(indWX,:,:)),squeeze(Pf(indWX,:,:)-PfC(indWX,:,:))), colorbar, axis image, title('Matlab-Cuda');
end

%% Show the stresses computed at the injection interval of CCS1
StateOfStress = 
[SxxC(indWX,indWY,indWZ)*Ks_phys TxyC(indWX,indWY,indWZ)*Ks_phys TxzC(indWX,indWY,indWZ)*Ks_phys;...
 TxyC(indWX,indWY,indWZ)*Ks_phys SyyC(indWX,indWY,indWZ)*Ks_phys TyzC(indWX,indWY,indWZ)*Ks_phys;...
 TxzC(indWX,indWY,indWZ)*Ks_phys TyzC(indWX,indWY,indWZ)*Ks_phys SzzC(indWX,indWY,indWZ)*Ks_phys;]
              

% xlabel(GBs)
delete *.dat *.exe *.exp *.lib