% clear
nMat    = 5;               % Number of materials
nProp   = 10;              % Number of properties
MatProp = nan(nMat,nProp);
% nMat == 1 - MT. Simon
% nMat == 2 - ARGenta
% nMat == 3 - intact PREcambrian
% nMat == 4 - Fractured PREcambrian
% nMat == 5 - Water
% nProp == 1 - Undrained bulk modulus/Bulk modulus of water                     [-]
% nProp == 2 - Shear modulus / Bulk modulus of water                            [-]
% nProp == 3 - Skempton's B coefficient                                         [-]
% nProp == 4 - Biot coefficient                                                 [-]
% nProp == 5 - Porosity                                                         [-]
% nProp == 6 - (Permeability/fluidviscosity)/(PermeabilityARG/fluidviscosity)   [-]
% nProp == 7 - rho*g                                                            [-]
% nProp == 8 - cohesion                                                         [MPa]
% nProp == 9 - friction angle                                                   [°]
% nProp == 10- power exponent in porosity-permeability relation                 [-]

% Physical scales
Ks_phys    = 2.1*1e3;     % Bulk modulus of water [MPa]
ketaf_phys = 1.652e-15;   % permeability/water viscosity for ARG

% Water properties
MatProp(5,1) = 1;
MatProp(5,2) = NaN;
MatProp(5,3) = NaN;
MatProp(5,4) = NaN;
MatProp(5,5) = NaN;
MatProp(5,6) = NaN;
MatProp(5,7) = 1e-2/Ks_phys*Lz_phys;
MatProp(5,8) = NaN;
MatProp(5,9) = NaN;
MatProp(5,10)= NaN;
% MT Simon Properties
MatProp(1,1) = 7.8; 
MatProp(1,2) = 8.5;
MatProp(1,3) = 0.33;
MatProp(1,4) = 0.64;
MatProp(1,5) = 0.22;
MatProp(1,6) = 268;
MatProp(1,7) = (MatProp(5,7)*2.7*(1-MatProp(1,5))+MatProp(5,7)*MatProp(1,5));
MatProp(1,8) = 16.5;
MatProp(1,9) = 49.5;
MatProp(1,10)= 3.90;
% ARGenta properties
MatProp(2,1) = 10.2;
MatProp(2,2) = 13.3;
MatProp(2,3) = 0.35;
MatProp(2,4) = 0.64;
MatProp(2,5) = 0.12;
MatProp(2,6) = 1.00;
MatProp(2,7) = (MatProp(5,7)*2.7*(1-MatProp(2,5))+MatProp(5,7)*MatProp(2,5));
MatProp(2,8) = 16.5;
MatProp(2,9) = 45.0;
MatProp(2,10)= 4.20;
% PREcambrian intact properties
MatProp(3,1) = 27.5;
MatProp(3,2) = 21.4;
MatProp(3,3) = 0.30;
MatProp(3,4) = 0.11;
MatProp(3,5) = 0.055;
MatProp(3,6) = 0.0146;
MatProp(3,7) = (MatProp(5,7)*2.7*(1-MatProp(3,5))+MatProp(5,7)*MatProp(3,5));
MatProp(3,8) = 0.00;
MatProp(3,9) = 61.0;
MatProp(3,10)= 15;
% PREcambrian fractured properties
MatProp(4,1) = 20.3;
MatProp(4,2) = 14.6;
MatProp(4,3) = 0.12;
MatProp(4,4) = 0.25;
MatProp(4,5) = 0.09;
MatProp(4,6) = 100.00;
MatProp(4,7) = (MatProp(5,7)*2.7*(1-MatProp(4,5))+MatProp(5,7)*MatProp(4,5));
MatProp(4,8) = 0.00;
MatProp(4,9) = 61.0;
MatProp(4,10)= 28.9;
% Use porosity permeability relation to update the permeability
if formations_run == 1
   ix        = 1:nx;
   iy        = 1:ny;
   iz        = 1:nz; 
    k_etaf = reshape(MatProp(MatIndex(ix,iy,iz),6),nx,ny,nz);
    phi0   = reshape(MatProp(MatIndex(ix,iy,iz),5),nx,ny,nz);
end
if phik_run == 1
ix        = 1:nx;
iy        = 1:ny;
iz        = 1:nz; 
    k_etaf = reshape(MatProp(MatIndex(ix,iy,iz),6),nx,ny,nz).*(phi(ix,iy,iz)./reshape(MatProp(MatIndex(ix,iy,iz), 5),nx,ny,nz)).^reshape(MatProp(MatIndex(ix,iy,iz),10),nx,ny,nz);
%     k_etaf = reshape(MatProp(MatIndex(ix,iy,iz),6),nx,ny,nz);
    k_etaf(k_etaf>400) = 400;
    k_etaf(k_etaf<5e-3)= 5e-3;
    save 'MaterialProperties.mat' MatProp nMat nProp Ks_phys ketaf_phys k_etaf
    clear KdARG KdPRE KdMTS KdPREF Kf
else
    phi  = phi0;
    save 'MaterialProperties.mat' MatProp nMat nProp Ks_phys ketaf_phys
    clear KdARG KdPRE KdMTS KdPREF Kf
end




