%% Convert to output type
% Read the full data file
load IBDP_HorizonsNEW.mat

nx           = 8*ngridX;
ny           = 8*ngridY;
nz           = 8*ngridZ;
% The coordinate system is shifted. (x,y)=(0,0) corresponds to
% (x,y) = (1.034745784000000e+05,3.556186964000000e+05)

% Change size of the numerical domain
xMAX = 2200;          yMAX = 3500;
ztop = -1500; zbot = -2150;
% xMAX = max(ARG(:,1)); yMAX = max(ARG(:,2));
dx  = xMAX/8/ngridX;
dy  = yMAX/8/ngridY;
x =  0:dx:xMAX-dx;
y =  0:dy:yMAX-dy;
[X,Y]=ndgrid(x,y);
% Create surface from grid (ARG(:,1),ARG(:,2),ARG(:,3)) on (X,Y)
ZARG = griddata(ARG(:,1),ARG(:,2),ARG(:,3), X,Y, 'natural');
ZPRE = griddata(PRE(:,1),PRE(:,2),PRE(:,3), X,Y, 'natural');
x = 0:dx:X(end,1);
y = 0:dy:Y(1,end);
dz  = (-zbot+ztop)/8/ngridZ
z = zbot:dz:ztop-dz;
% Find indWX, indWY, indWZ corresponding to CCS1 well
[~,indWX] = min((x-1025).^2);
[~,indWY] = min((y-870).^2);
[~,indWZ] = min((z+1929).^2);
% Physical length used for the dimensionless spatial coordinates
Lz_phys  = max(z(:))-min(z(:));
X1 = X; Y1 = Y;
[X,Y,Z] = ndgrid(x,y,z);
%% Initialize 3D heterogeneous porosity from seismic data
if phik_run==1
    load ("C:\Users\Nikita Bondarenko\Dropbox\Nikita work folder\Matlab Work Folder\IBDP Seismicity\3d Seismic volume\PorosityAndTraces.mat")
    phi     = 0.5*nan(nx,ny,nz);
    xx = x_Trace; yy = y_Trace;
    xx(end+1) = -50; yy(end+1)=-50;
    [xqPHI, yqPHI] = ndgrid(x,y);
    for i=1:nz
        [~,indPHI] = min(sqrt((z_Trace(:)-z(i)).^2));
        v   =  squeeze(porosity(indPHI,:));
        v(end+1) = mean(v(:));
        phi(:,:,i) = griddata(xx,yy,v,xqPHI,yqPHI,"natural");
    %     figure(6),clf
    %     contourf(xqPHI,yqPHI,squeeze(phi(:,:,i))), axis image
    %     colormap(flipud(jet)), colorbar, drawnow
    %     squeeze(phi(:,:,i));
    end
        phi(isnan(phi)) = 0.05;
    clear xx yy xqPHI yqPHI indPHI x_Trace y_Trace z_Trace porosity v
end
%% Find layers
    indARG = bsxfun(@lt, ZPRE, reshape(z,1,1,[])); % Find layers above ZPRE (Argenta + Mt Simon)
    indMTS = bsxfun(@lt, ZARG, reshape(z,1,1,[])); % Find layers above ZARG (Mt. Simon)
% Create flags for formation in MatIndex
if formations_run == 0
        MatIndex         =   ones(nx,ny,nz);
else
        MatIndex         = 3*ones(nx,ny,nz); % Initialize PREcambrian layer
    % Finding nearest grid point to faults
    if faults_run == 1
        xF  = Faults(:,1); yF  = Faults(:,2); zF  = Faults(:,3); indF  = Faults(:,4);
        indices = find((zF<-2150)&(zF>-1904)&(indF==3));   
        xF(indices) = []; yF(indices) = []; zF(indices) = []; % exclude any faults deeper than -2100 and shallower than -1904 m
        FaultThicknes = (5*min([dx dy dz]))^2;
        for i=1:length(xF)
              indFX = find(((x-xF(i)).^2)<FaultThicknes);
              indFY = find(((y-yF(i)).^2)<FaultThicknes);
              indFZ = find(((z-zF(i)).^2)<FaultThicknes);
            MatIndex(indFX,indFY,indFZ) = 4;
        end
    end
        MatIndex(indARG) = 2;                % Initialize ARGenta layer
        MatIndex(indMTS) = 1;                % Initialize MT Simon layer
end

% Convert from dimensional spatial coordinates to dimensionless
    z         = (Z - max(Z(:)));
    x         = X./(-min(z(:))); dx = dx/Lz_phys;
    y         = Y./(-min(z(:))); dy = dy/Lz_phys;   
    z         = z./(-min(z(:))); dz = dz/Lz_phys;


X_meshgrid = permute(x,[2,1,3]);
Y_meshgrid = permute(y,[2,1,3]);
Z_meshgrid = permute(z,[2,1,3]);
Ku_meshgrid = permute(k_etaf*ketaf_phys,[2,1,3]);
% s = slice(Lz_phys*X_meshgrid,Lz_phys*Y_meshgrid,Lz_phys*Z_meshgrid-1800, Ku_meshgrid, Lz_phys*[max(X_meshgrid(:))-0.05],Lz_phys*[max(Y_meshgrid(:))-0.5],Lz_phys*[-1]-1800); axis image
s = slice(Lz_phys*X_meshgrid,Lz_phys*Y_meshgrid,Lz_phys*Z_meshgrid-1800, Ku_meshgrid, [1777],[2977],Lz_phys*[-1]-1800); axis image
set(s,'EdgeAlpha', 0.1)
colormap(flipud("jet"))
% hold on
% % scatter3(xF,yF,zF,100,'filled','MarkerFaceAlpha',0.1)
% scatter3(CCS1(:,1),CCS1(:,2),CCS1(:,3), 50, 'bp','filled');
% % scatter3(CCS2(:,1),CCS2(:,2),CCS2(:,3), 50, 'bp','filled');
% % scatter3(VW1(:,1),VW1(:,2),VW1(:,3), 50, 'bp','filled');
% % scatter3(VW2(:,1),VW2(:,2),VW2(:,3), 50, 'bp','filled');
% scatter3(EQ(:,1),EQ(:,2),EQ(:,3), M, 'r','filled','MarkerEdgeColor','k')
% xlim([0 Lz_phys*max(X_meshgrid(:))])
% ylim([0 Lz_phys*max(Y_meshgrid(:))])
% zlim([-2120 -1800])
% surf(X1,Y1,ZARG,'FaceAlpha',0.3,'EdgeAlpha', 0,'FaceColor','r')
% surf(X1,Y1,ZPRE,'FaceAlpha',0.3,'EdgeAlpha', 0)
% ylim([0 3730])
% load CustomColorMap.mat
% % colorbar
% colormap(CustomColormap)

% light("Style","local","Position",[2000 1500 -1950]);
% view(-50,20)
% a1 = 45; a2 = 20;
% lightangle(a1, a2)
% lightangle(a1+90, a2)
% lightangle(a1+270, a2)
% lightangle(a1+360, a2)

if phik_run == 1
    save Geometry.mat x y z MatIndex phi nx ny nz dx dy dz ngridX ngridY ngridZ Lz_phys CCS1 CCS2 VW1 VW2 EQ M ZPRE ZARG indWX indWY indWZ phik_run injection_run formations_run
    clear
else
    save Geometry.mat x y z MatIndex nx ny nz dx dy dz ngridX ngridY ngridZ Lz_phys CCS1 CCS2 VW1 VW2 EQ M ZPRE ZARG indWX indWY indWZ phik_run injection_run formations_run
    clear
end

% figure(2)
% subplot(1,2,1), axis image, contourf(X,Y,ZPRE, 40), axis image, colorbar
% subplot(1,2,2), axis image, contourf(X,Y,ZARG, 40), axis image, colorbar

% figure(2),clf
% contourf(X1,Y1,ZPRE), axis image
% hold on
% scatter3(CCS1(:,1),CCS1(:,2),CCS1(:,3), 50, 'rp','filled');
% scatter3(CCS2(:,1),CCS2(:,2),CCS2(:,3), 50, 'rp','filled');
% scatter3(VW1(:,1),VW1(:,2),VW1(:,3), 50, 'rp','filled');
% scatter3(VW2(:,1),VW2(:,2),VW2(:,3), 50, 'rp','filled');
% drawnow
% 
% 

% scatter3(CCS1(:,1),CCS1(:,2),CCS1(:,3), 50, 'rp','filled');
% scatter3(EQ(:,1),EQ(:,2),EQ(:,3), M, 'k')
% scatter3(CCS1(:,1),CCS1(:,2),CCS1(:,3), 50, 'rp','filled');
% xlim([0 Lz_phys*max(X_meshgrid(:))])
% ylim([0 Lz_phys*max(Y_meshgrid(:))])
% zlim([-2120 -1800])
% view(45,10)







