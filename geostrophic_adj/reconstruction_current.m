clearvars -except ICS ICS_GEO;
close all;
addpath(genpath('/gpfs/f6/bil-coastal-gfdl/scratch/Biao.Zhao/matlab_toolbox/m_map'));
addpath(genpath('/gpfs/f6/bil-coastal-gfdl/scratch/Biao.Zhao/matlab_toolbox/othercolor'));

maskname='../../grid/C3200/ocean_mask.nc';
toponame='../../grid/C3200/topog.nc';
new_ini=strcat(ICS_GEO);disp(new_ini);
old_ini=strcat(ICS);disp(old_ini);
tini=0;


rho0=1025; % Bousinesq background density [kg.m-3]
g=9.8;     % Gravity acceleration [m.s-2]
Omega = 7.292115e-5;
%

ncid = netcdf.open(toponame, 'NC_NOWRITE');
varid = netcdf.inqVarID(ncid, 'depth');
depth = netcdf.getVar(ncid, varid,'double');
netcdf.close(ncid);

ncid = netcdf.open(maskname, 'NC_NOWRITE');
varidmask = netcdf.inqVarID(ncid, 'mask');
maskT = netcdf.getVar(ncid, varidmask,'double');
varidarea = netcdf.inqVarID(ncid, 'areaO');
area = netcdf.getVar(ncid, varidarea,'double');
netcdf.close(ncid);

% Original initial fileds from GLORYS 
ncid = netcdf.open(old_ini, 'NC_NOWRITE');
varidxh = netcdf.inqVarID(ncid, 'xh');
lon_xh = netcdf.getVar(ncid, varidxh,'double');
varidyh = netcdf.inqVarID(ncid, 'yh');
lat_yh = netcdf.getVar(ncid, varidyh,'double');
varidxq = netcdf.inqVarID(ncid, 'xq');
lon_xq = netcdf.getVar(ncid, varidxq,'double');
varidyq = netcdf.inqVarID(ncid, 'yq');
lat_yq = netcdf.getVar(ncid, varidyq,'double');
varidzl = netcdf.inqVarID(ncid, 'zl');
Z = netcdf.getVar(ncid, varidzl,'double');
Z = -1*Z;
varidu = netcdf.inqVarID(ncid, 'u');
U3d_org = netcdf.getVar(ncid, varidu,'double');
varidv = netcdf.inqVarID(ncid, 'v');
V3d_org = netcdf.getVar(ncid, varidv,'double');
varidt = netcdf.inqVarID(ncid, 'temp');
T3d = netcdf.getVar(ncid, varidt,'double');
varids = netcdf.inqVarID(ncid, 'salt');
S3d = netcdf.getVar(ncid, varids,'double');
varidssh = netcdf.inqVarID(ncid, 'ssh');
ssh_org = netcdf.getVar(ncid, varidssh,'double');
netcdf.close(ncid);

NK=length(Z);
NX=length(lon_xh);
NY=length(lat_yh);


U3d_org_T = 0.5*(U3d_org(1:NX, :, :)   + U3d_org(2:NX+1, :, :));  % (NX,NY,NK)
V3d_org_T = 0.5*(V3d_org(:, 1:NY, :)   + V3d_org(:, 2:NY+1, :));  % (NX,NY,NK)

[lat_yh2d,lon_xh2d]=meshgrid(lat_yh,lon_xh);
f = 2*Omega*sin(pi*lat_yh2d/180);
dy_T=3335.9525*ones(size(area,1),size(area,2));
dx_T=area./dy_T;

dy_u = zeros(NX+1,NY);
dy_u(2:NX, :) = 0.5 * ( dy_T(1:NX-1, :) + dy_T(2:NX, :) );
dy_u(1, :)     = dy_T(1, :);
dy_u(NX+1, :) = dy_T(NX, :);

dy_v = zeros(NX,NY+1);
dy_v(:,2:NY)   = 0.5*(dy_T(:,1:NY-1) + dy_T(:,2:NY));
dy_v(:,1)      = dy_T(:,1);
dy_v(:,NY+1)   = dy_T(:,NY);

dx_v = zeros(NX, NY+1);
dx_v(:, 2:NY) = 0.5 * ( dx_T(:, 1:NY-1) + dx_T(:, 2:NY) );
dx_v(:,1)     = dx_T(:,1);
dx_v(:,NY+1) = dx_T(:,NY);

dx_u = zeros(NX+1,NY);
dx_u(2:NX,:)   = 0.5*(dx_T(1:NX-1,:) + dx_T(2:NX,:));
dx_u(1,:)      = dx_T(1,:);
dx_u(NX+1,:)   = dx_T(NX,:);

umask = zeros(NX+1, NY);       % (NX+1, NY)
umask(2:NX,:)   = maskT(1:NX-1,:) .* maskT(2:NX,:);
umask(1,:)       = maskT(1,:);
umask(NX+1, :)   = maskT(NX, :);

vmask = zeros(NX, NY+1);       % (NX, NY+1)
vmask(:, 2:NY)   = maskT(:, 1:NY-1) .* maskT(:, 2:NY);
vmask(:, 1)       = maskT(:, 1);
vmask(:, NY+1)   = maskT(:, NY);


zT = repmat(reshape(Z,[1 1 NK]), [NX NY 1]);
zi = zeros(NX,NY,NK+1);
zi(:,:,2:NK) = 0.5*(zT(:,:,1:NK-1)+zT(:,:,2:NK));
zi(:,:,1)=0; 
zi(:,:,NK+1)=zT(:,:,NK)+(zT(:,:,NK)-zT(:,:,NK-1));

dz_geom = zi(:,:,1:NK) - zi(:,:,2:NK+1);     % >0
	
B3d = -repmat(depth,[1 1 NK]);
mask_full = (zT > B3d);      

dz_T3d = dz_geom .* mask_full;               

rho  = rho_eos(T3d, S3d, zT); 

[u_geo, v_geo] = make_geostrophic_current(ssh_org, rho, dx_T, dy_T, maskT, dz_T3d, f, g, rho0);

 level=1;
 figure;
 u_geo_T = 0.5*(u_geo(1:NX, :, :)   + u_geo(2:NX+1, :, :));  % (NX,NY,NK)
 v_geo_T = 0.5*(v_geo(:, 1:NY, :)   + v_geo(:, 2:NY+1, :));  % (NX,NY,NK)
 speed_geo=sqrt(u_geo_T.^2+v_geo_T.^2);
 set(gcf,'unit','normalized','position',[0.1 0.1 0.65 0.8]);
 set(gca,'unit','normalized','position',[0.1 0.1 0.85 0.75]);
 colormap_RDY=colormap(othercolor('RdYlBu11',20));colormap(flipud(colormap_RDY));
 m_proj('miller','lon',[-99 -35],'lat',[17.5 49.5]);
 %m_pcolor(lon_xh2d,lat_yh2d,ssh_org);shading flat; hold on;
 m_pcolor(lon_xh2d,lat_yh2d,squeeze(speed_geo(:,:,level))); shading flat; hold on;
 m_quiver(lon_xh2d(1:15:end,1:15:end),lat_yh2d(1:15:end,1:15:end),squeeze(u_geo_T(1:15:end,1:15:end,level)),squeeze(v_geo_T(1:15:end,1:15:end,level)),1.5, 'color','k', 'linewidth',1.2, 'MaxHeadSize',0.5);
 clim([0 2]);
 %clim([-1 0.6]);
 m_gshhs_c('color',[0.5 0.5 0.5]);hold on;
 m_grid('tickdir','out','xtick',-105:20:-35,'ytick',17:10:50,'fontsize',18,'linewidth',1.5,'linestyle','none');
 colorbar('fontsize',14);
 name=strcat('Geostrophic ocean currents at:',32,num2str(-1*Z(level)),' m');
 title(name,'fontsize',24);
 hold off;
 grid off;

 figure;
 speed_org=sqrt(U3d_org_T.^2+V3d_org_T.^2).*maskT;
 set(gcf,'unit','normalized','position',[0.1 0.1 0.65 0.8]);
 set(gca,'unit','normalized','position',[0.1 0.1 0.85 0.75]);
 colormap_RDY=colormap(othercolor('RdYlBu11',20));colormap(flipud(colormap_RDY));
 m_proj('miller','lon',[-99 -35],'lat',[17.5 49.5]);
 %m_pcolor(lon_xh2d,lat_yh2d,ssh_org);shading flat; hold on;
 m_pcolor(lon_xh2d,lat_yh2d,squeeze(speed_org(:,:,level))); shading flat; colorbar;hold on;
 m_quiver(lon_xh2d(1:15:end,1:15:end),lat_yh2d(1:15:end,1:15:end),squeeze(U3d_org_T(1:15:end,1:15:end,level).*maskT(1:15:end,1:15:end)),squeeze(V3d_org_T(1:15:end,1:15:end,level).*maskT(1:15:end,1:15:end)),1.5, 'color','k', 'linewidth',1.2, 'MaxHeadSize',0.5);
 clim([0 2]);
 %clim([-1 0.6]);
 m_gshhs_c('color',[0.5 0.5 0.5]);hold on;
 m_grid('tickdir','out','xtick',-105:20:-35,'ytick',17:10:50,'fontsize',18,'linewidth',1.5,'linestyle','none');
 colorbar('fontsize',14);
 name=strcat('original ocean currents at:',32,num2str(-1*Z(level)),' m');
 title(name,'fontsize',24);
 hold off;
 grid off;

 figure;
 speed_diff=squeeze(speed_org(:,:,level)-speed_geo(:,:,level)).*maskT;
 set(gcf,'unit','normalized','position',[0.1 0.1 0.65 0.8]);
 set(gca,'unit','normalized','position',[0.1 0.1 0.85 0.75]);
 colormap_RDY=colormap(othercolor('RdBu11',20));colormap(flipud(colormap_RDY));
 m_proj('miller','lon',[-99 -35],'lat',[17.5 49.5]);
 m_pcolor(lon_xh2d,lat_yh2d,speed_diff); shading flat; colorbar;hold on;
 %clim([0 2]);
 clim([-1 1]);
 m_gshhs_c('color',[0.5 0.5 0.5]);hold on;
 m_grid('tickdir','out','xtick',-105:20:-35,'ytick',17:10:50,'fontsize',18,'linewidth',1.5,'linestyle','none');
 colorbar('fontsize',14);
 name=strcat('original ocean currents at:',32,num2str(-1*Z(level)),' m');
 title(name,'fontsize',24);
 hold off;
 grid off;

 figure
 set(gcf,'unit','normalized','position',[0.01 0 1.2 0.36]);
 set(gca,'unit','normalized','position',[0.08 0.25 0.85 0.7]);
 u_geo_T(dz_T3d==0)=nan;
 section_u_geo=squeeze(u_geo_T(1295,:,:))';
 section_u_geo=smoothdata(section_u_geo,2,'movmean', 5);
 colormap_RDY=colormap(othercolor('Spectral10',20));colormap(flipud(colormap_RDY));
 contourf(squeeze(lat_yh2d(1,:)),Z,section_u_geo,-0.8:0.1:1.2,'linestyle','none');hold on;
 contour(squeeze(lat_yh2d(1,:)),Z,section_u_geo,-0.8:0.2:1.2,'color',[0.4 0.4 0.4],'linewidth',1,'linestyle','-');
 [c1,h1]=contour(squeeze(lat_yh2d(1,:)),Z,section_u_geo,[0.4 0.4],'color','k','linewidth',2);
 clabel(c1,h1);
 [c2,h2]=contour(squeeze(lat_yh2d(1,:)),Z,section_u_geo,[-0.4 -0.4],'color','k','linewidth',2,'linestyle','-.');
 clabel(c2,h2,'color','w');
 ylim([-2000 0]);
 clim([-0.8 1.2]);
 set(gca,'TickDir','out','Box','off','LineWidth',1,'FontSize',16,'Color',[0.2 0.2 0.2]);
 set(gca,'XTick',0:10:50,'XTickLabel',{"0°" "10°N" "20°N" "30°N" "40°N" "50°N"})
 ylabel('Depth (m)','FontSize',24); 
 xlabel('Latitude','FontSize',24); 
 cb = colorbar;
 cb.Label.String = 'zonal current (m s^{-1})';
 cb.Box = 'off';

 figure
 set(gcf,'unit','normalized','position',[0.01 0 1.2 0.36]);
 set(gca,'unit','normalized','position',[0.08 0.25 0.85 0.7]);
 U3d_org_T(dz_T3d==0)=nan;
 section_u_org=squeeze(U3d_org_T(1295,:,:))';
 section_u_org=smoothdata(section_u_org,2,'movmean', 5);
 colormap_RDY=colormap(othercolor('Spectral10',20));colormap(flipud(colormap_RDY));
 contourf(squeeze(lat_yh2d(1,:)),Z,section_u_org,-0.8:0.1:1.2,'linestyle','none');hold on;
 contour(squeeze(lat_yh2d(1,:)),Z,section_u_org,-0.8:0.2:1.2,'color',[0.4 0.4 0.4],'linewidth',1,'linestyle','-');
 [c3,h3]=contour(squeeze(lat_yh2d(1,:)),Z,section_u_org,[0.4 0.4],'color','k','linewidth',2);
 clabel(c3,h3);
 [c4,h4]=contour(squeeze(lat_yh2d(1,:)),Z,section_u_org,[-0.4 -0.4],'color','k','linewidth',2,'linestyle','-.');
 clabel(c4,h4,'color','w');
 ylim([-2000 0]);
 clim([-0.8 1.2]);
 set(gca,'TickDir','out','Box','off','LineWidth',1,'FontSize',16,'Color',[0.2 0.2 0.2]);
 set(gca,'XTick',15:5:50,'XTickLabel',{"15°" "20°N" "25°N" "30°N" "40°N" "45°N" "50°N"})
 ylabel('Depth (m)','FontSize',24); 
 xlabel('Latitude','FontSize',24);
 cb = colorbar;
 cb.Label.String = 'zonal current (m s^{-1})';
 cb.Box = 'off';

 % figure
 % set(gcf,'unit','normalized','position',[0.01 0 1.2 0.36]);
 % set(gca,'unit','normalized','position',[0.08 0.25 0.85 0.7]);
 % v_geo_T(dz_T3d==0)=nan;
 % section_v_geo=squeeze(v_geo_T(:,662,:))';
 % section_v_geo=smoothdata(section_v_geo,2,'movmean', 3);
 % colormap_RDY=colormap(othercolor('Spectral10',30));colormap(flipud(colormap_RDY));
 % contourf(squeeze(lon_xh2d(:,1)),Z,section_v_geo,-1.8:0.1:1.2,'linestyle','none');hold on;
 % contour(squeeze(lon_xh2d(:,1)),Z,section_v_geo,-1.8:0.2:1.2,'color',[0.4 0.4 0.4],'linewidth',1,'linestyle','-');
 % [c1,h1]=contour(squeeze(lon_xh2d(:,1)),Z,section_v_geo,[0.4 0.4],'color','k','linewidth',2);
 % clabel(c1,h1);
 % [c2,h2]=contour(squeeze(lon_xh2d(:,1)),Z,section_v_geo,[-0.4 -0.4],'color','k','linewidth',2,'linestyle','-.');
 % clabel(c2,h2,'color','w');
 % ylim([-2000 0]);
 % clim([-1.8 1.2]);
 % set(gca,'TickDir','out','Box','off','LineWidth',1,'FontSize',16,'Color',[0.2 0.2 0.2]);
 % set(gca,'XTick',-90:10:30,'XTickLabel',{"90W°" "80°W" "70°W" "60°W" "50°W" "40°W" "30°W"})
 % ylabel('Depth (m)','FontSize',24); 
 % xlabel('Longitude','FontSize',24); 
 % cb = colorbar;
 % cb.Label.String = 'meridional current (m s^{-1})';
 % cb.Box = 'off';
 % 
 % figure
 % set(gcf,'unit','normalized','position',[0.01 0 1.2 0.36]);
 % set(gca,'unit','normalized','position',[0.08 0.25 0.85 0.7]);
 % V3d_org_T(dz_T3d==0)=nan;
 % section_v_org=squeeze(V3d_org_T(:,662,:))';
 % section_v_org=smoothdata(section_v_org,2,'movmean', 3);
 % colormap_RDY=colormap(othercolor('Spectral10',30));colormap(flipud(colormap_RDY));
 % contourf(squeeze(lon_xh2d(:,1)),Z,section_v_org,-1.8:0.1:1.2,'linestyle','none');hold on;
 % contour(squeeze(lon_xh2d(:,1)),Z,section_v_org,-1.8:0.2:1.2,'color',[0.4 0.4 0.4],'linewidth',1,'linestyle','-');
 % [c3,h3]=contour(squeeze(lon_xh2d(:,1)),Z,section_v_org,[0.4 0.4],'color','k','linewidth',2);
 % clabel(c3,h3);
 % [c4,h4]=contour(squeeze(lon_xh2d(:,1)),Z,section_v_org,[-0.4 -0.4],'color','k','linewidth',2,'linestyle','-.');
 % clabel(c4,h4,'color','w');
 % ylim([-2000 0]);
 % clim([-1.8 1.2]);
 % set(gca,'TickDir','out','Box','off','LineWidth',1,'FontSize',16,'Color',[0.2 0.2 0.2]);
 % set(gca,'XTick',-90:10:30,'XTickLabel',{"90W°" "80°W" "70°W" "60°W" "50°W" "40°W" "30°W"})
 % ylabel('Depth (m)','FontSize',24); 
 % xlabel('Longitude','FontSize',24);
 % cb = colorbar;
 % cb.Label.String = 'meridional current (m s^{-1})';
 % cb.Box = 'off';


 % figure;
 % colormap_RDY=colormap(othercolor('RdYlBu11',30));colormap(flipud(colormap_RDY));
 % contourf(squeeze(lat_yh2d(1,:)),Z,squeeze(dz_T3d(1860,:,:))','linestyle','none'); shading flat; 
 % colorbar;
 %ylim([-1500 0]);

%============== baratropic adjustment =============================================
U3d_org=u_geo;
V3d_org=v_geo;

H_T = sum(dz_T3d,3);                % (NX, NY)
H_T(H_T<0) = 0;  

dz_u3d = zeros(NX+1, NY, NK);
dz_u3d(2:NX,:,:) = 0.5 * ( dz_T3d(1:NX-1,:,:) + dz_T3d(2:NX,:,:) );
dz_u3d(1,:,:)     = dz_T3d(1,:,:);
dz_u3d(NX+1,:,:) = dz_T3d(NX,:,:);

dz_v3d = zeros(NX, NY+1, NK);
dz_v3d(:,2:NY,:) = 0.5 * ( dz_T3d(:,1:NY-1,:) + dz_T3d(:,2:NY,:) );
dz_v3d(:,1,:)     = dz_T3d(:,1,:);
dz_v3d(:,NY+1,:) = dz_T3d(:,NY,:);

HU = squeeze(sum(U3d_org.* dz_u3d, 3,'omitnan'));   % (NX+1,NY)
HV = squeeze(sum(V3d_org.* dz_v3d, 3,'omitnan'));   % (NX, NY+1)
HU = HU .* umask;
HV = HV .* vmask;

% HU_T = 0.5*(HU(1:NX,:)   + HU(2:NX+1,:));  % (NY,NX)
% HV_T = 0.5*(HV(:,1:NY)   + HV(:,2:NY+1));  % (NY,NX)
% figure;
% m_proj('miller','lon',[-99 -35],'lat',[17.5 49.5]);
% m_contourf(lon_xh2d,lat_yh2d,sqrt(HU_T.^2+HV_T.^2),'linestyle','none');hold on;
% %clim([-2 33]);
% m_gshhs_c('color',[0.5 0.5 0.5]);hold on;
% m_grid('tickdir','out','xtick',-105:20:-35,'ytick',17:10:50,'fontsize',18,'linewidth',1.5,'linestyle','none');
% colorbar('fontsize',14);

Fx_w = HU(1:NX,:)     .* dy_u(1:NX,:);       % West  flux
Fx_e = HU(2:NX+1,:)   .* dy_u(2:NX+1,:);     % East  flux
Fy_s = HV(:,1:NY)     .* dx_v(:,1:NY);       % South flux
Fy_n = HV(:,2:NY+1)   .* dx_v(:,2:NY+1);     % North flux

RHS_T = (Fx_e - Fx_w + Fy_n - Fy_s) ./ area; % (NX,NY)
RHS_T = RHS_T .* maskT;                   

mask_eff = (maskT==1) & (H_T > 0);
[A, lin_id, c0_diag] = build_poisson_T(H_T, mask_eff, area, dx_u, dy_v, dy_u, dx_v);
b = RHS_T(lin_id); 
b = b - mean(b);
[chi_vec,flag,relres,iter] = pcg(A, b, 1e-8, 5000);
fprintf('pcg: flag=%d, relres=%g, iter=%d, ||x||=%g\n', flag, relres, iter, norm(chi_vec));
res = norm(A*chi_vec - b) / max(norm(b),eps);
fprintf('relative residual(Ax-b)=%.3e\n', res)

chi_T = zeros(NX, NY);
chi_T(lin_id) = chi_vec;

He = zeros(NX+1,NY);
He(2:NX,:)   = 0.5*(H_T(1:NX-1,:)+H_T(2:NX,:));
He(1,:)      = H_T(1,:);
He(NX+1,:)   = H_T(NX,:);

Hn = zeros(NX,NY+1);
Hn(:,2:NY)   = 0.5*(H_T(:,1:NY-1)+H_T(:,2:NY));
Hn(:,1)      = H_T(:,1);
Hn(:,NY+1)   = H_T(:,NY);

dchi_u = zeros(NX+1,NY);              % (i+1/2, j)
dchi_u(2:NX,:) = (chi_T(2:NX,:) - chi_T(1:NX-1,:)).* umask(2:NX,:);
dchi_u(1,:)    = 0;   
dchi_u(end,:) = 0;

dchi_v = zeros(NX,NY+1);              % (i, j+1/2)
dchi_v(:,2:NY) = (chi_T(:,2:NY) - chi_T(:,1:NY-1)).* vmask(:,2:NY);
dchi_v(:,1)    = 0;
dchi_v(:,end) = 0;

uc = (dchi_u ./ dx_u) .* umask;
vc = (dchi_v ./ dy_v) .* vmask;

uc3d = repmat(reshape(uc,[NX+1 NY 1]), [1 1 NK]);
vc3d = repmat(reshape(vc,[NX NY+1 1]), [1 1 NK]);
                  
uc_T = 0.5*(uc(1:NX,:)   + uc(2:NX+1,:));  % (NY,NX,NK)
vc_T = 0.5*(vc(:,1:NY)   + vc(:,2:NY+1));  % (NY,NX,NK)

u_new = U3d_org + uc3d;
v_new = V3d_org + vc3d;

u_new_T = 0.5*(u_new(1:NX,:,:) + u_new(2:NX+1,:,:)); % (NY,NX,NK) 
v_new_T = 0.5*(v_new(:,1:NY,:) + v_new(:,2:NY+1,:)); % (NY,NX,NK)

write_uv_to_newini(old_ini, new_ini, u_new, v_new, 1);


 figure;
 speed_new=sqrt(u_new_T.^2+v_new_T.^2);
 set(gcf,'unit','normalized','position',[0.1 0.1 0.65 0.8]);
 set(gca,'unit','normalized','position',[0.1 0.1 0.85 0.75]);
 colormap_RDY=colormap(othercolor('RdYlBu11',20));colormap(flipud(colormap_RDY));
 m_proj('miller','lon',[-99 -35],'lat',[17.5 49.5]);
 %m_contourf(lon_xh2d,lat_yh2d,ssh,'linestyle','none');hold on;
 %m_contourf(lon_xh2d,lat_yh2d,squeeze(speed_geo(:,:,30)),0:0.1:2,'linestyle','none');hold on;
 m_pcolor(lon_xh2d,lat_yh2d,squeeze(speed_new(:,:,level))); shading flat; colorbar;hold on;
 m_quiver(lon_xh2d(1:15:end,1:15:end),lat_yh2d(1:15:end,1:15:end),squeeze(u_new_T(1:15:end,1:15:end,level)),squeeze(v_new_T(1:15:end,1:15:end,level)),1.5, 'color','k', 'linewidth',1.2, 'MaxHeadSize',0.5);
 clim([0 2]);
 %clim([-1 0.6]);
 m_gshhs_c('color',[0.5 0.5 0.5]);hold on;
 m_grid('tickdir','out','xtick',-105:20:-35,'ytick',17:10:50,'fontsize',18,'linewidth',1.5,'linestyle','none');
 colorbar('fontsize',14);
 name=strcat('Sea Surface Height');
 title(name,'fontsize',24);
 hold off;
 grid off;
 %filename = sprintf('SSC_%03d.jpg', k);
 %print(gcf, filename, '-dtiff', '-r50');
 %saveas(gcf, filename, 'tif');

 figure;
 speed_chi=sqrt(uc_T.^2+vc_T.^2);
 set(gcf,'unit','normalized','position',[0.1 0.1 0.65 0.8]);
 set(gca,'unit','normalized','position',[0.1 0.1 0.85 0.75]);
 colormap_RDY=colormap(othercolor('RdYlBu11',20));colormap(flipud(colormap_RDY));
 m_proj('miller','lon',[-99 -35],'lat',[17.5 49.5]);
 %m_contourf(lon_xh2d,lat_yh2d,ssh,'linestyle','none');hold on;
 %m_contourf(lon_xh2d,lat_yh2d,speed_chi,0:0.01:0.1,'linestyle','none');hold on;
 m_pcolor(lon_xh2d,lat_yh2d,speed_chi); shading flat; colorbar;hold on;
 m_quiver(lon_xh2d(1:10:end,1:10:end),lat_yh2d(1:10:end,1:10:end),uc_T(1:10:end,1:10:end),vc_T(1:10:end,1:10:end),5, 'color','k', 'linewidth',1.2, 'MaxHeadSize',0.5);
 %clim([0 0.1]);
 %clim([-1 0.6]);
 m_gshhs_c('color',[0.5 0.5 0.5]);hold on;
 m_grid('tickdir','out','xtick',-105:20:-35,'ytick',17:10:50,'fontsize',18,'linewidth',1.5,'linestyle','none');
 colorbar('fontsize',14);
 name=strcat('Sea Surface Height');
 title(name,'fontsize',24);
 hold off;
 grid off;
 %filename = sprintf('SSC_%03d.jpg', k);
 %print(gcf, filename, '-dtiff', '-r50');
 %saveas(gcf, filename, 'tif');

%===================== Divergence ===========================================
HU_new = HU + (He .* dchi_u ./ dx_u);     % (NX+1,NY)
HV_new = HV + (Hn .* dchi_v ./ dy_v);     % (NX,NY+1)

Fx_w_n = HU_new(1:NX,:)     .* dy_u(1:NX,:);
Fx_e_n = HU_new(2:NX+1,:)   .* dy_u(2:NX+1,:);
Fy_s_n = HV_new(:,1:NY)     .* dx_v(:,1:NY);
Fy_n_n = HV_new(:,2:NY+1)   .* dx_v(:,2:NY+1);

DIV_after = (Fx_e_n - Fx_w_n + Fy_n_n - Fy_s_n) ./ area;
DIV_after = DIV_after .* mask_eff;  


Achi      = A*chi_vec;
DIV_theory          = zeros(NX,NY);
DIV_theory(lin_id)  = b - Achi;  
rb = std(RHS_T(mask_eff), 1);
rt = std(DIV_theory(mask_eff), 1);
ra = std(DIV_after(mask_eff), 1);
fprintf('RMS: before=%.3e  theory(after)=%.3e  reconstructed(after)=%.3e   ratio=%.3e\n', rb, rt, ra,ra/rb);

figure
subplot(2,1,1)
m_proj('miller','lon',[-99 -35],'lat',[17.5 49.5]);
m_pcolor(lon_xh2d,lat_yh2d,RHS_T); shading flat; colorbar;
m_gshhs_c('color','k'); m_grid('tickdir','out','fontsize',14);
clim([-0.01 0.01]);
title('Divergence before','fontsize',14);
subplot(2,1,2)
m_proj('miller','lon',[-99 -35],'lat',[17.5 49.5]);
m_pcolor(lon_xh2d,lat_yh2d,DIV_after); shading flat; colorbar;
m_gshhs_c('color','k'); m_grid('tickdir','out','fontsize',14);
clim([-0.01 0.01]);
title('Divergence after','fontsize',14);


