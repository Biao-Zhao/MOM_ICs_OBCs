function [u_geo, v_geo] = make_geostrophic_current(ssh_T, rho, dx_T, dy_T, maskT, dz_T3d, f2d, g, rho0)
%--------------------------------------------------------------------------
% Geostrophic velocity from SSH and density on a C-grid (z-levels).
% All derivatives & vertical integration are done on T cells first.
% Finally, interpolate once from T to U/V faces, enforcing 3-D wet masks.
%
% INPUTS
%   ssh_T    : sea surface height on T (NX,NY)
%   rho      : in-situ density on T (NX,NY,NK) at layer centers
%   dx_T     : grid spacing in x at T (NX,NY)   [meters]
%   dy_T     : grid spacing in y at T (NX,NY)   [meters]
%   umask    : 2-D U-face land/sea mask (NX+1,NY), 1=ocean, 0=land
%   vmask    : 2-D V-face land/sea mask (NX,NY+1)
%   maskT    : 2-D T-cell land/sea mask (NX,NY)
%   dz_T3d   : layer thickness at T (NX,NY,NK), >0 means wet at that level
%   f2d      : Coriolis parameter at T (NX,NY)
%   g        : gravity [m/s^2]
%   rho0     : reference density [kg/m^3]
%
% OUTPUTS
%   u_geo    : geostrophic U on U faces (NX+1,NY,NK)
%   v_geo    : geostrophic V on V faces (NX,NY+1,NK)
%
% Notes
%   - Robust near coasts/steep topography: no cross-land/bottom differencing.
%   - Strict zeros below topography via 3-D face masks.
%--------------------------------------------------------------------------

[NX,NY,NK] = size(dz_T3d);

%--------------------------------------------------------------------------%
% 0) 3-D wet masks on T, and corresponding U/V face masks
%--------------------------------------------------------------------------%
mask3d_T = (dz_T3d > 0);                  % T-cell is wet at level k if dz>0

% U-face is wet if both adjacent T cells are wet at that level
umask3d = false(NX+1,NY,NK);
umask3d(2:NX,:,:)   = mask3d_T(1:NX-1,:,:) & mask3d_T(2:NX,:,:);
umask3d(1,:,:)      = mask3d_T(1,:,:);     % conservative at outer edge
umask3d(NX+1,:,:)   = mask3d_T(NX,:,:);

% V-face is wet if both adjacent T cells are wet at that level
vmask3d = false(NX,NY+1,NK);
vmask3d(:,2:NY,:)   = mask3d_T(:,1:NY-1,:) & mask3d_T(:,2:NY,:);
vmask3d(:,1,:)      = mask3d_T(:,1,:);
vmask3d(:,NY+1,:)   = mask3d_T(:,NY,:);

%--------------------------------------------------------------------------%
% 1) SSH -> surface geostrophic velocity on T cells
%     u_ssh_T = -(g/f) * d(ssh)/dy ,  v_ssh_T = (g/f) * d(ssh)/dx
%     Finite differences only across ocean neighbors on T.
%--------------------------------------------------------------------------%
dssh_dx = zeros(NX,NY);
dssh_dy = zeros(NX,NY);

% centered x-derivative (i=2..NX-1), only if both sides are ocean
for i = 2:NX-1
  for j = 1:NY
    if maskT(i-1,j) && maskT(i+1,j)
      dssh_dx(i,j) = (ssh_T(i+1,j) - ssh_T(i-1,j)) / (dx_T(i,j) + dx_T(i-1,j));
    end
  end
end
% centered y-derivative (j=2..NY-1), only if both sides are ocean
for i = 1:NX
  for j = 2:NY-1
    if maskT(i,j-1) && maskT(i,j+1)
      dssh_dy(i,j) = (ssh_T(i,j+1) - ssh_T(i,j-1)) / (dy_T(i,j) + dy_T(i,j-1));
    end
  end
end
% one-sided at outer boundaries (only if neighbor is ocean)
for j = 1:NY
  if maskT(1,j) && maskT(2,j)
    dssh_dx(1,j) = (ssh_T(2,j) - ssh_T(1,j)) / dx_T(1,j);
  end
  if maskT(NX-1,j) && maskT(NX,j)
    dssh_dx(NX,j) = (ssh_T(NX,j) - ssh_T(NX-1,j)) / dx_T(NX,j);
  end
end
for i = 1:NX
  if maskT(i,1) && maskT(i,2)
    dssh_dy(i,1) = (ssh_T(i,2) - ssh_T(i,1)) / dy_T(i,1);
  end
  if maskT(i,NY-1) && maskT(i,NY)
    dssh_dy(i,NY) = (ssh_T(i,NY) - ssh_T(i,NY-1)) / dy_T(i,NY);
  end
end

% surface geostrophic (on T)
u_ssh_T = -(g ./ f2d) .* dssh_dy .* maskT;   % eastward
v_ssh_T =  (g ./ f2d) .* dssh_dx .* maskT;   % northward

%--------------------------------------------------------------------------%
% 2) Thermal-wind shear on T cells: du/dz, dv/dz from horizontal rho grads
%     Only across wet neighbors at each level (respect topography).
%--------------------------------------------------------------------------%
du_dz_T = zeros(NX,NY,NK);
dv_dz_T = zeros(NX,NY,NK);

for k = 1:NK
  % x-direction density gradient (centered where both sides wet)
  for i = 2:NX-1
    for j = 1:NY
      if mask3d_T(i-1,j,k) && mask3d_T(i+1,j,k)
        drho_dx = (rho(i+1,j,k) - rho(i-1,j,k)) / (dx_T(i,j) + dx_T(i-1,j));
        dv_dz_T(i,j,k) =  (g / (rho0 * f2d(i,j))) * drho_dx;
      end
    end
  end
  % y-direction density gradient (centered where both sides wet)
  for i = 1:NX
    for j = 2:NY-1
      if mask3d_T(i,j-1,k) && mask3d_T(i,j+1,k)
        drho_dy = (rho(i,j+1,k) - rho(i,j-1,k)) / (dy_T(i,j) + dy_T(i,j-1));
        du_dz_T(i,j,k) = -(g / (rho0 * f2d(i,j))) * drho_dy;
      end
    end
  end
end

%--------------------------------------------------------------------------%
% 3) Vertical integration on T cells to get relative velocity (surface ref 0)
%--------------------------------------------------------------------------%
u_rel_T = zeros(NX,NY,NK);
v_rel_T = zeros(NX,NY,NK);

for k = 2:NK
  add_u = du_dz_T(:,:,k-1) .* dz_T3d(:,:,k-1);
  add_v = dv_dz_T(:,:,k-1) .* dz_T3d(:,:,k-1);
  % add only where layer k-1 exists
  add_u(~mask3d_T(:,:,k-1)) = 0;
  add_v(~mask3d_T(:,:,k-1)) = 0;

  u_rel_T(:,:,k) = u_rel_T(:,:,k-1) + add_u;
  v_rel_T(:,:,k) = v_rel_T(:,:,k-1) + add_v;
end

% safety: zero out dry T cells at each level
u_rel_T(~mask3d_T) = 0;
v_rel_T(~mask3d_T) = 0;

%--------------------------------------------------------------------------%
% 4) Enforce SSH reference on T, form total geostrophic on T
%     u_T(k) = u_rel_T(k) + [u_ssh_T - u_rel_T(surface)]
%--------------------------------------------------------------------------%
u_ref_T = u_ssh_T;    % surface reference on T
v_ref_T = v_ssh_T;

u_ref_T3d = repmat(u_ref_T, [1 1 NK]);
v_ref_T3d = repmat(v_ref_T, [1 1 NK]);

% zero reference where T cell is dry at that level
u_ref_T3d(~mask3d_T) = 0;
v_ref_T3d(~mask3d_T) = 0;

u_T = u_rel_T + u_ref_T3d;             % total geostrophic on T (per level)
v_T = v_rel_T + v_ref_T3d;

%--------------------------------------------------------------------------%
% 5) Single T -> U/V interpolation at the end (per level), with 3-D masks
%--------------------------------------------------------------------------%
u_geo = zeros(NX+1,NY,NK);
v_geo = zeros(NX,NY+1,NK);

for k = 1:NK
  % U faces: average adjacent T cells only if both are wet at level k
  for i = 2:NX
    for j = 1:NY
      if mask3d_T(i-1,j,k) && mask3d_T(i,j,k)
        u_geo(i,j,k) = 0.5*(u_T(i-1,j,k) + u_T(i,j,k));
      end
    end
  end
  % V faces: average adjacent T cells only if both are wet at level k
  for i = 1:NX
    for j = 2:NY
      if mask3d_T(i,j-1,k) && mask3d_T(i,j,k)
        v_geo(i,j,k) = 0.5*(v_T(i,j-1,k) + v_T(i,j,k));
      end
    end
  end
end

% enforce strict zeros below topography on faces
u_geo(~umask3d) = 0;
v_geo(~vmask3d) = 0;

% final NaN-guard
u_geo(~isfinite(u_geo)) = 0;
v_geo(~isfinite(v_geo)) = 0;
end
