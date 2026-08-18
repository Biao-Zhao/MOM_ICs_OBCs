function [A, lin_id, c0_diag] = build_poisson_T(H_T, mask_eff, area, dx_u, dy_v, dy_u, dx_v)

[NX,NY] = size(H_T);


He = zeros(NX+1,NY);
He(2:NX,:)   = 0.5*(H_T(1:NX-1,:)+H_T(2:NX,:));
He(1,:)      = H_T(1,:);
He(NX+1,:)   = H_T(NX,:);

Hn = zeros(NX,NY+1);
Hn(:,2:NY)   = 0.5*(H_T(:,1:NY-1)+H_T(:,2:NY));
Hn(:,1)      = H_T(:,1);
Hn(:,NY+1)   = H_T(:,NY);


wet    = find(mask_eff);
nw     = numel(wet);
lin_id = wet;

ij2k         = -ones(NX,NY);
ij2k(wet)    = 1:nw;


I = zeros(5*nw,1); J = I; V = I;
c0_diag = zeros(nw,1);
ptr = 0;

[is,js] = ind2sub([NX,NY], wet);

for kk = 1:nw
    i = is(kk); j = js(kk);
    aij = area(i,j);
    c0  = 0;

    % East face (i+1/2, j)
    if i < NX && mask_eff(i+1,j)
        ce = He(i+1,j) * dy_u(i+1,j) / ( dx_u(i+1,j) * aij );
        c0 = c0 + ce;
        kn = ij2k(i+1,j);
        ptr=ptr+1; I(ptr)=kk; J(ptr)=kn; V(ptr) = -ce;
    end

    % West face (i-1/2, j)
    if i > 1 && mask_eff(i-1,j)
        cw = He(i,  j) * dy_u(i,  j) / ( dx_u(i,  j) * aij );
        c0 = c0 + cw;
        kn = ij2k(i-1,j);
        ptr=ptr+1; I(ptr)=kk; J(ptr)=kn; V(ptr) = -cw;
    end

    % North face (i, j+1/2)
    if j < NY && mask_eff(i,j+1)
        cn = Hn(i,  j+1) * dx_v(i,  j+1) / ( dy_v(i,  j+1) * aij );
        c0 = c0 + cn;
        kn = ij2k(i,j+1);
        ptr=ptr+1; I(ptr)=kk; J(ptr)=kn; V(ptr) = -cn;
    end

    % South face (i, j-1/2)
    if j > 1 && mask_eff(i,j-1)
        cs = Hn(i,  j) * dx_v(i,  j) / ( dy_v(i,  j) * aij );
        c0 = c0 + cs;
        kn = ij2k(i,j-1);
        ptr=ptr+1; I(ptr)=kk; J(ptr)=kn; V(ptr) = -cs;
    end

    c0_diag(kk) = c0;
    ptr=ptr+1; I(ptr)=kk; J(ptr)=kk; V(ptr) = c0;
end

A = sparse(I(1:ptr), J(1:ptr), V(1:ptr), nw, nw);
end
