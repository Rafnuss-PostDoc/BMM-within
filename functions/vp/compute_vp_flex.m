function [vp,err] = compute_vp_flex(F, data, g, varargin)
%% Method using MATLAB lsqr function
% Can set various regularization, with a more flexible framework (building
% the actual A matrix to inverse.
% use compute_vp for a more typical setting
p = inputParser;
addParameter(p,'identity',0);
addParameter(p,'smoothness',1);
addParameter(p,'water',true);  %
addParameter(p,'tol', 1e-6); % same default as lsqr
addParameter(p,'maxit', 100); % same default as lsqr
addParameter(p,'x0', []); % same default as lsqr
addParameter(p,'max0', true); % same default as lsqr

parse(p,varargin{:});

[D, id_nnz] = format_data(data(:,1));

% Build the matrix linking a pixel in 3D to its vertical index
U = repelem(speye(g.sz(3)),prod(g.sz(1:2)),1);

% If water, do not include pixels over water
if p.Results.water
    load("data/onwater/Fonwater.mat");
    onwater = logical(Fonwater({g.f_lat(g.y),g.f_lon(g.x)}))'; clear Fonwater
    % figure; hold on; imagesc(g.f_lon(g.x),g.f_lat(g.y),onwater'); tmp=axis;borders('states','w'); axis equal tight; axis(tmp);
    onwater3D=repmat(onwater,1,1,g.sz(3));
    U(onwater3D(:)==1,:)=0;
end

% Build the inverse system
Ainv = F(id_nnz,:)*U;
b = D;

% Add identity regularization
if p.Results.identity > 0
    Lid = speye(g.sz(3));
    Ainv=[Ainv; p.Results.identity*Lid];
    b=[b; zeros(g.sz(3),1)];
    % b=[b; mean(b) .* ones(g.sz(3),1)];
end

% Add smoothness regularization
if p.Results.smoothness > 0
    Lh = LaplacianMatrix1D(g.sz(3));
    Lh = Lh - diag(sum(Lh,2)); % making sure the sum is zero. sould only apply boundary condition.

    Ainv = [Ainv; p.Results.smoothness.*Lh];
    b = [b; sparse(g.sz(3),1)];
end

% Inverse
vp = lsqr(Ainv, b, p.Results.tol, p.Results.maxit,[],[],p.Results.x0);

if p.Results.max0 && any(vp<0)
    warning("Vertical profile below 0. limit the vp to 0")
    vp = max(0,vp);
else
    warning("Vertical profile below 0.")
end

if nargout>1
    err = mean((F(id_nnz,:)*U*vp - D).^2);
    if p.Results.identity > 0
        err = [err; mean((Lid*vp - zeros(g.sz(3),1)).^2)];
    end
    if p.Results.smoothness > 0
        err = [err; mean((Lh*vp).^2)];
    end
end

% Vertical profile of difference
% e = ones(g.sz(3),1);
% T3diff = spdiags([-e, e], [0,1], g.sz(3)-1, g.sz(3));
% L3Dv = kron(T3diff,speye(g.sz(1)*g.sz(2)));
% d_diff = repelem(diff(Dgridzreg),g.sz(1)*g.sz(2),1);
% 
% alpha1=1;
% 
% Ainv=[F(id_nnz,:);...
%     alpha1*L3Dv];
% b=[D; d_diff];
% 
% Dinv = lsqr(Ainv,b,[],100,[],[],U*Dgridzreg);
% Dinv = reshape(Dinv,g.sz);
% 
% figure;
% h = slice(gLON,gLAT,gZ,Dinv,radar_lon,radar_lat,[100]);
% set(h,'edgecolor','none')
% clim([0 zmax]);
% hold on;


% HPR Horizontal profile of reflectivity
% % Create the W matrix W*d_grid = g_grid,z. With W=FU
% [Xsq,Ysq]=meshgrid(df.x,df.y);
% U2 = sparse(X(:) == Xsq(:)'&Y(:) == Ysq(:)');
% W2 = F(id_nnz,:)*U2;
% g.sz(3)2 = size(U2,2);
% % W = W./sum(W,2,'omitnan');
% 
% % no regularization
% DgridxyI = reshape(W2\D,size(Xsq));
% 
% figure; 
% nexttile; hold on; box on; axis equal
% imagesc(DgridxyI)
