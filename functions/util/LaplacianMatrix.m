function L3D = LaplacianMatrix(g, varargin)

p = inputParser;
addParameter(p,'water',true); 
addParameter(p,'coeff_horver',0); % (+) more horizontal smoother, (-) more vertical smoothing, 0 equal smoothing. between 0 and .5
parse(p,varargin{:});


T1 = LaplacianMatrix1D(g.sz(1));
% T1 = T1 + .7*LaplacianMatrix1D(g.sz(1),3);
% T1 = T1 + .3*LaplacianMatrix1D(g.sz(1),20);
T1 = T1 - diag(sum(T1,2));

T2 = LaplacianMatrix1D(g.sz(2));
% T2 = T2 + .7*LaplacianMatrix1D(g.sz(2),3);
% T2 = T2 + .3*LaplacianMatrix1D(g.sz(2),20);
T2 = T2 - diag(sum(T2,2));

T3 = LaplacianMatrix1D(g.sz(3));
% T3 = T3 + 1*LaplacianMatrix1D(g.sz(3),3);
% T3 = T3 + .5*LaplacianMatrix1D(g.sz(3),20);
T3 = T3 - diag(sum(T3,2));

L2D = kron(T2, speye(g.sz(1))) + kron(speye(g.sz(2)), T1);

if p.Results.water
    % Load water mask
    load("data/onwater/Fonwater.mat");
    onwater = logical(Fonwater({g.f_lat(g.y),g.f_lon(g.x)}))'; clear Fonwater
    % figure; hold on; imagesc(g.f_lon(g.x),g.f_lat(g.y),onwater'); tmp=axis;borders('states','w'); axis equal tight; axis(tmp);
    
    L2Dwater=L2D;
    [i,j] = find(L2Dwater);
    id = onwater(i) ~= onwater(j);
    L2Dwater(i(id),j(id)) = 0;
    L2Dwater = L2Dwater - diag(sum(L2Dwater,2));
    % figure; hold on; imagesc(g.lon,g.lat,reshape(diag(L2Dwater),g.sz(1),g.sz(2))); tmp=axis;borders('states','w'); axis equal tight; axis(tmp);
    
    L2D = L2Dwater;
end
L3D = (kron(T3, speye(g.sz(1)*g.sz(2)))*(1-p.Results.coeff_horver*2) + kron(speye(g.sz(3)), L2D)*(1+p.Results.coeff_horver))/6;
end