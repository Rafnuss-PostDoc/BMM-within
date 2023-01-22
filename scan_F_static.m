addpath('/Users/raphael/Documents/GitHub/wsrlib');
addpath(genpath('functions/'))
wsrlib_setup();

%% Choise of radar
if false
    load('../BMM-US/data/density/inference-trans.mat')
    vidT = vidTS;
    vidT(~vidTS_id)=nan;

    i_r=radar.name=="KBYX";
    figure; plot(g.time, vidT(:,i_r))
end

%% Selection of exemple
scans = table( ...
    ["KOKX","KHGX","KLOT","KNQA","KDAX","KPAH"]',...
    datetime(["24-Oct-2020 3:47", "24-Apr-2020 2:44","16-May-2020 1:55","29-Sep-2020 4:06","25-Apr-2020 4:19","29-Sep-2020 0:17"])',...
    'VariableNames',["station" "time"]);


%% Download files

%% Process 

i_s=2;

% Select fild
scans_filepath = load_rd(scans.station(i_s), scans.time(i_s));

% Read and process file
rd = read_rd(scans_filepath, ...
    'removeRain', 'rho',... % "mistnet" slower
    'removeOuliars', 3,...
    'shape',true,...
    'smoothing', true,...
    'limitElev', 5,...
    'minRange', 5000,...
    'maxAlt', 3000,...
    'data_res', 1, ...
    'par',false);

% Load df
r_lim = 50000;
gridsize = num2str(r_lim/1000)+"km";
[df, rd] = loadDf(gridsize, rd);

%% Prepare the grid

% Grid
radar_lon = cellfun(@(x) x.lon, rd(:,1));
radar_lat = cellfun(@(x) x.lat, rd(:,1));

g.dlat = deg2km(distance(mean(radar_lat),mean(radar_lon),mean(radar_lat)+1,mean(radar_lon)))*1000;
g.dlon = deg2km(distance(mean(radar_lat),mean(radar_lon),mean(radar_lat),mean(radar_lon)+1))*1000;

xmax = range(radar_lon)*g.dlon+2*max(df{1}.x);
xmax = ceil(xmax/10000)*10000;
ymax = range(radar_lat)*g.dlat+2*max(df{1}.y);
ymax = ceil(ymax/10000)*10000;
dx = unique(diff(df{1}.x));
g.x=(0:ceil(xmax/dx)+3)*dx;
dy = unique(diff(df{1}.y));
g.y=(0:ceil(ymax/dy)+3)'*dy;
g.z=df{1}.z;
g.sz = [numel(g.x) numel(g.y) numel(g.z)];

g.lon0 = floor((min(radar_lon)-max(df{1}.x)/g.dlon)*10)/10; 
g.lat0 = floor((min(radar_lat)-max(df{1}.y)/g.dlat)*10)/10;

% get lat lon coordinates
g.f_lon = @(X) X./g.dlon+g.lon0;
g.f_lat = @(Y) Y./g.dlat+g.lat0;
f_X = @(lon) (lon-g.lon0)*g.dlon;
f_Y = @(lat) (lat-g.lat0)*g.dlat;

g.lat=g.f_lat(g.y);
g.lon=g.f_lon(g.x);

[~,radar_lon_id]=min(abs(g.lon-radar_lon));
[~,radar_lat_id]=min(abs(g.lat-radar_lat));

[gX,gY,~]=meshgrid(g.x,g.y,g.z);
[gLON,gLAT,gZ]=meshgrid(g.lon,g.lat,g.z);


%% 
[F, data, info] = rd2F_mosaic(rd, df, g);
szf=size(F);

%%

i_t=1;
% retrieve id of non-zeros
id_nnz = find(data(:,i_t));

% keep only non-zeros
D = full(data(id_nnz,i_t));

% find the optimal max density for plot
zmax = quantile(D,.999);

% find the index of the data in the dimension of radar, df, range and
% azimuth for the non-zeron data
[id_range, id_az, id_r, id_df]=ind2sub(info.sz_rd,id_nnz);

% 
df_range = cellfun(@(x) x.range,df,'UniformOutput',false);
tmp = df_range(id_df);
R = cellfun(@(x,y) x(y),tmp, num2cell(id_range));

df_az = cellfun(@(x) x.az,df,'UniformOutput',false);
tmp = df_az(id_df);
A = cellfun(@(x,y) x(y),tmp, num2cell(id_az));

df_e = cellfun(@(x) x.elev,df);
E = df_e(id_df);
clear tmp


% figure; 
tmp = reshape(full(sum(F,1)), g.sz);

figure;
h = slice(gLON,gLAT,gZ,tmp,median(g.lon),median(g.lat),[100]);
set(h,'edgecolor','none')

%% See data

figure; tiledlayout('flow','TileSpacing','tight','Padding','tight');
for i_elev=1:min(10,numel(rd{1}.dens.sweeps))
    nexttile;
    sw = rd{1}.dens.sweeps(i_elev);
    imagesc(sw.az, sw.range, sw.data, 'AlphaData',~isnan(sw.data))
    ylim([0 r_lim]); colorbar; 
    caxis([0 zmax]);
end


dr=4;
daz=2;
figure; tiledlayout('flow','TileSpacing','none','Padding','none');
for i_elev=1:min(3,numel(rd{1}.dens.sweeps))
    
    sw = rd{1}.dens.sweeps(i_elev);
    id_max =find(sw.range<r_lim,1,'last');
    [az,r]=meshgrid(sw.az(1:daz:end),sw.range(1:dr:id_max));
    dens = sw.data(1:dr:id_max, 1:daz:end );

    nexttile;
    polarscatter(deg2rad(az(:)),r(:),r(:)/8000,dens(:),'filled')
    set(gca,'ThetaZeroLocation',"top",'ThetaDir','clockwise'); 
    set(gca,"ThetaTickLabel","","rticklabels","")
    set(gca,"RLim",[0 80000])
    caxis([0 zmax]);
    title(round(sw.elev,1))
end


dr=10;
daz=5;
figure; tiledlayout('flow','TileSpacing','tight','Padding','tight'); hold on;
for i_elev=1:min(6,numel(df_e))
    sw = rd{1}.dens.sweeps(i_elev);
    id_max =find(sw.range<r_lim,1,'last');
    [az,r]=meshgrid(sw.az(1:daz:end),sw.range(1:dr:id_max));
    dens = sw.data(1:dr:id_max, 1:daz:end );

    [ s, h ] = slant2ground( r, sw.elev );
    [x1,x2] = pol2cart(cmp2pol(az),r);

    scatter3(x1(:),x2(:),h(:),10,dens(:),'filled')% r(id)/5000
    
end
caxis([0 zmax]); box on; grid on; zlim([0 3000])
xlim([-150 150]*1000);ylim([-150 150]*1000)
xticks([-150 0 150]*1000); yticks([-150 0 150]*1000)
zticks(0:1000:3000)
view(3)
% Ax=gca; Ax.XColor = 'w'; Ax.YColor = 'W'; Ax.ZColor = 'W'; Ax.Color = [50 50 50]/255;
% exportgraphics(gcf,'test-3.png','BackgroundColor',[50 50 50]/255)
% 
% figure; tiledlayout('flow','TileSpacing','tight','Padding','tight'); hold on;
% for i_elev=1:min(20,numel(rd.dz.sweeps))
%     dr=1;
%     daz=1;
%     data = squeeze(D(i_elev,1:dr:end,1:daz:end));
%     r = squeeze(R(i_elev,1:dr:end,1:daz:end));
%     az = squeeze(A(i_elev,1:dr:end,1:daz:end));
% 
%     [ s, h ] = slant2ground( r, rd.dz.sweeps(i_elev).elev );
%     [x1,x2] = pol2cart(cmp2pol(az),r);
% 
%     id = r<r_lim & x2<1000&x2>-1000;
%     scatter(x1(id),h(id),10*sqrt(r(id)/800),data(id),'filled')
% end
% axis tight; ylim([0 3000]); caxis([0 zmax]); grid on; colorbar;


%% Forward

dfz = squeeze(g.z);
vp = exp(-1/200*dfz);
% vp = 100*(dfz<1000);
figure; plot(vp,dfz,'-or'); xlabel('Density'); ylabel('Altitude'); title('vertical profile')

% bd=100*ones(g.sz); % bird/km^3
bd = repmat(permute(vp,[2 3 1]),numel(g.x),numel(g.y),1);
Dobs=reshape(F*bd(:),info.sz_rd);

figure; tiledlayout('flow','TileSpacing','tight','Padding','tight');
for i_elev=1:min(6,numel(rd.dz.sweeps))
    nexttile;
    imagesc(df{i_elev}.az, df{i_elev}.range, Dobs(:,:,:,i_elev), 'AlphaData',~isnan(Dobs(:,:,:,i_elev)))
    ylim([0 50000]); colorbar; 
    %caxis([0 zmax]);
end

dr=5;
daz=2;
figure; tiledlayout('flow','TileSpacing','tight','Padding','tight'); hold on;
for i_elev=1:min(6,numel(df_e))
    dens = Dobs(:,:,:,i_elev);
    id_max =find(df{i_elev}.range<r_lim,1,'last');
    [az,r]=meshgrid(df{i_elev}.az(1:daz:end),df{i_elev}.range(1:dr:id_max));
    dens = dens(1:dr:id_max, 1:daz:end );

    [ s, h ] = slant2ground( r, df{i_elev}.elev );
    [x1,x2] = pol2cart(cmp2pol(az),r);

    scatter3(x1(:),x2(:),h(:),10,dens(:),'filled')% r(id)/5000
    
end
caxis([0 zmax]); box on; grid on; zlim([0 3000])

%% interpolation

% 2.5D interpolation

[dens, x1, x2, x3] = radar2mat2(rd{1}, 'fields', {'dens'},'coords','cartesian','dim',numel(g.x),'r_max',g.x(end),'r_min',5000);

figure('position',[0 0 1200 750]); tiledlayout('flow','TileSpacing','tight','Padding','tight');
for i_elev=1:min(6,numel(rd{1}.dens.sweeps))
    nexttile;
    imagesc(x1, x2, dens.dens(:,:,i_elev),'alphadata',~isnan(dens.dens(:,:,i_elev))); set(gca,'ydir','normal')
    caxis([0 zmax]); 
    axis equal tight; title(round(rd{1}.dens.sweeps(i_elev).elev,2)+" °")
    xlim([-80 80]*1000); ylim([-80 80]*1000); xticklabels(''); yticklabels('')
end


% % 3d interpolation
% F=griddedInterpolant({x1 x2 x3},dens.dens,'linear','none');
% Dinterp=F(gX,gY,T); Dinterp=permute(Dinterp,[2 1 3]);
% figure('position',[0 0 1200 750]);  tiledlayout('flow','TileSpacing','tight','Padding','tight');
% for i_elev=1:2:12
%     nexttile;
%     imagesc(df.x, df.y, Dinterp(:,:,i_elev),'alphadata',~isnan(Dinterp(:,:,i_elev))); set(gca,'ydir','normal')
%     caxis([0 z_to_refl(idb(30))/11]); axis equal tight; title(g.z(i_elev)+" m")
%     xlim([-80 80]*1000); ylim([-80 80]*1000); xticklabels(''); yticklabels('')
% end
% 
% figure;
% h = slice(X,Y,Z,Dinterp,[0],[0],[100]);
% set(h,'edgecolor','none')
% caxis([0 30]);

















%% Identity

% initial condition
% x0=Dinterp(:);

% regulatization strength
alpha = 0.1;

% nothing
% Ainv=[F(id_nnz,:)];
% b=[D];

% Identity
Lid = speye(prod(g.sz));
Ainv=[F(id_nnz,:);alpha*Lid];
b=[D;zeros(prod(g.sz),1)];

Dinvtmp = lsqr(Ainv,b,[],1000);
Dinv=reshape(Dinvtmp,g.sz);


figure;
h = slice(gLON,gLAT,gZ,Dinv,radar_lon,radar_lat,[100]);
set(h,'edgecolor','none')
caxis([0 zmax]);

hold on;
[ s, h ] = slant2ground( R(id), E(id) );
[x1,x2] = pol2cart(cmp2pol(A(id)),R(id));
id2 = R(id)<r_lim & x2<1000&x2>-1000 & h<2000;
scatter3(g.f_lon(g.x(radar_lon_id)+x1(id2)),g.f_lat(g.y(radar_lat_id)+x2(id2)),h(id2),10,D(id(id2)),'filled')
id2 = R(id)<r_lim & x1<1000&x1>-1000 & h<2000;
scatter3(g.f_lon(g.x(radar_lon_id)+x1(id2)),g.f_lat(g.y(radar_lat_id)+x2(id2)),h(id2),10,D(id(id2)),'filled')



%% Smoothing Laplacian

% LaplacianMatrix1D => e = ones(n,1); L1 = spdiags([-e, 2*e, -e], [-1,0,1], n, n);

T1 = LaplacianMatrix1D(g.sz(1));
% T1 = T1 + .7*LaplacianMatrix1D(g.sz(1),3);
% T1 = T1 + .3*LaplacianMatrix1D(g.sz(1),20);
T1 = T1 - diag(sum(T1,2));
T2 = LaplacianMatrix1D(g.sz(2));
% T2 = T2 + .7*LaplacianMatrix1D(g.sz(2),3);
% T2 = T2 + .3*LaplacianMatrix1D(g.sz(2),20);
T2 = T2 - diag(sum(T2,2));
T3 = 2*LaplacianMatrix1D(g.sz(3));
% T3 = T3 + 1*LaplacianMatrix1D(g.sz(3),3);
% T3 = T3 + .5*LaplacianMatrix1D(g.sz(3),10);
T3 = T3 - diag(sum(T3,2));
coeff_horver = .4; % (+) more horizontal smoother, (-) more vertical smoothing, 0 equal smoothing. between 0 and .5
L2D = kron(T2, speye(g.sz(1))) + kron(speye(g.sz(2)), T1);
L3D = (kron(T3, speye(g.sz(1)*g.sz(2)))*(1-coeff_horver*2) + kron(speye(g.sz(3)), L2D)*(1+coeff_horver))/6;

% cut at coast
load("data/onwater/Fonwater.mat");
onwater = logical(Fonwater({g.f_lat(g.y),g.f_lon(g.x)}));
%figure; hold on; imagesc(g.f_lon(g.x),g.f_lat(g.y),onwater); tmp=axis;borders('states','w'); axis equal tight; axis(tmp);

L2Dwater=L2D;
[i,j] = find(L2Dwater);
id = onwater(i) ~= onwater(j);
L2Dwater(i(id),j(id)) = 0;
L2Dwater = L2Dwater - diag(sum(L2Dwater,2));
% figure; hold on; imagesc(f_lon(X(1,:,1))',f_lat(Y(:,1,1)),reshape(diag(L2Dsame),g.sz(1),g.sz(2))); tmp=axis;borders('states','w'); axis equal tight; axis(tmp);
L3Dwater = (kron(T3, speye(g.sz(1)*g.sz(2)))*(1-coeff_horver*2) + kron(speye(g.sz(3)), L2Dwater)*(1+coeff_horver))/6;


alpha=1;
Ainv=[F(id_nnz,:); alpha*L3Dwater];
b=[D; sparse(szf(2),1)];

Dinv = lsqr(Ainv,b,[],500,[],[],Dinv(:));
Dinv = reshape(Dinv,g.sz);

figure;
h = slice(gLON,gLAT,gZ,Dinv,radar_lon,radar_lat,[100]);
set(h,'edgecolor','none')
caxis([0 zmax]);

hold on;
[ s, h ] = slant2ground( R, E );
[x1,x2] = pol2cart(cmp2pol(A),R);
id = R<r_lim & x2<1000&x2>-1000 & h<2000;
scatter3(g.f_lon(g.x(radar_lon_id)+x1(id)),g.f_lat(g.y(radar_lat_id)+x2(id)),h(id),10,D(id),'filled')
id = R<r_lim & x1<1000&x1>-1000 & h<2000;
scatter3(g.f_lon(g.x(radar_lon_id)+x1(id)),g.f_lat(g.y(radar_lat_id)+x2(id)),h(id),10,D(id),'filled')

    
figure; imagesc(g.lon,g.lat,sum(Dinv,3)); set(gca,'ydir', 'normal')





%% Inverting the vertical profile

[ ~, H ] = slant2ground( R(:), E(:) );
[~,~,gr] = histcounts(H,[g.z(1)-diff(g.z(1:2)) ;  squeeze(g.z)]+diff(g.z(1:2))/2);
Dh=zeros(numel(g.z), 1);
[G,ID] = findgroups(gr);
Dh(ID)=splitapply(@nanmean, D, G);
% figure; plot(Dh,squeeze(g.z),'-o')

% Create the W matrix W*d_grid = g_grid,z. With W=FU
% U = sparse(gZ(:) == squeeze(g.z)');
U = repelem(speye(g.sz(3)),prod(g.sz(1:2)),1);
onwater3D=repmat(onwater,1,1,g.sz(3));
% Iland = speye(szf(2)); Iland(onwater3D(:),:)=0;
Uland = U; Uland(onwater3D(:)==1,:)=0;
W = F*Uland;

Lh = LaplacianMatrix1D(g.sz(3));
Lh = Lh - diag(sum(Lh,2));

alpha=5;
Ainv=[W(id_nnz,:);alpha^2*Lh];
b=[D;sparse(g.sz(3),1)];

Dgridzreg = lsqr(Ainv,b,[],100,[],[],Dh);
DgridzI = lsqr(W(id_nnz,:),D,[],100,[],[],Dh);

figure; hold on; box on; grid on;
plot(Dh,squeeze(g.z),'-o')
plot(DgridzI,squeeze(g.z),'-o')
plot(Dgridzreg,squeeze(g.z),'-o')
xlabel('density'); ylabel('Elevation');
legend('Mean of d_{ob}','Inversion','Inversion with Laplacian regularization')


%% HPR Horizontal profile of reflectivity

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

%% Vertical profile of difference
e = ones(g.sz(3),1);
T3diff = spdiags([-e, e], [0,1], g.sz(3)-1, g.sz(3));
L3Dv = kron(T3diff,speye(g.sz(1)*g.sz(2)));
d_diff = repelem(diff(Dgridzreg),g.sz(1)*g.sz(2),1);

alpha1=1;
alpha2=.1;
alpha3=0;

Ainv=[F(id_nnz,:);...
    alpha1*L3Dv];
b=[D; d_diff];

Dinv = lsqr(Ainv,b,[],100,[],[],U*Dgridzreg);
Dinv = reshape(Dinv,g.sz);

figure;
h = slice(gLON,gLAT,gZ,Dinv,radar_lon,radar_lat,[100]);
set(h,'edgecolor','none')
caxis([0 zmax]);
hold on;









%% VPR of ratio with 2-step
tmp = Dgridzreg;
tmp(tmp<0)=0;
Dgridratio=tmp./sum(tmp);

% V = sparse(X(:) == reshape(X(:,:,1),1,[]) & Y(:) == reshape(Y(:,:,1),1,[]));
% Q = repmat(V,1,numel(g.z)).*repelem(Dgridratio',1,prod(g.sz(1:2)))-speye(prod(g.sz));

V = repmat(speye(prod(g.sz(1:2))),g.sz(3),g.sz(3));
Q = V.*repelem(Dgridratio',1,prod(g.sz(1:2)))-speye(prod(g.sz));

alpha=1;
Ainv=[F(id_nnz,:);...
    alpha*Q];
b=[D; sparse(szf(2),1)];

Dinv = lsqr(Ainv,b,[],100,[],[],U*Dgridzreg);
Dinv = reshape(Dinv,g.sz);

figure;
h = slice(gLON,gLAT,gZ,Dinv,radar_lon,radar_lat,[100]);
set(h,'edgecolor','none')
caxis([0 zmax]);

hold on;
[ s, h ] = slant2ground( R, E );
[x1,x2] = pol2cart(cmp2pol(A),R);
id = R<r_lim & x2<1000&x2>-1000 & h<2000;
scatter3(g.f_lon(g.x(radar_lon_id)+x1(id)),g.f_lat(g.y(radar_lat_id)+x2(id)),h(id),10*sqrt(R(id)/800),D(id),'filled')
id = R<r_lim & x1<1000&x1>-1000 & h<2000;
scatter3(g.f_lon(g.x(radar_lon_id)+x1(id)),g.f_lat(g.y(radar_lat_id)+x2(id)),h(id),10*sqrt(R(id)/800),D(id),'filled')















%% All combined
alpha0=1;
alpha1=.1;%.25;
alpha2=1;%.5;

% 1-step, separate-abs, smoothness
% Ainv=[F(id_nnz,:) sparse(szf(1),g.sz(3)); ...
%     sparse(szf(1),szf(2),1) W; ...
%     alpha1.*speye(szf(2)) -alpha1'.*U;...
%     alpha2*L3D sparse(szf(2), g.sz(3))];
% b=[D(:); D(:); sparse(2*szf(2),1)];


% 2-step, seprate-abs, smoothness 
 Ainv=[F(id_nnz,:); ...
     alpha1.*speye(szf(2));...
     alpha2*L3Dwater];
 b=[D(:); alpha1'.*U*Dgridzreg ; sparse(szf(2),1)];


% 1-step, combined-abs, smoothness, 
% Ainv=[F(id_nnz,:) sparse(szf(1),g.sz(3)); ...
%     alpha1*speye(szf(2)) -alpha1*U; ...
%     alpha2*L3D sparse(szf(2), g.sz(3))];
% b=[D(:); sparse(2*szf(2),1)];


% 1-step, combined-abs, smoothness, water
% Ainv=[F(id_nnz,:) sparse(szf(1),g.sz(3)); ...
%     alpha1*Iland -alpha1*Uland; ...
%     alpha2*L3Dwater sparse(szf(2), g.sz(3))];
% b=[D(:); sparse(2*szf(2),1)];


% ratio
% Ainv=[F(id_nnz,:);...
%     alpha1*Q;...
%     alpha2*L3Dwater];
% b=[D; sparse(2*szf(2),1)];


% Solving
Dinvtmp = lsqr(Ainv,b,[],500,[],[],U*Dgridzreg);
Dinv=reshape(Dinvtmp(1:szf(2)),g.sz);
% Dinv = permute(Dinv,[2 1 3]);
% either... 
Dgridz = Dgridzreg;
%Dgridz = Dinvtmp(szf(2)+(1:g.sz(3)));

%%
%exportgraphics(gcf,filename+"_a1"+alpha1+alpha2+".png") 
    

figure;
h = slice(gX,gY,gZ,Dinv,g.x(radar_lon_id),g.y(radar_lat_id),[100]);
set(h,'edgecolor','none')
caxis([0 zmax]);

hold on;
[ s, h ] = slant2ground( R, E );
[x1,x2] = pol2cart(cmp2pol(A),R);
id = R<r_lim & x2<1000&x2>-1000 & h<2000;
scatter3(g.f_lon(g.x(radar_lon_id)+x1(id)),g.f_lat(g.y(radar_lat_id)+x2(id)),h(id),10*sqrt(R(id)/800),D(id),'filled')
id = R<r_lim & x1<1000&x1>-1000 & h<2000;
scatter3(g.f_lon(g.x(radar_lon_id)+x1(id)),g.f_lat(g.y(radar_lat_id)+x2(id)),h(id),10*sqrt(R(id)/800),D(id),'filled')


figure; imagesc(sum(Dinv,3))

%set(h,'edgecolor','w')
caxis([0 zmax]); view(3); xlim([-80 80]*1000); ylim([-80 80]*1000)
box on; grid on; zlim([0 3000])
xticks([-80 0 80]*1000); yticks([-80 0 80]*1000)
zticks(0:1000:3000)
view(3)
%Ax=gca; Ax.XColor = 'w'; Ax.YColor = 'W'; Ax.ZColor = 'W'; Ax.Color = [50 50 50]/255;
% exportgraphics(gcf,'test-4.png','BackgroundColor',[50 50 50]/255)



figure; hold on; 
imagesc(f_lon(X(1,25:135,1)),f_lat(Y(25:135,1,1)),sum(Dinv(25:135,25:135,:),3)'); axis tight equal; a_xis=axis;
borders('states','w'); axis(a_xis); %caxis([0 7*zmax]);
Ax=gca; Ax.XColor = 'w'; Ax.YColor = 'W'; Ax.ZColor = 'W'; Ax.Color = [50 50 50]/255;
% exportgraphics(gcf,'test-8.png','BackgroundColor',[50 50 50]/255)



%%

figure('position',[0 0 1400 800]); tiledlayout(3,5,'TileSpacing','tight','Padding','tight')
z=squeeze(g.z)/1000;
X2=X/1000; Y2=Y/1000; Z2=Z/1000;
[ xrd, yrd, zrd ] = radar2xyz( R2, A2, E2 );

nexttile; hold on; box on; grid on; axis equal
err1 = F(id_nnz,:) * Dinvtmp(1:szf(2)) - D(:);%sqrt(mean((F(id_nnz,:) * Dinv(:) - D(:)).^2)); 
plot(D(:), F(id_nnz,:) * Dinvtmp(1:szf(2)), '.k'); plot([0 max(D(:))],[0 max(D(:))],'r','linewidth',2)
% xline(mean(D(:)),'w'); yline(mean(F(id_nnz,:) * Dinv(:)),'w')
% xline(quantile(D(:),[.5 .9 .99]),'w'); yline(quantile(F(id_nnz,:) * Dinv(:),[.5 .9 .99]),'w')
% histogram(err1); set(gca,'yscale','log')
title("||Fd_{true}-d_{obs}|| = "+round(sqrt(sum(err1.^2))/1000) +"K"); 
xlabel('d_{obs}'); ylabel('Fd_{true}');

nexttile; hold on; box on; grid on; colorbar;
title("||Fd_{true}-d_{obs}||" +" avg over az."); 
eu = unique(E2);
for i_elev=1:numel(eu)
    id = find(E2==eu(i_elev)&R2<80000);
    R2b = round(R2(id),-3);
    [G,ID]=findgroups(R2b);
    [ s, h ] = slant2ground( ID, eu(i_elev) );
    tmp = splitapply(@mean,err1(id),G);
    scatter(s,h,[],tmp,'filled')
end
axis tight; ylim([0 3000]); xlim([0 80000]); caxis([-100 100]);
xlabel('R'); ylabel('Z'); xticklabels(''); yticklabels('')

nexttile; 
id = find(E2==min(E2)&R2<80000);
R2b = round(R2(id),-3);
A2b = round(A2(id)/2)*2;
[G,IDR,IDA]=findgroups(R2b,A2b);
tmp = splitapply(@mean,err1(id),G);
polarscatter(deg2rad(IDA),IDR,IDR/8000,tmp,'filled')
set(gca,'ThetaZeroLocation','top','Thetadir','clockwise','ThetaTick',[],'RTick',[])
caxis([-100 100]); title("||Fd_{true}-d_{obs}||" + " lowest elevation"); 


nexttile; hold on; box on; grid on; box on; 
err2 = Dinvtmp(1:szf(2)) - U*Dgridz(:); %sqrt(mean((Dinv(:) - U*Dgridz(:)).^2));
%plot(U*Dgridz(:), Dinvtmp(1:szf(2)), '.k'); plot([0 max(D(:))],[0 max(D(:))],'color',[.7 .7 .7])
% histogram(err2); set(gca,'yscale','log')
tmp = mean(reshape(err2,size(X)),3);
tmp(onwater)=nan;
imagesc(tmp,'alphadata',~isnan(tmp))
% boxplot(reshape(Dinv,[],size(Dinv,3)),'PlotStyle','compact','Orientation','horizontal')%,'Whisker',10)
% plot(Dgridz,1:31,'-or','linewidth',2);
title("||d_{true}-Ud_{VPR}|| = "+round(sqrt(sum(err2.^2))/1000) +"K"+ " | \alpha_1="+full(mean(alpha1))); 
% xlabel('Density'); ylabel('Z'); xlim([0 zmax]);
xlabel('X'); ylabel('Y'); caxis([-100 100]); yticklabels(''); xticklabels(''); axis equal tight;

nexttile; hold on; box on; grid on;
err3 = L3D*Dinvtmp(1:szf(2));% sqrt(mean((L3D*Dinv(:)).^2));
% plot(U*Dgridz(:), Dinv(:), '.k'); plot([0 max(D(:))],[0 max(D(:))],'color',[.7 .7 .7])
histogram(err3); set(gca,'yscale','log'); xline(0,'r','linewidth',2)
title("||L_{3D}d_{true}|| = "+round(sqrt(sum(err3.^2))/1000) +"K"+ " | \alpha_2="+alpha2); 
xlabel('Density'); ylabel('Histogram'); yticklabels('')

nexttile([2 1]); hold on; box on; grid on; title('VPR')
tmp = reshape(Dinv,[],size(Dinv,3));
id = randsample(size(tmp,1),500);      
p1=plot(tmp(id,:),z,'-o','color',[.7 .7 .7]);
p5=plot(mean(tmp),z,'-s','color',[.3 .3 .3],'linewidth',2);
p2=plot(Dh,(50:100:2950)/1000,'-o','linewidth',2);
p3=plot(DgridzI,z,'-o','linewidth',2);
p4=plot(Dgridz,z,'--k','linewidth',2);
legend([p1(1) p5 p2 p3 p4],'Radom VPR', 'Avg VPR of inv. data','Avg of raw data','Inverted VPR only','Inverted VPR combined');
xlabel('Bird density'); ylabel('Z'); xlim([0 zmax]); yticklabels('')

nexttile([1 2]); hold on; box on; grid on; title('Y=0')
imagesc(squeeze(X2(1,:,1))',squeeze(Z2(1,1,:))',squeeze(Dinv(:,round(end/2),:))')
id = R2<r_lim & zrd<3000 & abs(yrd)<1000;
scatter(xrd(id)/1000, zrd(id)/1000,20,D(id),'filled')
caxis([0 zmax]); xlim([-80 80]); ylim([0 3]); %xticklabels(''); yticklabels('')
ylabel('Z'); xlabel('X');xticklabels(''); yticklabels('')

nexttile([1 2]); hold on; box on; grid on; title('X=0')
imagesc(squeeze(Y2(:,1,1))',squeeze(Z2(1,1,:))',squeeze(Dinv(round(end/2),:,:))')
id = R2<r_lim & zrd<3000 & abs(xrd)<1000;
scatter(yrd(id)/1000, zrd(id)/1000,20,D(id),'filled')
caxis([0 zmax]); xlim([-80 80]); ylim([0 3]); %xticklabels(''); yticklabels('')
ylabel('Z'); xlabel('Y');xticklabels(''); yticklabels('')

nexttile; hold on; box on; grid on; title('Z=200m')
imagesc(squeeze(X2(1,:,1))',squeeze(Y2(:,1,1))',squeeze(Dinv(:,:,2))'); 
caxis([0 zmax]); xlim([-80 80]); ylim([-80 80]);xticklabels(''); yticklabels('')
ylabel('Y'); xlabel('X'); axis square

nexttile; hold on; box on; grid on; title('Z=400m')
imagesc(squeeze(X2(1,:,1))',squeeze(Y2(:,1,1))',squeeze(Dinv(:,:,4))')
caxis([0 zmax]); xlim([-80 80]); ylim([-80 80]); xticklabels(''); yticklabels('')
ylabel('Y'); xlabel('X'); axis square

nexttile; hold on; box on; grid on; title('Z=800m')
imagesc(squeeze(X2(1,:,1))',squeeze(Y2(:,1,1))',squeeze(Dinv(:,:,8))')
caxis([0 zmax]); xlim([-80 80]); ylim([-80 80]); xticklabels(''); yticklabels('')
ylabel('Y'); xlabel('X'); axis square


nexttile; hold on; box on; grid on; title('Z=1000m')
imagesc(squeeze(X2(1,:,1))',squeeze(Y2(:,1,1))',squeeze(Dinv(:,:,10))')
caxis([0 zmax]); xlim([-80 80]); ylim([-80 80]); xticklabels(''); yticklabels('')
ylabel('Y'); xlabel('X'); axis square



%%

figure; hold on
h = slice(f_lon(X),f_lat(Y),Z,Dinv,f_lon([0]),f_lat([0]),[100]);
set(h,'edgecolor','none')

dr=3;
daz=1;
for i_elev=1:min(6,numel(rd.dz.sweeps))

    data = squeeze(D(i_elev,1:dr:end,1:daz:end));
    r = squeeze(R(i_elev,1:dr:end,1:daz:end));
    az = squeeze(A(i_elev,1:dr:end,1:daz:end));

    [s, h] = slant2ground( r, rd.dz.sweeps(i_elev).elev );
    [x1,x2] = pol2cart(cmp2pol(az), r);

    id = r<r_lim & x2<1000&x2>-1000;
    scatter3(f_lon(x1(id)),f_lat(x2(id)),h(id),10*sqrt(r(id)/800),data(id),'filled')
    
end
axis tight; zlim([0 3000]); caxis([0 zmax]); view(3)



