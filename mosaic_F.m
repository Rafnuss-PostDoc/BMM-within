addpath('/Users/raphael/Documents/GitHub/wsrlib');
wsrlib_setup();
addpath(genpath('functions/'))
load('../BMM-US/data/density/inference-trans.mat','radar')

%% Selection of case

if false
    i_r = ismember(radar.name,scan_radar);
    figure; tiledlayout('flow','TileSpacing','none','Padding','none');
    plot(radar.lon(i_r), radar.lat(i_r),'or'); hold on;
    rectangle('Position',[grid_extend(1) grid_extend(3) grid_extend([2,4])-grid_extend([1,3])])
    axis equal tight;  a_xis=axis;
    plot(radar.lon(~i_r), radar.lat(~i_r),'sb')
    borders('states','k'); axis(a_xis);
end


scans = table( ...
    ["KOKX","KHGX","KLOT","KNQA","KDAX","KPAH"]',...
    datetime(["24-Oct-2020 3:47", "24-Apr-2020 2:44","16-May-2020 1:55","29-Sep-2020 4:06","25-Apr-2020 4:19","29-Sep-2020 0:17"])',...
    'VariableNames',["station" "time"]);
i_s=2;
scan_radar = scans.station(i_s);
scan_time = scans.time(i_s);

% Large Midwest
scan_radar = [ "KARX" "KDMX" "KDVN" "KEAX" "KFSD" "KGRB" "KGRR" "KILX" "KIND" "KIWX" "KLOT" "KLSX" "KMKX" "KMPX" "KOAX" "KTWX" "KVWX"];
grid_extend = [-96, -86, 38.5, 45];
scan_time =  datetime("16-May-2020") + (2:(.5):10)/24;

% Narrow midwest
scan_radar = ["KMKX","KLOT","KIWX","KDVN","KGRB","KILX","KIND","KARX"]; %,"KGRX","KFSD","KGRR","KMKX","KARX","KDMX","KOAX","KEAX","KLSX","KEAX"];
grid_extend = [-92.5, -85.5, 40, 45];
scan_time =  datetime("16-May-2020") + (2:(.5):10)/24;


% Ithaca
scan_radar = ["KBGM","KTYX","KENX","KCXX", "KCCX"]; %"KBUF",
scan_time =  datetime("29-Sep-2022") + (-3:(.5):9)/24;
grid_extend = [-79, -73, 40.5, 45];


% Cape May
scan_radar = ["KDIX","KDOX"]; %"KBUF",
scan_time =  datetime("9-Oct-2022") + (-1:(.25):10)/24;
grid_extend = [-76, -74, 38, 40];

%% Download
scans_filepath = load_rd(scan_radar, scan_time, download=true);

%% Load
scans_filepath = load_rd(scan_radar, scan_time);

%% Read and process file
rd = read_rd(scans_filepath, ...
    'removeRain','rho',... % "mistnet" slower
    'removeOuliars',4,...
    'smoothing',true,...
    'limitElev',5,...
    'minRange',6000,...
    'maxAlt',3000,...
    'data_res',4);

% save('temp.mat','rd','scans_filepath','scan_radar','scan_time','-v7.3')
% load("temp2_3")

% Load df
r_lim = 150000;
gridsize = num2str(r_lim/1000)+"km";
[df, rd] = load_F_elev("data/" + gridsize + "/", rd);

% save a bit of memory
for i_r=1:numel(rd)
    rd{i_r} = rmfield(rd{i_r},{'vcp','vr'});
end


%% Prepare the grid
g = create_grid(df, rd, grid_extend(1), grid_extend(2), grid_extend(3), grid_extend(4));

[gX,gY,~]=meshgrid(g.x,g.y,g.z);
[gLON,gLAT,gZ]=meshgrid(g.lon,g.lat,g.z);

radar_lon = cellfun(@(x) x.lon, rd(:,1));
radar_lat = cellfun(@(x) x.lat, rd(:,1));
[~,radar_lon_id]=min(abs(g.lon(:)'-radar_lon(:)),[],2);
[~,radar_lat_id]=min(abs(g.lat(:)'-radar_lat(:)),[],2);

%% Compute the F matrix and reshape the data in sparse
[F, data, info] = construct_F(rd, df, g);
szf=size(F);
zmax = quantile(nonzeros(data),.999);

%% Compute geometry stuff and prepare everything for the inversion


% Compute the VPR inversion matrix W
U=repelem(speye(g.sz(3)),prod(g.sz(1:2)),1);
onwater3D=repmat(g.mask_water,1,1,g.sz(3));
Uland = U; Uland(onwater3D(:)==1,:)=0;
W = F*Uland;

% Compute the lapacian smoothing for the VPR
Lh = LaplacianMatrix1D(g.sz(3));
Lh = Lh - diag(sum(Lh,2));
Lh = kron(speye(height(rd)),Lh);

% one VPR per radar
[row, col]=find(W);
[~, ~, col2, ~]=ind2sub(info.sz_rd, row);
Ws = sparse(row,col+g.sz(3)*(col2-1),nonzeros(W),szf(1),g.sz(3)*height(rd));
szWs=size(Ws,2);

L3D = LaplacianMatrix(g, ...
    "water", true,...
    "coeff_horver", 0);


% U forward acconding for Inverse Distance Weight
% compute distance from grid to radar
dgr = sqrt((gX(:,:,1)-permute(g.f_X(radar_lon),[3 2 1])).^2 + (gY(:,:,1)-permute(g.f_Y(radar_lat),[3 2 1])).^2)/1000;
% Create weight of the inverse distance wieghting scheme
p=3;
IDW_w=1./dgr.^p;
IDW_w=IDW_w./sum(IDW_w,3);
IDW_w(IDW_w<.1)=0;
IDW_w=IDW_w./sum(IDW_w,3);
assert(all(sum(IDW_w,3)~=0,'all'))
% figure; tiledlayout('flow','TileSpacing','none','Padding','none')
% for i_s=1:size(IDW_w,3)
%     nexttile; hold on;
%     imagesc(IDW_w(:,:,i_s))
% end


tmp = repmat(permute(IDW_w,[2 1 4 3]),1,1,g.sz(3),1);
tmp = reshape(tmp,[],height(rd));
tmp = repelem(sparse(tmp),1,g.sz(3));
Us = repmat(U,1,height(rd));
Us = Us .* tmp;
Usland = repmat(Uland,1,height(rd));
Usland = Usland .* tmp;
% Us = sparse(Us);

% V = repmat(speye(prod(g.sz(1:2))),g.sz(3),g.sz(3));
clear tmp W col col2 row

%%
DgridzI=nan(g.sz(3)*height(rd),numel(scan_time));
Dinv = nan(g.sz(2),g.sz(1),g.sz(3),numel(scan_time));
Ainv = cell(numel(scan_time),1);
b=cell(numel(scan_time),1);
for i_t=1:numel(scan_time)
    % retrieve id of non-zeros
    id_nnz = find(data(:,i_t));

    % keep only non-zeros0
    D = data(id_nnz,i_t);

    % Inversion of the VPR
    alpha=10;
    Avpr=[Ws(id_nnz,:);alpha*Lh];
    bvpr=[D;sparse(g.sz(3)*height(rd),1)];
    DgridzI(:,i_t) = lsqr(Avpr,bvpr,[],100);

    % Ratio of VPR
%     tmp = reshape(DgridzI(:,i_t),g.sz(3),[]);
%     tmp(tmp<0)=0;
%     tmp = tmp ./ sum(tmp);
%     AgridzI = tmp(:);
%     % V = repmat(speye(prod(g.sz(1:2))),g.sz(3),g.sz(3));
%     Q = repmat(speye(prod(g.sz(1:2))),g.sz(3),g.sz(3)).*(Us*AgridzI)'-speye(prod(g.sz));

    alpha2=1;

    Ainv{i_t}=[F(id_nnz,:); ...
        ...% alpha1*Q;...
        alpha2*L3D];
    b{i_t}=[D; sparse(1*szf(2),1)];

    %     if i_t==1
    %         init = Us*DgridzI(:,i_t);
    %         niter=100;
    %     else
    %         init = Dinvtmp;
    %         niter=100;
    %     end
    %
    %     Dinvtmp = lsqr(Ainv{i_t},b{i_t},[],niter,[],[],init(:));
    %     Dinv(:,:,:,i_t)=permute(reshape(Dinvtmp,g.sz),[2 1 3 4]);

end


%%

i_t=16;

Dinvtmp = lsqr(Ainv{i_t},b{i_t},[],100,[],[],Usland*DgridzI(:,i_t));
Dinv=permute(reshape(Dinvtmp,g.sz),[2 1 3 4]);

% retrieve id of non-zeros
id_nnz = find(data(:,i_t));

% keep only non-zeros
D = full(data(id_nnz,i_t));

% find the optimal max density for plot
zmax = quantile(D(:),.999);

%
[id_range, id_az, id_r, id_df]=ind2sub(info.sz_rd,id_nnz);


df_range = cellfun(@(x) x.range,df,'UniformOutput',false);
tmp = df_range(id_df);
R = cellfun(@(x,y) x(y),tmp, num2cell(id_range));

df_az = cellfun(@(x) x.az,df,'UniformOutput',false);
tmp = df_az(id_df);
A = cellfun(@(x,y) x(y),tmp, num2cell(id_az));

df_e = cellfun(@(x) x.elev,df);
E = df_e(id_df);

[ S, H ] = slant2ground( R, E );
[x1,x2] = pol2cart(cmp2pol(A),R);


figure; hold on;
h = slice(gLON,gLAT,gZ,Dinv,radar_lon,radar_lat,[500]);
set(h,'edgecolor','none')
caxis([0 zmax]); axis tight;

id = R<r_lim & x2<1000&x2>-1000 & H<2000;
scatter3(g.f_lon(g.x(radar_lon_id(id_r(id)))'+x1(id)),g.f_lat(g.y(radar_lat_id(id_r(id)))+x2(id)),H(id),10,D(id),'filled')
id = R<r_lim & x1<1000&x1>-1000 & H<2000;
scatter3(g.f_lon(g.x(radar_lon_id(id_r(id)))'+x1(id)),g.f_lat(g.y(radar_lat_id(id_r(id)))+x2(id)),H(id),10,D(id),'filled')
view(3)
title(datestr(scan_time(i_t)))





%%
Aspeed_prev = cell(numel(scan_time),1);
Aspeed_next = cell(numel(scan_time),1);
Aspeed_coef_prev = cell(numel(scan_time),1);
Aspeed_coef_next = cell(numel(scan_time),1);

uv_u = nan([g.sz numel(scan_time)]);
uv_v = nan([g.sz numel(scan_time)]);

% Extract all speed from radar in a matrix
uv_r = permute(cell2mat(cellfun(@(x) reshape(x.uv,1,1,[],2), rd, 'UniformOutput', false)), [ 1 3 2 4]);  

for i_t=1:numel(scan_time)

    i_t

    % Create speed map for each timestep and altitude. 
    for i_a=1:g.sz(3)
        if height(rd)<=2
            uv_u(:,:,i_a,i_t) = mean(uv_r(:,i_a,i_t,1));
            uv_v(:,:,i_a,i_t) = mean(uv_r(:,i_a,i_t,2));
        else
            tmp = scatteredInterpolant(radar_lon, radar_lat, uv_r(:,i_a,i_t,1), 'natural','nearest');
            uv_u(:,:,i_a,i_t) = tmp({g.lon, g.lat});
            uv_u(:,:,i_a,i_t) = min(uv_u(:,:,i_a),max(uv_r(:,i_a,i_t,1)));
            uv_u(:,:,i_a,i_t) = max(uv_u(:,:,i_a),min(uv_r(:,i_a,i_t,1)));
        
            tmp = scatteredInterpolant(radar_lon, radar_lat,uv_r(:,i_a,i_t,2),'natural');
            uv_v(:,:,i_a,i_t) = tmp({g.lon, g.lat});
            uv_v(:,:,i_a,i_t) = min(uv_v(:,:,i_a),max(uv_r(:,i_a,i_t,2)));
            uv_v(:,:,i_a,i_t) = max(uv_v(:,:,i_a),min(uv_r(:,i_a,i_t,2)));
        end
    end

    if i_t>1

        tmp_u = nanmean(uv_u(:,:,:,i_t-1:i_t),4) ./ g.dx .* seconds((diff(g.time(i_t-1:i_t))));
        tmp_v = nanmean(uv_v(:,:,:,i_t-1:i_t),4) ./ g.dy .* seconds((diff(g.time(i_t-1:i_t))));

        tmp_u = round(reshape(tmp_u,[],g.sz(3)));
        tmp_v = round(reshape(tmp_v,[],g.sz(3)));

        % i_t  <- i_t-1
        [row, col] = ind2sub(g.sz(1:2),(1:prod(g.sz(1:2)))');
        row = row+tmp_u;
        row(row<1)=1;
        row(row>g.sz(1))=g.sz(1);
        col = col+tmp_v;
        col(col<1)=1;
        col(col>g.sz(2))=g.sz(2);
        id = sub2ind(g.sz,row(:),col(:),repelem((1:g.sz(3))',prod(g.sz(1:2))));
        id_full=1:prod(g.sz);
        Aspeed_prev{i_t} = sparse(id_full(~isnan(id)),id(~isnan(id)),1,prod(g.sz),prod(g.sz));
        
        dest_id = sub2ind(g.sz(1:2),row,col);
        Aspeed_coef_prev{i_t} = ones(size(dest_id));
        % Aspeed_coef_prev{i_t}( g.mask_water(:) & ~g.mask_water(dest_id)) = 0.5; % Going from water to land
        % Aspeed_coef_prev{i_t}(~g.mask_water(:) &  g.mask_water(dest_id)) = 0.5; % Going from land to water
        

        % i_t-1 -> i_t
        [row, col] = ind2sub(g.sz(1:2),(1:prod(g.sz(1:2)))');
        row = row-tmp_u;
        row(row<1)=1;
        row(row>g.sz(1))=g.sz(1);
        col = col-tmp_v;
        col(col<1)=1;
        col(col>g.sz(2))=g.sz(2);
        id = sub2ind(g.sz,row(:),col(:),repelem((1:g.sz(3))',prod(g.sz(1:2))));
        id_full=1:prod(g.sz);
        Aspeed_next{i_t-1} = sparse(id_full(~isnan(id)),id(~isnan(id)),1,prod(g.sz),prod(g.sz));

        orig_id = sub2ind(g.sz(1:2),row,col);
        Aspeed_coef_next{i_t-1} = ones(size(orig_id));
        % Aspeed_coef_next{i_t-1}(~g.mask_water(orig_id) &  g.mask_water(:)) = 0.5; % Going from land to water
        % Aspeed_coef_next{i_t-1}( g.mask_water(orig_id) & ~g.mask_water(:)) = 0.5; % Going from water to land
    end
end



%%
Dinvt = nan(prod(g.sz),numel(scan_time));
st = true;
nb_loop=10;
nb_iter=round(linspace(20,100,nb_loop));
wb = waitbar(0,"Inversion: "+ num2str(0) + "/"+num2str(nb_loop));
for i=1:nb_loop
    for i_t=1:numel(scan_time)
        if st
            init = Usland*DgridzI(:,i_t);
            Dinvt(:,i_t) = lsqr(Ainv{i_t},b{i_t},[],20,[],[],init(:));
        else
            if i_t==1
                Dinvt(:,i_t) = lsqr(...
                    [Ainv{i_t}; Aspeed_coef_next{i_t}(:).*Aspeed_next{i_t}],...
                    [b{i_t}; Aspeed_coef_next{i_t}(:).*Dinvt(:,i_t+1)],...
                    [],nb_iter(i),[],[], Dinvt(:,i_t));
            elseif i_t==numel(scan_time)
                Dinvt(:,i_t) = lsqr(...
                    [Ainv{i_t}; Aspeed_coef_prev{i_t}(:).*Aspeed_prev{i_t}],...
                    [b{i_t}; Aspeed_coef_prev{i_t}(:).*Dinvt(:,i_t-1)],...
                    [],nb_iter(i),[],[], Dinvt(:,i_t));
            else
                Dinvt(:,i_t) = lsqr(...
                    [Ainv{i_t}; Aspeed_coef_prev{i_t}(:).*Aspeed_prev{i_t}; Aspeed_coef_next{i_t}(:).*Aspeed_next{i_t}],...
                    [b{i_t}; Aspeed_coef_prev{i_t}(:).*Dinvt(:,i_t-1); Aspeed_coef_next{i_t}(:).*Dinvt(:,i_t+1)],...
                    [],nb_iter(i),[],[], Dinvt(:,i_t));
            end
        end
    end
    st=false;
    waitbar(i/nb_loop,wb,"Inversion: "+ num2str(i) + "/"+num2str(nb_loop));
end
close(wb)
%% 

Dinv = reshape(Dinvt,[g.sz numel(scan_time)]);

rho = squeeze(sum(Dinv,3)); % bird/km^3 -> bird/m^3

tmp = reshape(DgridzI,height(rd),g.sz(3),numel(scan_time));
uv_rg = squeeze(sum(uv_r .* tmp,2) ./ sum(tmp,2));

tmp = reshape(Us*DgridzI,g.sz(1),g.sz(2),g.sz(3),numel(scan_time));
vx = squeeze(sum(tmp.*uv_u,3)./sum(tmp,3));
vy = squeeze(sum(tmp.*uv_v,3)./sum(tmp,3));

[takingoff, landing, entering, leaving, gext, MVT] = sinksource(g,rho/(1000^3),vx,vy);



%%

dlq=10;

figure('position',[0 0 800 800]); hold on;
im = imagesc(g.lon,g.lat,rho(:,:,i_t)');
axis equal tight;  a_xis=axis;
% qm = quiver(g.lon(1:dlq:end,1:dlq:end),g.lat(1:dlq:end,1:dlq:end),vx(1:dlq:end,1:dlq:end,i_t),vy(1:dlq:end,1:dlq:end,i_t),1,'k');
plot(radar_lon,radar_lat,'or');
qr = quiver(radar_lon,radar_lat,uv_rg(:,i_t,1),uv_rg(:,i_t,2),.25,'r');

borders('states','w'); axis(a_xis);
% set(gca,'ColorScale','log')
clim([0 3*zmax]);
set(gca,"XColor","w","YColor","w")

for i_t=1:numel(scan_time)
    title(datestr(scan_time(i_t)),FontSize=20,Color="w")
    im.CData = rho(:,:,i_t)';
    %qm.UData = vx(1:dlq:end,1:dlq:end,i_t);
    %qm.VData = vy(1:dlq:end,1:dlq:end,i_t);
    qr.UData = uv_rg(:,i_t,1);
    qr.VData = uv_rg(:,i_t,2);
    pause(.2)
    % exportgraphics(gcf,"test3.gif",'Append',true);
    exportgraphics(gcf,strjoin(scan_radar,"_")+"_"+datestr(mean(scan_time),'YYYYMMDD')+".gif",'Append',true,'BackgroundColor',[48 48 48]/256);
end








%%
c_lim = 100;
clmap = brewermap([],'Spectral');

Fsm = takingoff+landing;
for i_a=1:g.sz(3)
    Fsm(:,:,i_a) = imgaussfilt(Fsm(:,:,i_a),2,'FilterSize',5); 
end
takingoff_s = Fsm;
takingoff_s(takingoff_s<0) = 0;
landing_s = Fsm;
landing_s(landing_s>0) = 0;

figure('position',[0 0 1200 600]); tiledlayout('flow','TileSpacing','tight','Padding','tight');
colormap(clmap)

nexttile; hold on;
im1 = imagesc(g.lon,g.lat,takingoff_s(:,:,1)');
axis equal tight;  a_xis=axis;
borders('states','k'); axis(a_xis); 
% set(gca,'ColorScale','log')
clim([-c_lim c_lim]);
title("Taking-off "+datestr(scan_time(i_t)),'FontSize',24,'Color','w'); box off
xticklabels(''); yticklabels('');

nexttile; hold on;
im2 = imagesc(g.lon,g.lat,landing_s(:,:,1)');
axis equal tight;  a_xis=axis;
borders('states','k'); axis(a_xis); 
% set(gca,'ColorScale','log')
clim([-c_lim c_lim]);
title("Landing "+datestr(scan_time(i_t)),'FontSize',24,'Color','w'); box off
xticklabels(''); yticklabels('');


for i_t=1:numel(scan_time)-1
    title(datestr(scan_time(i_t)))
    im1.CData = sum(takingoff_s(:,:,1:i_t),3)';
    im2.CData = sum(landing_s(:,:,1:i_t),3)';
    pause(.2)
   % exportgraphics(gcf,"takeoff-landing_"+strjoin(scan_radar,"_")+"_"+datestr(mean(scan_time),'YYYYMMDD')+".gif",'Append',true,'BackgroundColor',[48 48 48]/256);
end
















%% Inversion of the vertical profile
[ ~, H ] = slant2ground( R(:), E(:) );
[~,~,gr] = histcounts(H,[g.z(1)-diff(g.z(1:2)) ;  squeeze(g.z)]+diff(g.z(1:2))/2);
Dh=nan(numel(g.z), height(rd));
for i_r = 1:height(rd)
    [G,ID] = findgroups(gr(id_r==i_r));
    Dh(ID,i_r)=splitapply(@nanmean, D(id_r==i_r), G);
end
figure; plot(Dh,squeeze(g.z),'-o')



% DgridzI = W\D2;
% figure; plot(DgridzI,squeeze(g.z)','-o')


alpha=10;
Ainv=[Ws(id_nnz,:);alpha*Lh];
b=[D;sparse(g.sz(3)*height(rd),1)];

% DgridzI = Ainv\b;
DgridzI = lsqr(Ainv,b,[],100,[],[],Dh(:));

figure; hold on
plot(Dh,squeeze(g.z)/1000,'-o')
set(gca,'ColorOrderIndex',1)
plot(reshape(DgridzI,g.sz(3),[]),squeeze(g.z)'/1000,'--o')
legend(scan_radar)



Dinit = permute(reshape(Us*DgridzI,g.sz),[2 1 3]);
figure;
h = slice(gLON,gLAT,gZ,Dinit,[radar_lon],[radar_lat],[100]);
set(h,'edgecolor','none')
axis tight;caxis([0 zmax]);

% figure; imagesc(sum(Dinit,3))

%% Compute VPR ratio
tmp = reshape(DgridzI,g.sz(3),[]);
tmp = tmp ./ sum(tmp);
AgridzI = tmp(:);

figure; hold on
plot(reshape(AgridzI,g.sz(3),[]),squeeze(g.z)'/1000,'--o')

Ainit = permute(reshape(Us*AgridzI,g.sz),[2 1 3]);
figure;
h = slice(gLON,gLAT,gZ,Ainit,[radar_lon],radar_lat,[100]);
set(h,'edgecolor','none'); axis tight;

Q = repmat(speye(prod(g.sz(1:2))),g.sz(3),g.sz(3)).*(Us*AgridzI)'-speye(prod(g.sz));


%%

alpha1=1;%.25*height(scans);
alpha2=2;%.5*height(scans);

% Ainv=[Fsp2 sparse(szf(1),g.sz(3)); ...
%     alpha1*speye(szf(2)) -alpha1*U; ...
%     alpha2*L3D sparse(szf(2), g.sz(3))];
% b=[D2(:); sparse(2*szf(2),1)];

% one vpr per radar
% Ainv=[Fsp2 sparse(szf(1),szWs); ...
%     sparse(szf(1),szf(2),1) Ws; ...
%     alpha1*speye(szf(2)) -alpha1*Us;...
%     alpha2*L3D sparse(szf(2), szWs)];
% b=[D2(:); D2(:); sparse(2*szf(2),1)];

Ainv=[F(id_nnz,:); ...
    ...%alpha1*Q;...
    alpha2*L3Dwater];
b=[D; sparse(szf(2),1)];

Dinv = lsqr(Ainv,b,[],100,[],[],Us*DgridzI);
Dinv = permute(reshape(Dinv,g.sz),[2 1 3]);
% Dgridz = reshape(Dinvtmp((szf(2)+1):end),g.sz(3),[]);


figure;
h = slice(gLON,gLAT,gZ,Dinv,radar_lon,radar_lat,[300]);
set(h,'edgecolor','none')
caxis([0 zmax]); %set(gca,'ColorScale','log')
axis tight;title(datestr(scan_time(i_t)))
% Ax=gca; Ax.XColor = 'w'; Ax.YColor = 'W'; Ax.ZColor = 'W'; Ax.Color = [50 50 50]/255;
% exportgraphics(gcf,'test-6.png','BackgroundColor',[50 50 50]/255)



figure; hold on;
imagesc(gLON(1,:,1),gLAT(:,1,1),(sum(Dinv,3)))
axis equal tight;  a_xis=axis;
borders('states','w'); axis(a_xis); set(gca,'ColorScale','log')
% Ax=gca; Ax.XColor = 'w'; Ax.YColor = 'W'; Ax.ZColor = 'W'; Ax.Color = [50 50 50]/255;
% exportgraphics(gcf,'test-7.png','BackgroundColor',[50 50 50]/255)


%%

figure('position',[0 0 1400 800]); tiledlayout(3,5,'TileSpacing','tight','Padding','tight')
z=squeeze(g.z)/1000;
X2=gX/1000; Y2=gY/1000; Z2=gZ/1000;
[ xrd, yrd, zrd ] = radar2xyz( R2, A2, E2 );

nexttile; hold on; box on; grid on; axis equal
err1 = Fsp2 * Dinvtmp(1:szf(2)) - D2(:);%sqrt(mean((Fsp2 * Dinv(:) - D2(:)).^2));
plot(D2(:), Fsp2 * Dinvtmp(1:szf(2)), '.k'); plot([0 max(D2(:))],[0 max(D2(:))],'r','linewidth',2)
% xline(mean(D2(:)),'w'); yline(mean(Fsp2 * Dinv(:)),'w')
% xline(quantile(D2(:),[.5 .9 .99]),'w'); yline(quantile(Fsp2 * Dinv(:),[.5 .9 .99]),'w')
% histogram(err1); set(gca,'yscale','log')
title("||Fd_{true}-d_{obs}|| = "+round(sqrt(sum(err1.^2))/1000) +"K");
xlabel('d_{obs}'); ylabel('Fd_{true}');

nexttile; hold on; box on; grid on; colorbar;
title("||Fd_{true}-d_{obs}||" +" avg over az.");
eu = unique(E2);
for i_elev=1:numel(eu)
    id = find(E2==eu(i_elev)&R2<80000);
    R2b = round(R2(id),-3);
    [G,ID] = findgroups(R2b);
    [ s, h ] = slant2ground( ID, eu(i_elev) );
    tmp = splitapply(@mean, err1(id), G);
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
err2 = Dinvtmp(1:szf(2)) - Us*Dinvtmp((szf(2)+1):end); %sqrt(mean((Dinv(:) - U*Dgridz(:)).^2));
%plot(U*Dgridz(:), Dinvtmp(1:szf(2)), '.k'); plot([0 max(D2(:))],[0 max(D2(:))],'color',[.7 .7 .7])
% histogram(err2); set(gca,'yscale','log')
tmp = mean(reshape(err2,g.sz),3);
% tmp(onwater)=nan;
imagesc(tmp,'alphadata',~isnan(tmp))
% boxplot(reshape(Dinv,[],size(Dinv,3)),'PlotStyle','compact','Orientation','horizontal')%,'Whisker',10)
% plot(Dgridz,1:31,'-or','linewidth',2);
title("||d_{true}-Ud_{VPR}|| = "+round(sqrt(sum(err2.^2))/1000) +"K"+ " | \alpha_1="+alpha1);
% xlabel('Density'); ylabel('Z'); xlim([0 zmax]);
xlabel('X'); ylabel('Y'); caxis([-100 100]); yticklabels(''); xticklabels(''); axis equal tight;

nexttile; hold on; box on; grid on;
err3 = L3D*Dinvtmp(1:szf(2));% sqrt(mean((L3D*Dinv(:)).^2));
% plot(U*Dgridz(:), Dinv(:), '.k'); plot([0 max(D2(:))],[0 max(D2(:))],'color',[.7 .7 .7])
histogram(err3); set(gca,'yscale','log'); xline(0,'r','linewidth',2)
title("||L_{3D}d_{true}|| = "+round(sqrt(sum(err3.^2))/1000) +"K"+ " | \alpha_2="+alpha2);
xlabel('Density'); ylabel('Histogram'); yticklabels('')

nexttile([2 1]); hold on; box on; grid on; title('VPR')
tmp = reshape(Dinv,[],size(Dinv,3));
id = randsample(size(tmp,1),50);
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
scatter(xrd(id)/1000, zrd(id)/1000,20,D2(id),'filled')
caxis([0 zmax]); xlim([-80 80]); ylim([0 3]); %xticklabels(''); yticklabels('')
ylabel('Z'); xlabel('X');xticklabels(''); yticklabels('')

nexttile([1 2]); hold on; box on; grid on; title('X=0')
imagesc(squeeze(Y2(:,1,1))',squeeze(Z2(1,1,:))',squeeze(Dinv(round(end/2),:,:))')
id = R2<r_lim & zrd<3000 & abs(xrd)<1000;
scatter(yrd(id)/1000, zrd(id)/1000,20,D2(id),'filled')
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
szf=size(Fspns);

alpha=8;
Ainv=[Fspns sparse(szf(1),size(W,2)); ...
    sparse(numel(id_rd),szf(2)) W; ...
    alpha*speye(szf(2)) -alpha*U];
b=[D(id_rd);D(id_rd); sparse(szf(2),1)];

tic; Dinvtmp = lsqr(Ainv,b,[],20); toc;
Dinv=Dinvtmp(1:size(Fspns,2));
Dgridz = Dinvtmp(size(Fspns,2)+1:end);

figure; hold on; box on; grid on;
plot(Dh,squeeze(z),'-o','color',[.7 .7 .7])
plot(mean(Dh,2),squeeze(z),'-ok')
plot(Dgridz,squeeze(z),'-or')
xlabel('DZ (?)'); ylabel('Elevation');
legend('Mean of d_{ob}','Inversion','Inversion')

% tmp=repmat(permute(Dgridz,[2 3 1]),g.sz(1),g.sz(2));
tmp=reshape(Dinv,g.sz);
tmp = permute(tmp,[2 1 3]);
% tmp(tmp==0)=nan;
figure; tiledlayout('flow','Padding','none');
h=slice(gLON,gLAT,gZ,tmp,g.lon([stations(:).lon_id]),g.lat([stations(:).lat_id]),100);
set(h,'edgecolor','none')
axis vis3d; caxis([0 30]);

figure; imagesc(g.lon,g.lat,sum(tmp,3,'omitnan')); caxis([300 500]);



%% Solve alltogether
alpha1=4; alpha2=5;

szf=size(Fspns);

Ainv=[Fspns sparse(szf(1),size(W,2)); ...
    sparse(numel(id_rd),szf(2)) W; ...
    alpha1*speye(szf(2)) -alpha1*U;...
    alpha2*L3D sparse(szf(2),size(W,2))];
b=[D(id_rd);D(id_rd); sparse(szf(2),1); sparse(szf(2),1)];

tic;Dinvtmp = lsqr(Ainv,b,[],100);toc
Dinv=Dinvtmp(1:size(Fspns,2));
Dgridz = Dinvtmp(size(Fspns,2)+1:end);



tmp=reshape(Dinv,g.sz);
tmp = permute(tmp,[2 1 3]);
figure; tiledlayout('flow','Padding','none');
h=slice(gLON,gLAT,gZ,tmp,g.lon([stations(:).lon_id]),g.lat([stations(:).lat_id]),100);
set(h,'edgecolor','none')
axis vis3d tight; caxis([0 30]);

figure; imagesc(g.lon,g.lat,sum(tmp,3,'omitnan')); caxis([300 500]);


%% polynomial vertical profile

% same as before
p_n=4;

V = (gZ(:)/1000).^(0:p_n);
WV = (Fsp(id_rd,:)'./sum(Fsp(id_rd,:),2)')';
p = (WV*V) \ D(id_rd);

figure; hold on; box on; grid on;
plot(Dh,squeeze(z),'-o','color',[.7 .7 .7])
plot(mean(Dh,2),squeeze(z),'-ok')
plot((squeeze(z/1000).^(0:p_n)) * p,squeeze(z),'-or')

%% Space time polynomial vertical profile
p_n=[0 0 4];

% V = [ones(prod(g.sz),1) (gY(:)./1000000).^(1:p_n(1)) (gX(:)./1000000).^(1:p_n(2))  (gZ(:)./1000).^(1:p_n(3))];
V = (gZ(:)/1000).^(0:p_n(3));
WV = (Fsp(id_rd,:)'./sum(Fsp(id_rd,:),2)')';
% p = (WV*V) \ D(id_rd);
p = lsqr(WV*V,D(id_rd));

figure; hold on; box on; grid on;
plot(Dh,squeeze(z),'-o','color',[.7 .7 .7])
plot(mean(Dh,2),squeeze(z),'-ok')
plot((squeeze(z/1000).^(0:p_n(3))) * [p(1); p((sum(p_n(1:2))+2):end)],squeeze(z),'-or')

M = permute(reshape(V*p,[g.sz(2) g.sz(1) g.sz(3)]),[2 1 3]);

figure; plot(squeeze(mean(M,[1 2])))
nb_prof=3;
figure; tiledlayout(nb_prof,nb_prof,'TileSpacing','tight','Padding','none')
for i_x=1:nb_prof
    for i_y=1:nb_prof
        nexttile;
        plot(squeeze(M(round(numel(g.x)/(nb_prof+1)*i_x),round(numel(g.y)/(nb_prof+1)*i_y),:)),squeeze(z),'-ok' );
        xlim([0 25])
    end
end

% Solve alltogether
alpha=1; alpha2=5;
Ainv=[Fspns sparse(szf(1),size(V,2)); ...
    sparse(numel(id_rd),szf(2)) WV*V; ...
    alpha*speye(szf(2)) -alpha*V(id_gr,:)];%;...
%alpha2*L3D(id_gr,id_gr) sparse(numel(id_gr),size(V,2))];
b=[D(id_rd);D(id_rd); zeros(szf(2),1)];%; zeros(numel(id_gr),1)];

tic;Dinvtmp = lsqr(Ainv,b,[],100);toc
Dinv=Dinvtmp(1:size(Fspns,2));
Dparm = Dinvtmp(size(Fspns,2)+1:end);


tmp = repmat(permute(Dgridz,[2 3 1]),g.sz(1),g.sz(2));
tmp(id_gr)=Dinv;
tmp = permute(tmp,[2 1 3]);
tmp(tmp==0)=nan;
figure; tiledlayout('flow','Padding','none');
h=slice(gLON,gLAT,gZ,tmp,g.lon([stations(:).lon_id]),g.lat([stations(:).lat_id]),100);
set(h,'edgecolor','none')
axis vis3d tight; caxis([0 30]);

figure; imagesc(g.lon,g.lat,sum(tmp,3,'omitnan'))