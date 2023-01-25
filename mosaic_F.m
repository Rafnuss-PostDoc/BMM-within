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
grid_extend = [-97, -85, 38, 45];
scan_time =  datetime("16-May-2020") + (1.5:(.5):12)/24;

% Narrow midwest
scan_radar = ["KMKX","KLOT","KIWX","KDVN","KGRB","KILX","KIND","KARX"]; %,"KGRX","KFSD","KGRR","KMKX","KARX","KDMX","KOAX","KEAX","KLSX","KEAX"];
grid_extend = [-91.5, -85.5, 40, 45];
scan_time =  datetime("16-May-2020") + (2:(.5):11)/24;

% very narrow midwest
scan_radar = ["KMKX","KLOT","KDVN","KARX"]; %,"KGRX","KFSD","KGRR","KMKX","KARX","KDMX","KOAX","KEAX","KLSX","KEAX"];
grid_extend = [-92.5, -87, 40.5, 44];
scan_time =  datetime("16-May-2020") + (1.5:(.5):11)/24;


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
% Read and process file
scans_filepath = load_rd(scan_radar, scan_time);
rd = read_rd(scans_filepath, ...
    'removeRain','rho',... % "mistnet" slower
    'limitElev',5,...
    'minRange',6000,...
    'maxAlt',3000,...
    'removeOuliars',4,...
    'shape',false,...
    'smoothing',true,...
    'data_res',4);

% Load df
r_lim = 150000;
[df, rd] = load_F_elev("data/" + num2str(r_lim/1000)+"km/", rd);

% save a bit of memory
for i_r=1:numel(rd)
    rd{i_r} = rmfield(rd{i_r},{'vcp','vr'});
end

% save('temp.mat','df','rd','scans_filepath','scan_radar','scan_time','grid_extend','-v7.3')
% load("temp")


%% Prepare the grid
g = create_grid(df, rd, grid_extend(1), grid_extend(2), grid_extend(3), grid_extend(4),2000);

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

% Compute the W matric for the VPR inversion: only consider above water
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
if false
    figure; tiledlayout('flow','TileSpacing','none','Padding','none')
    for i_s=1:size(IDW_w,3)
        nexttile; hold on;
        imagesc(IDW_w(:,:,i_s))
    end
end

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

%% Build the main system (without advection)

DgridzI=nan(g.sz(3)*height(rd),numel(scan_time));
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
end


%% Check what the vpr look like
% use as inital model

if false

    figure; tiledlayout('flow','TileSpacing','tight','Padding','none')
    tmp = reshape(DgridzI,g.sz(3),height(rd),[]);
    for i_s=1:size(IDW_w,3)
        nexttile; plot(squeeze(tmp(:,i_s,:)),squeeze(g.z))
        xlim([0 max(DgridzI(:))]); title(rd{i_s,1}.station)
    end

    % Option: All
%     Dinv = nan(g.sz(2),g.sz(1),g.sz(3),numel(scan_time));
%     for i_t=1:numel(scan_time)
%         if i_t==1
%             init = Usland*DgridzI(:,i_t);
%             niter=100;
%         else
%             init = Dinvtmp;
%             niter=100;
%         end
%         Dinvtmp = lsqr(Ainv{i_t},b{i_t},[],niter,[],[],init(:));
%         Dinv(:,:,:,i_t)=permute(reshape(Dinvtmp,g.sz),[2 1 3 4]);
%     end
    % Option: Just one

    i_t=15;
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

    figure('position',[0 0 800 800]); hold on; hold on;

    % h = slice(gLON,gLAT,gZ,permute(reshape(Usland*DgridzI(:,i_t),g.sz),[2 1 3 4]),radar_lon,radar_lat,[500]); 
     h = slice(gLON,gLAT,gZ,Dinv,radar_lon,radar_lat,[500]);

    set(h,'edgecolor','none')
    clim([0 zmax]); axis tight;
    
    id = R<max(df{1}.x) & x2<1000&x2>-1000 & H<2000;
    scatter3(g.f_lon(g.x(radar_lon_id(id_r(id)))'+x1(id)),g.f_lat(g.y(radar_lat_id(id_r(id)))+x2(id)),H(id),50,D(id),'filled')
    id = R<max(df{1}.x) & x1<1000&x1>-1000 & H<2000;
    scatter3(g.f_lon(g.x(radar_lon_id(id_r(id)))'+x1(id)),g.f_lat(g.y(radar_lat_id(id_r(id)))+x2(id)),H(id),50,D(id),'filled')
    view(3)
    title(datestr(scan_time(i_t)))

end



%% 

uv_u = nan([g.sz numel(scan_time)]);
uv_v = nan([g.sz numel(scan_time)]);

% Extract all speed from radar in a matrix
uv_r = permute(cell2mat(cellfun(@(x) reshape(x.uv,1,1,[],2), rd, 'UniformOutput', false)), [ 1 3 2 4]);  

for i_t=1:numel(scan_time)
    % Create speed map for each timestep and altitude. 
    for i_a=1:g.sz(3)
        if height(rd)<=2
            uv_u(:,:,i_a,i_t) = mean(uv_r(:,i_a,i_t,1));
            uv_v(:,:,i_a,i_t) = mean(uv_r(:,i_a,i_t,2));
        else
            % interpolate speed over space
            tmp = scatteredInterpolant(radar_lon, radar_lat, uv_r(:,i_a,i_t,1), 'natural');
            uv_u(:,:,i_a,i_t) = tmp({g.lon, g.lat});
            uv_u(:,:,i_a,i_t) = min(uv_u(:,:,i_a,i_t),max(uv_r(:,i_a,i_t,1)));
            uv_u(:,:,i_a,i_t) = max(uv_u(:,:,i_a,i_t),min(uv_r(:,i_a,i_t,1)));
        
            tmp = scatteredInterpolant(radar_lon, radar_lat,uv_r(:,i_a,i_t,2),'natural');
            uv_v(:,:,i_a,i_t) = tmp({g.lon, g.lat});
            uv_v(:,:,i_a,i_t) = min(uv_v(:,:,i_a,i_t),max(uv_r(:,i_a,i_t,2)));
            uv_v(:,:,i_a,i_t) = max(uv_v(:,:,i_a,i_t),min(uv_r(:,i_a,i_t,2)));
        end
    end
end


if false
    figure; tiledlayout('flow','TileSpacing','tight','Padding','none')
    uv_r = permute(cell2mat(cellfun(@(x) reshape(x.uv,1,1,[],2), rd, 'UniformOutput', false)), [ 1 3 2 4]);  
    for i_s=1:size(IDW_w,3)
        nexttile; quiver(datenum(repmat(g.time',1,g.sz(3)))*24,repmat(squeeze(g.z)'/100,width(rd),1),squeeze(uv_r(i_s,:,:,1))',squeeze(uv_r(i_s,:,:,2))',1,'k')
       title(rd{i_s,1}.station)
    end

    figure; tiledlayout('flow','TileSpacing','tight','Padding','none')
    tmp = reshape(DgridzI,height(rd),g.sz(3),numel(scan_time));
    uv_rg = squeeze(sum(uv_r(:,1:g.sz(3),:,:) .* tmp,2) ./ sum(tmp,2));
    dlq=25;
    [X,Y]=meshgrid(g.lon(1:dlq:end,1:dlq:end),g.lat(1:dlq:end,1:dlq:end));
    for i_t=1:5:numel(g.time)
        for i_a=1:8:g.sz(3)
           nexttile; hold on; 
            borders('states','k'); 
           quiver(X,Y,uv_u(1:dlq:end,1:dlq:end,i_a,i_t)',uv_v(1:dlq:end,1:dlq:end,i_a,i_t)',1,'k');
           plot(radar_lon,radar_lat,'or');
           qr = quiver(radar_lon,radar_lat,uv_rg(:,i_t,1),uv_rg(:,i_t,2),.25,'r');
           qr = quiver(radar_lon,radar_lat,uv_r(:,i_a,i_t,1),uv_r(:,i_a,i_t,2),.25,'k');
           title(num2str(g.z(i_a)) + " - " + datestr(g.time(i_t))); axis equal tight;axis(grid_extend);
        end
    end
end

%% Build the advection regularization matrix

Aspeed_prev = cell(numel(scan_time),1);
Aspeed_next = cell(numel(scan_time),1);
Aspeed_coef_prev = cell(numel(scan_time),1);
Aspeed_coef_next = cell(numel(scan_time),1);

coeff_flux = interp1(0:0.1:1,[.2 0.9 1.2 1.6 1.8 1.9 1.8 1.6 1.2 0.9 .2],linspace(0,1,numel(scan_time)),'cubic');

for i_t=2:numel(scan_time)

        % Convert speed into distance in the grid. 
        tmp_u = nanmean(uv_u(:,:,:,i_t-1:i_t),4) ./ g.dx .* seconds((diff(g.time(i_t-1:i_t))));
        tmp_v = nanmean(uv_v(:,:,:,i_t-1:i_t),4) ./ g.dy .* seconds((diff(g.time(i_t-1:i_t))));
        tmp_u = round(reshape(tmp_u,[],g.sz(3)));
        tmp_v = round(reshape(tmp_v,[],g.sz(3)));

        % build Aspeed_prev corresponding to the transition backward
        % i_t  <- i_t-1
        [row, col] = ind2sub(g.sz(1:2),(1:prod(g.sz(1:2)))'); % find the index in x and y of the entire grid
        % add the transition 
        row = row+tmp_u; 
        col = col+tmp_v;

        % for all transition leading to oustide the grid, use the nearest
        % grid cell.
        id_edge = row<1|row>g.sz(1)|col<1|col>g.sz(2);
        row(row<1)=1;
        row(row>g.sz(1))=g.sz(1);
        col(col<1)=1;
        col(col>g.sz(2))=g.sz(2);

        % Build the transition matrix from origin (id_full) and destination
        % in 3d (id)
        id = sub2ind(g.sz,row(:),col(:),repelem((1:g.sz(3))',prod(g.sz(1:2))));
        id_full=1:prod(g.sz);
        Aspeed_prev{i_t} = sparse(id_full(~isnan(id)),id(~isnan(id)),1,prod(g.sz),prod(g.sz));
        id_same = id(:)==id_full(:);
        
        % compute coefficient of the transition. And adjust if necessary
        dest_id = sub2ind(g.sz(1:2),row,col);
        Aspeed_coef_prev{i_t} = coeff_flux(i_t).*ones(size(dest_id));
        % Aspeed_coef_prev{i_t}( g.mask_water(:) & ~g.mask_water(dest_id)) = 0.5; % Going from water to land
        % Aspeed_coef_prev{i_t}(~g.mask_water(:) &  g.mask_water(dest_id)) = 0.5; % Going from land to water
        %Aspeed_coef_prev{i_t}(g.mask_water(:) &  g.mask_water(dest_id)) = 2; % Going from water to water
        Aspeed_coef_prev{i_t}(id_edge) = 0; % contrains on the edge of the domain
        Aspeed_coef_prev{i_t}(id_same) = 0;
        

        % i_t-1 -> i_t
        [row, col] = ind2sub(g.sz(1:2),(1:prod(g.sz(1:2)))');
        row = row-tmp_u;
        col = col-tmp_v;
        id_edge = row<1|row>g.sz(1)|col<1|col>g.sz(2);
        row(row<1)=1;
        row(row>g.sz(1))=g.sz(1);
        col(col<1)=1;
        col(col>g.sz(2))=g.sz(2);
        id = sub2ind(g.sz,row(:),col(:),repelem((1:g.sz(3))',prod(g.sz(1:2))));
        id_full=1:prod(g.sz);
        id_same = id(:)==id_full(:);
        Aspeed_next{i_t-1} = sparse(id_full(~isnan(id)),id(~isnan(id)),1,prod(g.sz),prod(g.sz));
        orig_id = sub2ind(g.sz(1:2),row,col);
        Aspeed_coef_next{i_t-1} = coeff_flux(i_t).*ones(size(orig_id));
        % Aspeed_coef_next{i_t-1}(~g.mask_water(orig_id) &  g.mask_water(:)) = 0.5; % Going from land to water
        % Aspeed_coef_next{i_t-1}( g.mask_water(orig_id) & ~g.mask_water(:)) = 0.5; % Going from water to land
        %Aspeed_coef_next{i_t-1}(g.mask_water(:) &  g.mask_water(dest_id)) = 2; % Going from water to water
        Aspeed_coef_next{i_t-1}(id_edge) = 0;
        Aspeed_coef_next{i_t}(id_same) = 0;
end



%%
if ~(exist("Dinvt","var") && all(size(Dinvt)==[prod(g.sz),numel(scan_time)]))
    Dinvt = nan(prod(g.sz),numel(scan_time));
    st = true;
end
nb_loop=10;
nb_iter=round(linspace(20,100,nb_loop));
wb = waitbar(0,"Inversion: "+ num2str(0) + "/"+num2str(nb_loop));
alphat=1;
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
                    [Ainv{i_t}; alphat.*Aspeed_coef_prev{i_t}(:).*Aspeed_prev{i_t}],...
                    [b{i_t}; alphat.*Aspeed_coef_prev{i_t}(:).*Dinvt(:,i_t-1)],...
                    [],nb_iter(i),[],[], Dinvt(:,i_t));
            else
                Dinvt(:,i_t) = lsqr(...
                    [Ainv{i_t}; alphat.*Aspeed_coef_prev{i_t}(:).*Aspeed_prev{i_t}; Aspeed_coef_next{i_t}(:).*Aspeed_next{i_t}],...
                    [b{i_t}; alphat.*Aspeed_coef_prev{i_t}(:).*Dinvt(:,i_t-1); Aspeed_coef_next{i_t}(:).*Dinvt(:,i_t+1)],...
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


% tmp = reshape(Us*DgridzI,g.sz(1),g.sz(2),g.sz(3),numel(scan_time));
% vx = squeeze(sum(tmp.*uv_u,3)./sum(tmp,3));
% vy = squeeze(sum(tmp.*uv_v,3)./sum(tmp,3));
% 
% [takingoff, landing, entering, leaving, gext, MVT] = sinksource(g,rho/(1000^3),vx,vy);



figure('position',[0 0 800 800]); tiledlayout(1,1,'TileSpacing','none','Padding','tight'); nexttile; hold on;
im = imagesc(g.lon,g.lat,rho(:,:,i_t)');
axis equal tight;  a_xis=axis;
% qm = quiver(g.lon(1:dlq:end,1:dlq:end),g.lat(1:dlq:end,1:dlq:end),vx(1:dlq:end,1:dlq:end,i_t),vy(1:dlq:end,1:dlq:end,i_t),1,'k');
plot(radar_lon,radar_lat,'or');
qr = quiver(radar_lon,radar_lat,uv_rg(:,i_t,1),uv_rg(:,i_t,2),.25,'r');

borders('states','w'); axis(a_xis);
% set(gca,'ColorScale','log')
clim([0 9*zmax]); colormap("jet")
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
    % exportgraphics(gcf,strjoin(scan_radar,"_")+"_"+datestr(mean(scan_time),'YYYYMMDD')+".gif",'Append',true,'BackgroundColor',[48 48 48]/256);
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



figure;
sum(takingoff_s(:,:,1:i_t),3)'


%%

thr_im=quantile(rho(:),.75) ;

figure('position',[0 0 800 1000]);  tiledlayout(4,1,'TileSpacing','none','Padding','tight')


for i_t=1:5:18
    nexttile; hold on;
    tmp = rho(:,:,i_t)';
    scatter(radar_lon,radar_lat, 50,"w",'filled','MarkerEdgeColor','k');
    imagesc(g.lon,g.lat,rho(:,:,i_t)','Alphadata', (tmp-thr_im)/5000.*(tmp>thr_im));
    axis equal tight; a_xis=axis;
    clim([thr_im quantile(rho(:),.999)]); colormap("jet")
    set(gca,"XColor","none","YColor","none")
    plot_google_map(MapType='satellite',Refresh=0, AutoAxis=0); % [lonVec, latVec, imag] = 
    axis([-97 -86 38.5 45]);
 end



%%
figure(Position=[0 0 1400 550]);  hold on;
% 
Ax=gca; % Ax.XColor = 'w'; Ax.YColor = 'W'; Ax.ZColor = 'W';  
Ax.Color = [230 230 230]/255;
Ax.DataAspectRatio = [1 1 800];
Ax.XTick=[];Ax.YTick=[];
zlim([0 2000])
col=colormap("jet");
axis(a_xis); 
[lonVec, latVec, imag] = plot_google_map(MapType='satellite');
Fdem = griddedInterpolant({g.lon,g.lat},imgaussfilt(g.dem,'FilterSize',3)','nearest');
warp(lonVec,latVec,Fdem({lonVec, latVec})',imag);
axis(a_xis); Ax.YDir="normal";
s=scatter3(radar_lon,radar_lat, cellfun(@(x) x.height, rd(:,1)), 100,"w",'filled','MarkerEdgeColor','k');

% [caz,cel] = view();
view([16,13])
view([-27,23])
c_lim=700;
c_lim2 = linspace(c_lim,2000,5);

for i_t=1:numel(g.time)
    tmp = smooth3(Dinv(:,:,:,i_t));
    for i=1:numel(c_lim2)
        p(i)=patch(isosurface(permute(gLON,[2 1 3]), permute(gLAT,[2 1 3]),permute( g.dem+gZ,[2 1 3]),tmp,c_lim2(i)),...
            'EdgeColor','none','FaceColor',col(round(i*255/numel(c_lim2)),:), "facealpha",0.4);
    end
    % pause(.2);
    exportgraphics(gcf,'testAnimated2.gif','Append',true);
    delete(p)
end

clim([c_lim max(Dfor2(:))])
p2=patch(isocaps(gX(:,end/2:end,:),gY(:,end/2:end,:),gZ(:,end/2:end,:),Dfor2,c_lim),...
    'FaceColor','interp','EdgeColor','none', 'FaceAlpha',0.8);

view(3)
box on; grid on;
xticks(g.x(1:100:end)); yticks(g.y(1:100:end));zticks(0:1000:3000)
xticklabels([-50 0 50]); yticklabels([-50 0 50]);
Ax=gca; % Ax.XColor = 'w'; Ax.YColor = 'W'; Ax.ZColor = 'W';  
Ax.Color = [230 230 230]/255;
Ax.DataAspectRatio = [1 1 1/15];
xlim([0 100000]); ylim([0 100000]); zlim([-100 2000])

view(-80,9)














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