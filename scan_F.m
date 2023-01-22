
addpath('/Users/raphael/Documents/GitHub/wsrlib');
wsrlib_setup();
addpath(genpath('functions/'))


% Florida
scan_radar = ["KBYX"];
scan_day = datetime(["24-Apr-2019", "19-Apr-2019", "14-Apr-2019","13-Apr-2019",  ...
    "7-May-2018", "23-Apr-2018", "22-Apr-2018", "14-Apr-2018", "7-Apr-2018",...
    "12-May-2017","10-May-2017", "30-Apr-2017", "27-Apr-2017", "25-Apr-2017"]);
i_d = 11;
scan_time =  scan_day(i_d) +  (1:(.25):10)/24;

% Cape May
% scan_radar = ["KDOX"]; %"KBUF",
% scan_time =  datetime("9-Oct-2022") + (-1:(.25):10)/24;

%% Download
scans_filepath = load_rd(scan_radar, scan_time, download=true);
 
%% Load file
scans_filepath = load_rd(scan_radar, scan_time);

%% Read and process file
rd = read_rd(scans_filepath, ...
    'removeRain', 'rho',... % "mistnet" slower
    'removeOuliars', 3,...
    'shape',true,...
    'smoothing', true,...
    'limitElev', 5,...
    'minRange', 5000,...
    'maxAlt', 3000,...
    'data_res', 2);

% save('temp2_3.mat','rd','scans_filepath','scan_radar','scan_time','-v7.3')
% load("temp2_3")

%% Load df
r_lim = 150000; 
gridsize = num2str(r_lim/1000)+"km";
[df, rd] = load_F_elev("data/" + gridsize + "/", rd);

%% Prepare the grid
g = create_grid(df, rd);

[gX,gY,~]=meshgrid(g.x,g.y,g.z);
[gLON,gLAT,gZ]=meshgrid(g.lon,g.lat,g.z);

%% Compute the F matrix and reshape the data in sparse
[F, data, info] = construct_F(rd, df, g);
szf=size(F);
zmax = quantile(nonzeros(data),.999);

%% Compute geometry stuff and prepare everything for the inversion

% Load water mask
load("data/onwater/Fonwater.mat");
onwater = logical(Fonwater({g.f_lat(g.y),g.f_lon(g.x)}))'; clear Fonwater
% figure; hold on; imagesc(g.f_lon(g.x),g.f_lat(g.y),onwater'); tmp=axis;borders('states','w'); plot(radar_lon, radar_lat,".r"); axis equal tight; axis(tmp);


% Compute the VPR inversion matrix W
U = repelem(speye(g.sz(3)),prod(g.sz(1:2)),1);
W = F*U;

% Compute the lapacian smoothing for the VPR
Lh = LaplacianMatrix1D(g.sz(3));
Lh = Lh - diag(sum(Lh,2));

L3D = LaplacianMatrix(g, ...
    "water", false,...
    "coeff_horver", 0);


%%

DgridzI=nan(g.sz(3),numel(scan_time));
Ainv = cell(numel(scan_time),1);
b=cell(numel(scan_time),1);
Aspeed = repmat({sparse(prod(g.sz),prod(g.sz))},numel(scan_time),numel(scan_time));

uv_u = nan([g.sz numel(scan_time)]);
uv_v = nan([g.sz numel(scan_time)]);

for i_t=1:numel(scan_time)
    % retrieve id of non-zeros
    id_nnz = find(data(:,i_t));

    % keep only non-zeros0
    D = data(id_nnz,i_t);

    % VPR
    alpha=10;
    Avpr=[W(id_nnz,:);alpha*Lh];
    bvpr=[D;sparse(g.sz(3)*height(rd),1)];
    DgridzI(:,i_t) = lsqr(Avpr,bvpr,[],100);

    % Ratio
    % AgridzI = DgridzI(:,i_t) ./ sum(DgridzI(:,i_t));
    % Q = repmat(speye(prod(g.sz(1:2))),g.sz(3),g.sz(3)).*(U*AgridzI)'-speye(prod(g.sz));

    % Inverse system
    alpha1=1;
    alpha2=10;
    Ainv{i_t}=[F(id_nnz,:); ...
        ...% alpha1*Q;...
        ...% speye(prod(g.sz));...
        alpha2*L3D];
    b{i_t}=[D; ...
        ...% U*DgridzI(:,i_t) ;
        sparse(szf(2),1)];

    %     Dinvtmp = lsqr(Ainv{i_t}, b{i_t},[],100);
    %     Dinv(:,:,:,i_t)=reshape(Dinvtmp,g.sz);


    % Speed interpolation
    uv_u(:,:,:,i_t) = repmat(permute(rd{i_t}.uv(:,1),[2,3,1]),g.sz(1),g.sz(2),1);
    uv_v(:,:,:,i_t) = repmat(permute(rd{i_t}.uv(:,2),[2,3,1]),g.sz(1),g.sz(2),1);


    if i_t>1

        tmp_u = nanmean(uv_u(:,:,:,i_t-1:i_t),4) ./ g.dx .* seconds((diff(g.time(i_t-1:i_t))));
        tmp_v = nanmean(uv_v(:,:,:,i_t-1:i_t),4) ./ g.dy .* seconds((diff(g.time(i_t-1:i_t))));


        tmp_u = round(reshape(tmp_u,[],g.sz(3)));
        tmp_v = round(reshape(tmp_v,[],g.sz(3)));

        [row, col] = ind2sub(g.sz(1:2),(1:prod(g.sz(1:2)))');

        row = row+tmp_u;
        col = col+tmp_v;
        id_edge = row<1|row>g.sz(1)|col<1|col>g.sz(2);
        row(row<1)=1;
        row(row>g.sz(1))=g.sz(1);
        col(col<1)=1;
        col(col>g.sz(2))=g.sz(2);

        % Penilize more change in number over water
        dest_id = sub2ind(g.sz(1:2),row,col);
        coeff = ones(size(dest_id));
        coeff(onwater(:) & ~onwater(dest_id)) = 0.3; % Going from water to land
        % coeff(~onwater(:) & onwater(dest_id)) = 0.1; % Going from land to water
        coeff(~onwater(:) & ~onwater(dest_id)) = 0.3; % Going from land to land

        % Allow change at the boundard
        coeff(id_edge) = 0.1;

        id = sub2ind(g.sz,row(:),col(:),repelem((1:g.sz(3))',prod(g.sz(1:2))));
        id_full=1:prod(g.sz);
        Aspeed{i_t-1,i_t} = sparse(id_full(~isnan(id)),id(~isnan(id)),-coeff(~isnan(id)),prod(g.sz),prod(g.sz));
        Aspeed{i_t-1,i_t-1} = coeff(:).*speye(prod(g.sz));
        Aspeed{i_t-1,i_t-1}(isnan(id),:)=0;
    end
end

%% Inversion
alphat=5;
Dinit = U*DgridzI;
Dinit(:) = 0;
%Dinit = Dinvtmp;
Dinvtmp = lsqr([blkdiag(Ainv{:});alphat*cell2mat(Aspeed)],[cell2mat(b);sparse(prod(g.sz)*numel(scan_time),1)],[],200,[],[],Dinit(:));
Dinv=reshape(Dinvtmp,[g.sz numel(scan_time)]);

%% Flow

rho = squeeze(sum(Dinv,3)/(1000^3)); % bird/km^3 -> bird/m^3
vx = squeeze(sum(Dinv.*uv_u,3))./rho;
vy = squeeze(sum(Dinv.*uv_v,3))./rho;

[takingoff, landing, entering, leaving, gext, MVT] = sinksource(g,rho,vx,vy);

%% Save
% save("data/invfullnight"+scan_radar+scan_time, "Dinv", "data", "info", "df")


%%
url = "https://api.mapbox.com/styles/v1/mapbox/satellite-v9/static/["+g.lon(1)+","+g.lat(1)+","+g.lon(end)+","+g.lat(end)+"]/1280x"+round(1280*range(g.lat)/range(g.lon))+"@2x?access_token=pk.eyJ1IjoicmFmbnVzcyIsImEiOiIzMVE1dnc0In0.3FNMKIlQ_afYktqki-6m0g";
mapbox_background = imread(url);
% figure; imagesc(g.lon,g.lat,url_img); axis tight equal

%%

% find the optimal max density for plot
zmax = quantile(nonzeros(data),.99);

figure('position',[0 0 1200 1000],'Color',[48 48 48]/256);
t=tiledlayout(3,2,'TileSpacing','tight','Padding','none');

% nexttile;
% speed = cellfun(@(x) x.uv_dxy, rd(2:end), 'UniformOutput', false);
% imagesc(scan_time, squeeze(g.z),speed)

n1 = nexttile([2 1]);
i_elev=1;
dr=1;
daz=1;
sw = rd{1}.dens.sweeps(i_elev);
id_max =find(sw.range<r_lim,1,'last');
[az,r]=meshgrid(sw.az(1:daz:end),sw.range(1:dr:id_max));
dens = sw.data(1:dr:id_max, 1:daz:end );
ps = polarscatter(deg2rad(az(:)),r(:),r(:)/1000,dens(:),'filled');
hold on;
uvp = sum(rd{1}.uv.* [1 1j],2)*diff(df{1}.x(1:2));
q1 = polarplot([zeros(size(uvp)) pi/2-angle(uvp)]', [zeros(size(uvp)) abs(uvp)]',"color",[.3 .3 .3]);
q2 = polarplot([0 mean(pi/2-angle(uvp))], [0 mean(abs(uvp))],"k","linewidth",3);
set(gca,'ThetaZeroLocation',"top",'ThetaDir','clockwise'); 
set(gca,"ThetaTickLabel","","rticklabels","")
% set(gca,"Color",'none');
set(gca,"RLim",[0 r_lim])
clim([0 zmax]);
n1 =gca;

n2 = nexttile([2 1]); box on;
set(gca,"Color",'none');
imagesc(g.lon,g.lat,mapbox_background);
axis equal tight; a_xis=axis;
set(gca,"XColor","w","YColor","w")
%n2.XTick = [];
%n2.YTick = [];
% borders('states','w');borders("cuba",'k'); axis(a_xis)
% set(gca,'ydir','normal')

ax2 = axes(t);
set(gca,"Color",'none');
ax2.Layout=n2.Layout;
thr_im = 150;
im = imagesc(g.lon,g.lat,sum(Dinv(:,:,:,1),3),'Alphadata', (sum(Dinv(:,:,:,1),3))>thr_im);
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
axis equal tight; a_xis=axis;
set(gca,'ydir','normal')
clim([0 8*zmax]); colormap("jet")
set(gca,"XColor","w","YColor","w")


n3=nexttile;
set(gca,"Color",'none');
im2 = imagesc(g.lon,squeeze(g.z),squeeze(Dinv(:,end/2,:,1))');
set(gca,"ydir","normal")
clim([0 zmax]);
set(gca,"XColor","w","YColor","w")

nexttile;
set(gca,"Color",'none');
im3 = imagesc(g.lat,squeeze(g.z),squeeze(Dinv(end/2,:,:,1))');
set(gca,"ydir","normal")
clim([0 zmax]);
set(gca,"XColor","w","YColor","w")


for i_t=(4*3):numel(scan_time)
    title(n1,datestr(scan_time(i_t)),FontSize=20,Color="w")
    tmp = sum(Dinv(:,:,:,i_t),3)';
    im.CData = tmp;
    im.AlphaData = (tmp-thr_im)/200.*(tmp>thr_im);

    dens = rd{i_t}.dens.sweeps(i_elev).data(1:dr:id_max, 1:daz:end );
    ps.CData = dens(:);
    uvp = sum(rd{i_t}.uv.* [1 1j],2)*diff(df{1}.x(1:2));
    for i_q=1:numel(q1)
        q1(i_q).ThetaData = [0 pi/2-angle(uvp(i_q))];
        q1(i_q).RData = [0 abs(uvp(i_q))];
    end
    q2.ThetaData = [0 mean(pi/2-angle(uvp))];
    q2.RData = [0 mean(abs(uvp))];

    im2.CData = squeeze(Dinv(:,end/2,:,i_t))';

    im3.CData = squeeze(Dinv(end/2,:,:,i_t))';

   
   pause(.1)
   % exportgraphics(gcf,scan_radar+"_"+datestr(mean(scan_time),'YYYYMMDD')+".gif",'Append',true,'BackgroundColor',[48 48 48]/256);
end



%% 
% i_t=20;
% figure('position',[0 0 2400 800]);
% tiledlayout(1,3,'TileSpacing','tight','Padding','none')
% for i=1:3
%     nexttile;
%     im = imagesc(g.lon,g.lat,sum(Dinv(:,:,:,i_t+(i-2)*2),3)');
%     axis equal tight; a_xis=axis;
%     borders('states','w','linewidth',2); axis(a_xis)
%     set(gca,'ydir','normal')
%     clim([0 8*zmax]); xticklabels('');yticklabels('')
%     title(datestr(scan_time(i_t+(i-2)*2)),'FontSize',24,'Color','w'); box off
%     set(gca,'Visible','off')
% end

figure('position',[0 0 800 800],'Color',[48 48 48]/256);
ax1 = axes; box on;
imagesc(g.lon,g.lat,mapbox_background);
axis equal tight; a_xis=axis;
set(gca,"XColor","w","YColor","w")
% borders('states','w');borders("cuba",'k'); axis(a_xis)
% set(gca,'ydir','normal')

ax2 = axes;linkaxes([ax1 ax2])
hold on; clim([100 8*zmax]); colormap("jet")
i_t=[10 20 27];
thr_im = 150;
 for i=1:numel(i_t)
     tmp = sum(Dinv(:,:,:,i_t(i)),3)';
     tmp2=quantile(tmp(:),.95);
     im = imagesc(g.lon,g.lat,tmp,'Alphadata', (tmp-tmp2)/200.*(tmp>tmp2));
     % contour(g.lon,g.lat,tmp,[-1000 quantile(tmp(:),.99)],'r');
 end
 ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
axis equal tight; axis(a_xis);
set(gca,'ydir','normal')
set(gca,"XColor","w","YColor","w")


%%
c_lim = 1e10;
clmap = brewermap([],'Spectral');

figure('position',[0 0 1600 800]); colormap(clmap)

nexttile; hold on;
im1 = imagesc(g.lon,g.lat,takingoff(:,:,1)');
axis equal tight;  a_xis=axis;
borders('states','k'); axis(a_xis); 
% set(gca,'ColorScale','log')
clim([-c_lim c_lim]);
title("Taking-off "+datestr(scan_time(i_t)),'FontSize',24,'Color','w'); box off
xticklabels(''); yticklabels('');

nexttile; hold on;
im2 = imagesc(g.lon,g.lat,landing(:,:,1)');
axis equal tight;  a_xis=axis;
borders('states','k'); axis(a_xis); 
% set(gca,'ColorScale','log')
clim([-c_lim c_lim]);
title("Landing "+datestr(scan_time(i_t)),'FontSize',24,'Color','w'); box off
xticklabels(''); yticklabels('');

for i_t=1:numel(scan_time)-1
    title(datestr(scan_time(i_t)))
    im1.CData = sum(takingoff(:,:,1:i_t),3)';
    im2.CData = sum(landing(:,:,1:i_t),3)';
    pause(.2)
   % exportgraphics(gcf,"test_landing.gif",'Append',true);
end