addpath('/Users/raphael/Documents/GitHub/wsrlib');
wsrlib_setup();
addpath(genpath('functions/'))

scans_filepath = load_rd("KBYX", datetime("19-Apr-2019 01:00"));
rd = read_rd(scans_filepath, 'maxAlt', 2000);

%%

r_lim = 50000; 
gridsize = num2str(r_lim/1000)+"km";
[df, rd] = load_F_elev("data/" + gridsize + "/", rd);
g = create_grid(df, rd);
[gX,gY,~]=meshgrid(g.x,g.y,g.z);
[gLON,gLAT,gZ]=meshgrid(g.lon,g.lat,g.z);
[F, ~, info] = construct_F(rd, df, g);
szf=size(F);


%% Load scan pattern

% Set grid size
nx=120;
ny=40;
%nz=10;
nt=12*60/5;

% Generate Departure
addpath('/Users/raphael/Documents/Github/FastGaussianSimulation')
sim.s = g.sz;
sim.n = 1;
covar(1).model = 'k-bessel';
covar(1).range = [40 40 1]./[diff(g.x(1:2)) diff(g.y(1:2)) diff(g.z(1:2))].*1000; % meter
covar(1).azimuth = [0 0];
covar(1).var = 1;
covar(1).alpha = 1.2;
% covar(2).model = 'k-bessel';
% covar(2).range = [1 1 1];
% covar(2).azimuth = [0 0];
% covar(2).var = 0.01;
% covar(2).alpha = 5;


% vertical profile
dfz = squeeze(g.z);
vpr = exp(-1*(2+dfz/1000));
vpr = pdf(makedist('Lognormal','mu',log(dfz(12)),'sigma',1), dfz);
vpr = vpr./sum(vpr,"all",'omitnan');
% vp = 100*(dfz<1000);
figure; plot(vpr,dfz,'-or'); xlabel('VPR'); ylabel('Altitude'); title('vertical profile')


sim.seed=28;
Dfor0 = FGS(sim,covar);
Dfor = Dfor0{1}.*5+20;
Dfor = z_to_refl(idb(Dfor))/11;
Dfor = Dfor.*reshape(vpr,size(g.z));
max(Dfor(:))

% figure; hold on;
% plot(Dfor(:),Z(:),'.k')
% plot(squeeze(mean(Dfor,[1 2])),dfz,'-or'); 
% xlabel('VPR'); ylabel('Altitude'); title('vertical profile')

c_lim = quantile(Dfor(:),.9);
%% 

figure;
Dfor2=Dfor; Dfor2(Dfor<c_lim)=nan;
h = slice(gX,gY,gZ,Dfor2,g.x(end/2),g.y(end/2),[150]);
clim([c_lim max(Dfor(:))])
set(h,'facealpha',.95)
set(h,'edgecolor','none')
%caxis([0 z_to_refl(idb(30))/11]);
view(3); %caxis([0 5])
box on; grid on; zlim([0 2000])
xticks(g.x(1:100:end)); yticks(g.y(1:100:end));
zticks(0:1000:2000)
Ax=gca; Ax.XColor = 'w'; Ax.YColor = 'W'; Ax.ZColor = 'W'; Ax.Color = 'none';%[50 50 50]/255;
axis tight;
colormap("jet")

figure; hold on;
Dfor2=Dfor(end/2:end,:,:);
col=colormap("jet");
c_lim2 = linspace(c_lim,max(Dfor2(:)),5);
for i=1:numel(c_lim2)
    p=patch(isosurface(gX(end/2:end,:,:),gY(end/2:end,:,:),gZ(end/2:end,:,:),Dfor2,c_lim2(i)),...
        'EdgeColor','none','FaceColor',col(round(i*255/numel(c_lim2)),:), "facealpha",0.4);
end
clim([c_lim max(Dfor2(:))])
p2=patch(isocaps(gX(end/2:end,:,:),gY(end/2:end,:,:),gZ(end/2:end,:,:),Dfor2,c_lim),...
    'FaceColor','interp','EdgeColor','none', 'FaceAlpha',0.8);
view(3)
box on; grid on; zlim([0 3000])
xticks(g.x(1:100:end)); yticks(g.y(1:100:end));
zticks(0:1000:3000)
Ax=gca; Ax.XColor = 'w'; Ax.YColor = 'W'; Ax.ZColor = 'W'; Ax.Color = 'none';%[50 50 50]/255;
Ax.DataAspectRatio = [1 1 1/15];
axis tight;
% camlight; lighting gouraud

% cdata = smooth3(rand(size(Dfor)),'box',7);
% isonormals(gX,gY,gZ,Dfor,p)
% isocolors(gX,gY,gZ,Dfor,p)
% p.FaceColor = 'interp';
% p.EdgeColor = 'none';

exportgraphics(gcf,'test-1.png','BackgroundColor',[50 50 50]/255)
%% Forward
obs=nan(info.sz_rd);
obs(:)=F*Dfor(:);

figure; tiledlayout('flow','TileSpacing','tight','Padding','tight'); hold on;
dr=4; daz=2; max_z=5000; min_r=5000;
for i_elev=1:1:min(6,numel(df))

    [~, z] = slant2ground(df{i_elev}.range, df{i_elev}.elev);
    az = df{i_elev}.az;

     id_max_z = find(max_z<z,1);
     id_min_r = find(min_r<df{i_elev}.range,1);
    id_rz = id_min_r:dr:id_max_z;
    id_az = 1:daz:numel(az);

    [az, range, z] = expand_coords(az(id_az), df{i_elev}.range(id_rz), z(id_rz));

    [x,y,z] = pol2cart(az/180*pi,range,z);

    data = obs(id_rz,id_az,1,i_elev);

    id = data>0 & y>0;
    scatter3(x(id),y(id),z(id),10,data(id),'filled')% r(id)/5000
    
end
box on; grid on; zlim([0 2000])
xlim([-50 50]*1000);ylim([0 50]*1000)
xticks([-50 0 50]*1000); yticks([-50 0 50]*1000)
zticks(0:1000:2000)
view(3)
Ax=gca; Ax.XColor = 'w'; Ax.YColor = 'W'; Ax.ZColor = 'W'; Ax.Color = [50 50 50]/255;
Ax.DataAspectRatio = [1 1 1/15];
% exportgraphics(gcf,'test-2.png','BackgroundColor',[50 50 50]/255)


