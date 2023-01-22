
load('../BMM-US/data/density/inference-trans.mat','radar');
load('../BMM-US/data/DEM/dem');
Fdem = griddedInterpolant({dem_lat,dem_lon},double(dem),'linear');
a_xis = [-125 -67 24.5 50];

i_r=radar.name =="KLWX";

g.lon=radar.lon(i_r)+[-50 50]/111;
g.lat=radar.lat(i_r)+[-50 50]/111;
url = "https://api.mapbox.com/styles/v1/mapbox/satellite-v9/static/["+g.lon(1)+","+g.lat(1)+","+g.lon(end)+","+g.lat(end)+"]/800x800?access_token=pk.eyJ1IjoicmFmbnVzcyIsImEiOiIzMVE1dnc0In0.3FNMKIlQ_afYktqki-6m0g&logo=false";
[mapbox_background,map] = imread(url);
% figure; imagesc(g.lon,g.lat,mapbox_background); axis tight equal


dem_radar = Fdem({linspace(g.lat(1), g.lat(2),2560), linspace(g.lon(1), g.lon(2),2560)});
dem_radar = imgaussfilt(dem_radar,3,'FilterSize',11);

surf(linspace(-50,50,2560), linspace(-50,50,2560),dem_radar,'EdgeColor','none')
demcmap(dem_radar); zlim([0 2000])
[X,Y]=meshgrid(linspace(-50000,50000,2560), linspace(-50000,50000,2560));


%%

figure(Position=[0 0 900 450]);  hold on;
warp(X+50000,Y+50000,dem_radar'-100,permute(mapbox_background,[2,1,3]));

Dfor2=Dfor(:,end/2:end,:);
col=colormap("jet");
c_lim2 = linspace(c_lim,max(Dfor2(:)),5);
for i=1:numel(c_lim2)
    p=patch(isosurface(gX(:,end/2:end,:),gY(:,end/2:end,:),gZ(:,end/2:end,:),Dfor2,c_lim2(i)),...
        'EdgeColor','none','FaceColor',col(round(i*255/numel(c_lim2)),:), "facealpha",0.4);
end
clim([c_lim max(Dfor2(:))])
p2=patch(isocaps(gX(:,end/2:end,:),gY(:,end/2:end,:),gZ(:,end/2:end,:),Dfor2,c_lim),...
    'FaceColor','interp','EdgeColor','none', 'FaceAlpha',0.8);
scatter3(50000,50000, 20,300,"w",'filled','MarkerEdgeColor','k')
view(3)
box on; grid on;
xticks(g.x(1:100:end)); yticks(g.y(1:100:end));zticks(0:1000:3000)
xticklabels([-50 0 50]); yticklabels([-50 0 50]);
Ax=gca; % Ax.XColor = 'w'; Ax.YColor = 'W'; Ax.ZColor = 'W';  
Ax.Color = [230 230 230]/255;
Ax.DataAspectRatio = [1 1 1/15];
xlim([0 100000]); ylim([0 100000]); zlim([-100 2000])

view(-80,9)
exportgraphics(gcf,'test-1.png','BackgroundColor',[50 50 50]/255)

%% 
figure(Position=[0 0 900 450]); hold on;
warp(X+50000,Y+50000,dem_radar'-100,permute(mapbox_background,[2,1,3]));

obs=nan(info.sz_rd);
obs(:)=F*Dfor(:);

dr=4; daz=4; max_z=5000; min_r=5000;
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

    id = data>0 & x>0;
    scatter3(x(id)+50000,y(id)+50000, z(id),30,data(id),'filled','MarkerFaceAlpha',0.8)% r(id)/5000
    
end
scatter3(50000,50000, 20,300,"w",'filled','MarkerEdgeColor','k')
colormap("jet");
clim([c_lim max(Dfor2(:))])
box on; grid on;
xticks(g.x(1:100:end)); yticks(g.y(1:100:end));zticks(0:1000:3000)
xticklabels([-50 0 50]); yticklabels([-50 0 50]);
Ax=gca; % Ax.XColor = 'w'; Ax.YColor = 'W'; Ax.ZColor = 'W'; 
Ax.Color = [230 230 230]/255;
Ax.DataAspectRatio = [1 1 1/15];
xlim([0 100000]); ylim([0 100000]); zlim([-100 2000])

view(-60,11)
view(-80,9)
exportgraphics(gcf,'test-2.png','BackgroundColor',[50 50 50]/255)
