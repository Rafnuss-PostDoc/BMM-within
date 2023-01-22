function g = create_grid(df, rd, varargin)
% create_grid() return a structure with grid information computed from the
% radar position and Felev

radar_lon = cellfun(@(x) x.lon, rd(:,1));
radar_lat = cellfun(@(x) x.lat, rd(:,1));

tmp_dlat = deg2km(distance(mean(radar_lat),mean(radar_lon),mean(radar_lat)+1,mean(radar_lon)))*1000;
tmp_dlon = deg2km(distance(mean(radar_lat),mean(radar_lon),mean(radar_lat),mean(radar_lon)+1))*1000;

radar_dlat = max(df{1}.y)/tmp_dlat; 
radar_dlon = max(df{1}.x)/tmp_dlon; 

lonmin=min(radar_lon)-radar_dlon;
lonmax=max(radar_lon)+radar_dlon;
latmin=min(radar_lat)-radar_dlat;
latmax=max(radar_lat)+radar_dlat;

% overwrite the grid if provided. 
if nargin>2
    lonmin=varargin{1};
end
if nargin>3
    lonmax=varargin{2};
end
if nargin>4
    latmin=varargin{3};
end
if nargin>5
    latmax=varargin{4};
end

% Find how much 1km=degree
lat_mean=(latmax+latmin)/2;
lon_mean=(lonmax+lonmin)/2;
g.dlat = deg2km(distance(lat_mean,lon_mean,lat_mean+1,lon_mean))*1000;
g.dlon = deg2km(distance(lat_mean,lon_mean,lat_mean,lon_mean+1))*1000;

% 
xmax = (lonmax-lonmin)*g.dlon;
xmax = ceil(xmax/10000)*10000;
ymax = (latmax-latmin)*g.dlat;
ymax = ceil(ymax/10000)*10000;

g.dx = diff(df{1}.x(1:2));
g.x=(0:ceil(xmax/g.dx)+3)*g.dx;
g.dy = unique(diff(df{1}.y));
g.y=(0:ceil(ymax/g.dy)+3)'*g.dy;
g.z=df{1}.z;

g.sz = [numel(g.x) numel(g.y) numel(g.z)];

% get lat lon coordinates
g.f_lon = @(X) X./g.dlon+lonmin;
g.f_lat = @(Y) Y./g.dlat+latmin;
g.f_X = @(lon) (lon-lonmin)*g.dlon;
g.f_Y = @(lat) (lat-latmin)*g.dlat;

g.lat=g.f_lat(g.y);
g.lon=g.f_lon(g.x);

% for i_r=1:numel(rd)
%     [~, rd{i_r}.radar_lon_id] = min(abs(g.lon(:)' - rd{i_r}.lon(:)),[],2);
%     [~, rd{i_r}.radar_lat_id] = min(abs(g.lat(:)' - rd{i_r}.lat(:)),[],2);
% end

g.time = cellfun(@(x) x.datetime, rd(1,:));

% Load water mask
load("data/onwater/Fonwater.mat");
g.mask_water = logical(Fonwater({g.f_lat(g.y),g.f_lon(g.x)}))'; clear Fonwater
% figure; hold on; imagesc(g.f_lon(g.x),g.f_lat(g.y),onwater'); tmp=axis;borders('states','w'); plot(radar_lon, radar_lat,".r"); axis equal tight; axis(tmp);


end