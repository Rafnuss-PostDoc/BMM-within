
% Run from server with
% module load matlab/R2021a
% matlab -nodisplay -nosplash -nodesktop -nojvm -r "elev_iter = $1 ;run('generate_F_elev.m');exit;"

%% coarse grid
dx=2000; 
dy=dx; 
xy_max=150000;
dz=100;
z_0=dz/2;
z_max=3000-dz/2;
x=reshape(-xy_max:dx:xy_max,1,[]);
y=reshape(-xy_max:dy:xy_max,[],1);
z=reshape(z_0:dz:z_max,1,1,[]);
nxx=40;
nyy=40;
nzz=5;

%%
dx=500; 
dy=dx; 
xy_max=50000;
dz=100;
z_0=dz/2;
z_max=3000-dz/2;
x=reshape(-xy_max:dx:xy_max,1,[]);
y=reshape(-xy_max:dy:xy_max,[],1);
z=reshape(z_0:dz:z_max,1,1,[]);



%%
load("data/scan_pattern.mat")
%load("scan_pattern/new.mat")

nxx=10;
nyy=10;
nzz=5;

compute_F_elev(x,y,z, nxx, nyy, nzz, dfnew, 0.1)