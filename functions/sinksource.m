function [takingoff, landing, entering, leaving, gext, MVT] = sinksource(g,rho,vx,vy)

%%

g.ny = numel(g.y);
g.nx = numel(g.x);
g.area = repmat(g.dx*g.dy,g.nx,g.ny);

assert(all(size(rho)==[g.nx g.ny numel(g.time)]))
assert(all(size(vx)==[g.nx g.ny numel(g.time)]))
assert(all(size(vy)==[g.nx g.ny numel(g.time)]))

% *Build extended grid*
%
% Build a grid |gext| with a cell more on each side of the domain (ext=extanded)
% in order to compute the fluxes in/out of the existing grid

gext.y = fillmissing([nan(1,1) ; g.y ;nan(1,1)],'linear');
gext.x = fillmissing([nan(1,1) ; g.x'; nan(1,1)],'linear');
gext.ny = numel(gext.y);
gext.nx = numel(gext.x);
% [gext.LAT, gext.LON] = ndgrid(gext.y,gext.x);
gext.mask_water = padarray(g.mask_water,[1 1],true);

%%%%%%%%
gext.mask_water = false(size(g.mask_water));

gext.mask_out=movmean(movmean(~gext.mask_water,[1 1],2),[1 1],1)>0;
gext.mask_out(~gext.mask_water)=false;

% use regular grid

gext.time = g.time;
round_min = 15;
gext.time.Minute = round((gext.time.Minute+(gext.time.Second)/60)/round_min)*round_min;
gext.time.Second = 0;
dt = unique(seconds(diff(gext.time)));
assert(numel(dt)==1)
gext.nt = numel(gext.time);


%% *Compute the flux* $\Phi =\mathit{\mathbf{v}}\rho$ *at +/- 1/2 grid cell.*
% First, add a nan layer in y or x direction
Phiy_pad = padarray(rho .* vy .* g.dx,[0 1 0],nan); % bird/km^2 * km/h *km ->  bird/h
Phix_pad = padarray(rho .* vx .* g.dy,[1 0 0],nan);

% and then, compute the flux at +/- 1/2 even if either previous or next cells
% is nan.
Phiy_h = movmean(Phiy_pad,[0 1],2,'omitnan','Endpoints','discard');
Phix_h = movmean(Phix_pad,[0 1],1,'omitnan','Endpoints','discard');
% clear Phix_pad Phiy_pad

% *Compute the delta flux / delta distance* $\frac{\Delta \Phi }{\Delta \left(\textrm{y},\textrm{x}\right)}$
% First, add 0 padding for the outer zone (allows to compute boundary cell)
Phiy_h_0=padarray(Phiy_h,[1 1 0],0);
Phix_h_0=padarray(Phix_h,[1 1 0],0);
% clear Phiy_h Phix_h

% Then, replace the nan by zero to compute the change of fluxes even no fluxes
% on the other side (i.e. boundary cell)
Phiy_h_0(isnan(Phiy_h_0))=0;
Phix_h_0(isnan(Phix_h_0))=0;

% Finally, compute the delta flux over delta distance.
dPhiy = diff(Phiy_h_0,1,2);
dPhix = diff(Phix_h_0,1,1);
% clear Phix_h_0 Phiy_h_0

% *Compute the variation of bird in/out each cell*
F = (dPhiy + dPhix ).*dt; % bird/h * hr -> bird
% clear dPhiy dPhix

% Crank-nicolson Move F at t+1/2:
% F(F==0)=nan;
% F = movmean( F ,[0 1],3,'omitnan');
% F(isnan(F))=0;

% F comprise both the intral change of bird and the bird going out of the area
% Inner flux: remove the boundary cells
Fin = F;
%Fin(~gext.mask_water)=0;
Fin(repmat(gext.mask_water,1,1,gext.nt)) = 0;

% Outer Flux. Note that Phiy_h_0 and Phiy_h_0 are 0 for the cell outside
% the domain, such at dPhiydy is equal to the flux at the outer cell (with
% the sign corresponding to the direction).
Fout = F;
%Fout.Fout(gext.mask_water)=0;
Fout(repmat(~gext.mask_water,1,1,gext.nt)) = 0;
% clear F

% *Compute the landing/derpature* $W$
% |rho0| is used (NaN replaced by 0) to account for landing as sunset and sunrise.
% Unit of W is [bird =bird (*15min/15min)], when (+) Birds entering the airspace
% (takingoff) and (-) Birds leaving the raispace (landing)
rho0 = rho; rho0(isnan(rho)) = 0; % bird/km^2
W = diff(rho0 .* repmat(g.area,1,1,gext.nt),1,3) + Fin(2:end-1,2:end-1,1:end-1); % bird/km^2 * km^2  + bird -> bird
% clear rho0 Fin

% *Check methodology*
%
% Check for a step
if false
    i=20;
    gext.time(i);
    figure('position',[0 0 1600 600]); tiledlayout('flow','TileSpacing','tight','Padding','tight');
    nexttile; imagesc(g.lon,g.lat,rho(:,:,i)','AlphaData',~isnan(rho(:,:,i)')); set(gca,'ydir','normal'); colorbar; title('\rho'); axis equal tight;  a_xis=axis; borders('states','w'); axis(a_xis);
    nexttile; imagesc(g.lon,g.lat,vy(:,:,i)','AlphaData',~isnan(vy(:,:,i)')); set(gca,'ydir','normal'); colorbar; title('vy'); axis equal tight; borders('states','w'); axis(a_xis);
    nexttile; imagesc(g.lon,g.lat,vx(:,:,i)','AlphaData',~isnan(vx(:,:,i)')); set(gca,'ydir','normal'); colorbar; title('vx'); axis equal tight; borders('states','w'); axis(a_xis);
    
    figure('position',[0 0 1200 600]); tiledlayout('flow','TileSpacing','tight','Padding','tight');
    nexttile; imagesc(g.lon,g.lat,Phix_pad(:,:,i)','AlphaData',~isnan(Phix_pad(:,:,i)')); set(gca,'ydir','normal'); colorbar; title('\Phi_{x}');axis equal tight; borders('states','w'); axis(a_xis);
    nexttile; imagesc(g.lon,g.lat,Phiy_pad(:,:,i)','AlphaData',~isnan(Phiy_pad(:,:,i)')); set(gca,'ydir','normal'); colorbar; title('\Phi_{y}'); axis equal tight; borders('states','w'); axis(a_xis);
    nexttile; imagesc(g.lon,g.lat,dPhix(:,:,i)','AlphaData',~(0==dPhix(:,:,i)')); set(gca,'ydir','normal'); colorbar; title('\Delta\Phi_{x}/\Delta x'); axis equal tight; borders('states','w'); axis(a_xis);
    nexttile; imagesc(g.lon,g.lat,dPhiy(:,:,i)','AlphaData',~(0==dPhiy(:,:,i)')); set(gca,'ydir','normal'); colorbar; title('\Delta\Phi_{y}/\Delta y'); axis equal tight; borders('states','w'); axis(a_xis);
    nexttile; imagesc(g.lon,g.lat,F(:,:,i)','AlphaData',~(0==F(:,:,i)')); set(gca,'ydir','normal'); colorbar; title('F=(\Delta\Phi_{y}+\Delta\Phi_{x})\Delta t'); axis equal tight; borders('states','w'); axis(a_xis);
    
    % Check the spatial aggregated mass balence.
    %     a = reshape(sum(sum(Fout,1),2),1,[]);
    %     b = reshape(sum(sum(rho0 .* repmat(g.area,1,1,gext.nt),1),2),1,[]);
    %     c = reshape(sum(sum(W,1),2),1,[]);
    %     figure('position',[0 0 1000 200]); hold on; xlabel('time');ylabel('Error')
    %     plot(b(2:end) - b(1:end-1) - c - a(1:end-1))
    
    % See a single timestep.
    figure('position',[0 0 1300 600]);tiledlayout('flow','TileSpacing','tight','Padding','tight');
    nexttile; imagesc(g.lon,g.lat,rho0(:,:,i)','AlphaData',~(0==rho0(:,:,i)')); title('\rho(t)');set(gca,'ydir','normal'); colorbar;axis equal tight; borders('states','w'); axis(a_xis);
    nexttile; imagesc(g.lon,g.lat,rho0(:,:,i+1)','AlphaData',~(0==rho0(:,:,i+1)')); title('\rho(t+1)'); set(gca,'ydir','normal'); colorbar;axis equal tight; borders('states','w'); axis(a_xis);
    nexttile; imagesc(g.lon,g.lat,rho0(:,:,i+1)'-rho0(:,:,i)','AlphaData',~(0==rho0(:,:,i)')); title('\rho(t+1) - \rho(t)'); set(gca,'ydir','normal'); colorbar;axis equal tight; borders('states','w'); axis(a_xis);
    nexttile; imagesc(g.lon,g.lat,Fin(:,:,i)','AlphaData',~(0==Fin(:,:,i)'));  title('F_{in}'); set(gca,'ydir','normal'); colorbar;axis equal tight; borders('states','w'); axis(a_xis);
    nexttile; imagesc(g.lon,g.lat,W(:,:,i)','AlphaData',~(0==W(:,:,i)'));  title('W^{t\rightarrow t+1}=\rho(t+1) - \rho(t)+F_{in}'); set(gca,'ydir','normal'); colorbar;axis equal tight; borders('states','w'); axis(a_xis);
    nexttile; imagesc(g.lon,g.lat,Fout(:,:,i)','AlphaData',~(0==Fout(:,:,i)'));  title('F_{out}'); set(gca,'ydir','normal'); colorbar;axis equal tight; borders('states','w'); axis(a_xis);

end
%% *Migratory Processes*
% *Seperature the landing/landing (-) and takingoff/take-off (+)*
landing=W;
landing(landing>=0)=0;
% landing2
%     landing2 = landing;
%     landing2(gext.mask_water(2:end-1,2:end-1,1:end-1))=0;
takingoff=W;
takingoff(takingoff<=0)=0;
leaving=Fout;
leaving(leaving>=0)=0;
entering=Fout;
entering(entering<=0)=0;

% MTR
MVT = rho .* sqrt(vy.^2 + vx.^2) .* dt; % bird/km^2 .* km/hr .* hr -> bird/km

end