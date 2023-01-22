function compute_F_elev(x, y, z, nxx, nyy, nzz, dfe)
% compute_F_elev compute the matrix F for a specific radar elevation scan
% df and a grid specified by x,y,z. nxx, nyy, nzz are the downscaling of  
% the grid on which the power function is approximated. 

disp(dfe.name)

%% Define the grid
assert(numel(unique(diff(x)))==1)
dx = diff(x(1:2));
assert(numel(unique(diff(y)))==1)
dy = diff(y(1:2));
assert(numel(unique(diff(z)))==1)
dz = diff(z(1:2));

% Size of the coarse grid
sz_coarse_grid = [numel(x) numel(y) numel(z)];

% Fine grid
dxx=dx/nxx; dyy=dy/nyy; dzz=dz/nzz;
xx=(x(1)-dx/2+dxx/2):dxx:(x(end)+dx/2-dxx/2);
xx=single(reshape(xx,1,[]));
yy=(y(1)-dy/2+dyy/2):dyy:(y(end)+dy/2-dyy/2);
yy=single(reshape(yy,[],1));
zz=(z(1)-dz/2+dzz/2):dzz:(z(end)+dz/2-dzz/2);
zz=single(reshape(zz,1,1,[]));

% [XX,YY,ZZ] = meshgrid(xx,yy,zz);
sz_fine_grid = [numel(xx) numel(yy) numel(zz)];


%% Convert the cartisian grid to polar coordinate accounting for the radius of the earth
% xx,yy,zz -> SS,TT,PP

% range along the ground (great circle distance)
SS = sqrt(xx.^2 + yy.^2);

% set cst
earth_radius = 6371200; % from NARR GRIB file
multiplier = 4/3; % correct for atmospheric refraction
r_e = earth_radius * multiplier; % earth effective radius
% h height above ground is zz here

% document with equation to see the computation
sigma = SS/r_e; clear SS
h_p =(r_e+zz).*cos(sigma)-r_e;
A = (r_e+zz).*sin(sigma);
RR = sqrt(h_p.^2 + A.^2);
TT = rad2deg(asin(h_p ./ RR ));

PP = angle(xx + 1i*yy);%atan(yy./xx);% + pi.*(xx<0).*sign(yy)+pi/2;
PP = pol2cmp( PP );
PP = repmat(PP,1,1,numel(zz));

% coarse scale
% sigma = sqrt(x.^2 + y.^2 )/r_e;
% h_p =(r_e+z).*cos(sigma)-r_e;
% A = (r_e+z).*sin(sigma);
% R = sqrt(h_p.^2 + A.^2);
% T = rad2deg(asin(h_p ./ R ));
% P = angle(x + 1i*y);%atan(yy./xx);% + pi.*(xx<0).*sign(yy)+pi/2;
% P = pol2cmp( P );
% P = repmat(P,1,1,numel(z));

clear sigma h_p A


%% Compute weighting matrix

% compute the std dev of the beam along phi and theta
hpbw = 0.925; % half power beamwiths
simga_tp = hpbw/2/sqrt(2*log(2));


theta0 = dfe.elev;
id_theta = find(TT>(theta0-3*simga_tp) & TT<(theta0+3*simga_tp));


% we can compute only for the datapoint which fall within our grid
i_range_max = sum(dfe.range<=max(RR(id_theta)));
simga_r = 0.35*median(diff(dfe.range));

Felevf = cell(i_range_max,numel(dfe.az));
Felevid = cell(i_range_max,numel(dfe.az));

% Loop through all the ranges
for i_range = 1:i_range_max
    t2=tic;
    r0 = dfe.range(i_range);
    id_range = id_theta(abs(RR(id_theta)-r0)<=3*simga_r);

    for i_azim = 1:numel(dfe.az)

        phi0 = dfe.az(i_azim);

        diff_phi = mod(PP(id_range)-phi0,360);
        diff_phi = min(360-diff_phi, diff_phi);
        id_range_phi = diff_phi<3*simga_tp;
        id_phi = id_range(id_range_phi);

        if numel(id_phi)==0
            continue
        end

        % compute the weight and normalize them on the grid
        f = exp(-(RR(id_phi)-r0).^2/2/simga_r^2 -(TT(id_phi)-theta0).^2/simga_tp^2 - (diff_phi(id_range_phi)).^2/simga_tp^2 );

        % find the index on the fine grid
        [id_xx,id_yy, id_zz] = ind2sub(sz_fine_grid,id_phi);
        % figure; polarscatter(deg2rad(phi0),r0,100)
        % hold on;
        % polarscatter(deg2rad(PP(id_f)),RR(id_f),10,f,'filled');
        % rlim([46500 47200]); thetalim([15 25])
        % figure; scatter3(squeeze(xx(id_xx)),squeeze(yy(id_yy)),squeeze(zz(id_zz)),100,f,'filled')

        % convert the index to the large grid scale
        id_x=floor((id_xx-1)/(dx/dxx))+1;
        id_y=floor((id_yy-1)/(dy/dyy))+1;
        id_z=floor((id_zz-1)/(dz/dzz))+1;
        id_coarse = sub2ind(sz_coarse_grid,id_x,id_y,id_z);

        [G,id_g]=findgroups(id_coarse);

        % sum all f within the same coarse grid scale.
        f_g=splitapply(@sum,f,G);
        % count_g=splitapply(@sum,f>.01,G);

        % sort f by coarse grid scale
        [f_gs,id_s]=sort(f_g);

        % Keep grid scale which cumulative sum of f amount for more
        % than 10% of the sum
        id_ss=cumsum(f_gs)./sum(f_gs)>.1;

        % Normalize the sum of f to 1 for this observation.
        Felevf{i_range,i_azim} = f_gs(id_ss)./sum(f_gs(id_ss));
        Felevid{i_range,i_azim} = uint32([i_range*ones(sum(id_ss),1) i_azim*ones(sum(id_ss),1) id_g(id_s(id_ss)) ]);

        % Illustration
        if false
            figure; tiledlayout(2,1,'TileSpacing','tight','Padding','tight')
            nexttile; hold on;
            imagesc(xx,yy,tmp(:,:,1)); plot(0,0,'.k');
            yline(0);xline(0)
            set(gca,'ydir','normal'); % colormap(flipud(bone));
            axis tight equal; ylim([0 50000])

            nexttile; hold on; imagesc(xx,squeeze(zz),squeeze(tmp(:,end/2+.5,:))')
            axis tight equal; xlim([0 50000])
        end
    end

    disp(i_range+"/"+num2str(i_range_max)+" | "+num2str(round(toc(t1)))+"sec |"+num2str(round(toc(t2)))+"sec")
end

% Indices of F
tmp = reshape(Felevid,[],1);
id=cellfun(@(x) size(x,2),tmp)==3;
Felevidc=cell2mat(tmp(id));

% Value of F
tmp = reshape(Felevf,[],1);
Felevfc=cell2mat(tmp(id));

% Save
save("F_"+dfe.name+"_"+num2str(max(x)/1000)+"km.mat",'Felevidc','Felevfc','x','y','z','dfe')


%% Illustration

if false
    figure; hold on
    xline(0:2:6,'k','linewidth',2)
    yline(0:2:6,'k','linewidth',2)
    for i=2:2:6
        p = nsidedpoly(1000, 'Center', [0 0], 'Radius', i);
        plot(p, 'edgeColor', [.7 .7 .7],'facealpha',0)
    end
    for i=-180:45:180
        fplot( @(x) x*sind(i),'Color',[.7 .7 .7])
    end

    axis equal tight
    xlim([0 6]); ylim([0 4])

    gridsub=[2 0];
    for i=0:.2:2%0:.05:2
        plot(gridsub(1)+i*[1 1], gridsub(2)+[0 2],'Color',[.7 .7 .7])
        plot(gridsub(1)+[0 2], gridsub(2)+i*[1 1],'Color',[.7 .7 .7])
    end
end

end
