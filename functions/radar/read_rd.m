function rd = read_rd(scans_filepath, varargin)
% read_rd() read and proccess the scans from an array of string scans_filepath

p = inputParser;
addParameter(p,'removeRain','rho'); % method to remove rain 'rh'
addParameter(p,'removeOuliars',3); % MAD value which considered outlier
addParameter(p,'shape',true); % 
addParameter(p,'smoothing',true); % perform a smoothing on the density
addParameter(p,'limitElev', 5); % max elevation scan to consider in degree
addParameter(p,'minRange', 5000); % minimium range to consider from the radar in m
addParameter(p,'maxAlt', 3000); % max altitude to consider in meter
addParameter(p,'data_res', 2); % resolution of the data to keep as a coefficient of the original resolution. 

parse(p,varargin{:});
parms = p.Results;

rd = cell(size(scans_filepath));

q = parallel.pool.DataQueue;
afterEach(q,@UpdateWaitbar);
i_s_wb = 0;
N = numel(scans_filepath);
wb = waitbar(0,"Clean and load weather radar data: "+ num2str(0) + "/"+num2str(N));

% check that all file exist
for i_s = 1:numel(scans_filepath)
    if(~exist(scans_filepath(i_s),'file'))
        error("scans_filepath does not exist")
    end
end

parfor i_s = 1:numel(scans_filepath)
    rdt = rsl2mat_wrap(scans_filepath(i_s), struct('max_elev', parms.limitElev));

    % get unique sweeps
    rdt.dz.sweeps = unique_elev_sweeps(rdt, 'dz');
    rdt.rh.sweeps = unique_elev_sweeps(rdt, 'rh');

    % limit elevation angle
    rdt.dz.sweeps = rdt.dz.sweeps([rdt.dz.sweeps.elev]< parms.limitElev);
    rdt.rh.sweeps = rdt.rh.sweeps([rdt.rh.sweeps.elev]< parms.limitElev);

    % Removal of rain
    if parms.removeRain=="rho"
        for i_elev=1:numel(rdt.dz.sweeps)
            [min_elev,i_elev_rh] = min(abs(rdt.dz.sweeps(i_elev).elev-vertcat(rdt.rh.sweeps.elev)));
            assert(min_elev<.01)
            tmp = true(size(rdt.dz.sweeps(i_elev).data));
            min_nbins=min(rdt.rh.sweeps(i_elev_rh).nbins,rdt.dz.sweeps(i_elev).nbins);
            tmp(1:min_nbins,:) = rdt.rh.sweeps(i_elev_rh).data(1:min_nbins,:)>.95;
            rdt.dz.sweeps(i_elev).data(tmp)=nan;
        end
    elseif parms.removeRain=="misnet"
        warning("NOT WORKING. OFFSET")
        [preds, classes] = mistnet_polar_radar(rdt);
        for i_elev=1:numel(rdt.dz.sweeps)
            rdt.dz.sweeps(i_elev).data(preds{i_elev}==classes.rain)=nan;
        end
    else
        warning("no rain removal")
    end

    % remove unused field
    rdt = rmfield(rdt,{'sw','dr','ph','rh','spulse','lpulse','constants','params'});

    % create bird density
    rdt.dens = struct();
    rdt.dens.type = "Bird density";
    rdt.dens.sweeps = struct([]);

    % Process each elevation
    for i_elev=1:numel(rdt.dz.sweeps)

        % Extract reflectivity factor data
        [data, sweep_range, sweep_az] = sweep2mat(rdt.dz.sweeps(i_elev));

        % Make a copy for plot
        data0=data;

        if false
            [z, x, y] = sweep2cart(rdt.dz.sweeps(i_elev),100000,1000);
            figure; imagescnan(x, y, z); axis xy equal;
        end

        if parms.removeOuliars>0

            % 1. Threashold high value
            % density of bird too high
            dens = z_to_refl(idb(data))/11; % inverse decibel transform --> bird density [bird/km^3]
            thr_dens = z_to_refl(idb(32))/11; % 1 bird  36000 eta z_to_refl(36000)/11 32 dbz
            data(dens>thr_dens)=nan;


            % 2. Threashold on marginal distribution.
            % grubbs marginal distribution
            % tmp = isoutlier(data(:),'grubbs','ThresholdFactor',1);
            % n = sum(~isnan(data(:)));
            % alpha=.05;
            % tval = tinv(alpha./n/2,repmat(n-2,1,numel(alpha))).^2;
            % G = (n-1)./sqrt(n) .* sqrt(tval./(n-2+tval));
            % thr = mean(data(:),"omitnan")+[-1 1]*G.*std(data(:),"omitnan")

            mdata = median(data(:),'omitnan');
            MAD=-1/(sqrt(2)*erfcinv(3/2))*median(abs(data(:)-mdata),'omitnan');
            thr = mdata + [-1 1]* parms.removeOuliars*MAD;
            id = data<thr(1)| data>thr(2);
            % figure; hold on; histogram(data); histogram(data(id))
            data(id)=nan;

            % 3. Spatial continuity
            % spatial continiuty with moving median
            data(end-1:end,:)=nan;
            dd = repmat(data,1,1,5*5); u=1;
            for i=-2:2
                for j=-2:2
                    dd(:,:,u) = circshift(data,[i j]);
                    u=u+1;
                end
            end
            data_median = median(dd,3,'omitnan');
            data_n = sum(~isnan(dd),3);
            data_median(end-1:end,:)=nan;
            data_n(end-1:end,:)=0;

            MAD = -1/(sqrt(2)*erfcinv(3/2))*median(abs(dd-data_median),3,'omitnan');
            MAD = nanmean(MAD(:));
            thr = data_median + permute([-1 1],[3 1 2]).* 3.*MAD;
            id = data<thr(:,:,1) | data>thr(:,:,2);
            if false
                figure; hold on; axis tight; box on;
                imagesc(data,'AlphaData',~isnan(data));
                [x,y]=ind2sub(size(id),find(id)); plot(y,x,'.r')

                figure; hold on;
                histogram(data)
                histogram(data(id))
            end

            data(id) = nan;

            % Fill
            % data2=fillmissing(data,"linear");
            % remove datapoint which are too far from neighboring value.
            % data2(data_n<9) = nan;
            % figure; imagesc(data2,"AlphaData",~isnan(data2));colorbar; ylim([0 900])
            % data=data2;

            %         d = data;
            %         % d(isnan(d)) = -inf;
            %         dens = z_to_refl(idb(d))/11;
            %         figure; imagesc(dens, "AlphaData",~isnan(dens)); colorbar; ylim([0 900])

        end

        % Set nodata pulse volumes to 0 reflectivity (-inf on decibel scale)
        % data(isnan(data)) = -inf;

        % Convert to reflectivity
        data    = idb(data);             % inverse decibel transform --> reflectivity factor
        data = z_to_refl(data);          % reflectivity factor (mm^6/m^3) --> reflectivity (cm^2/km^3)
        dens = data/11;                  % reflectivity --> bird density [bird/km^3]

        if parms.shape
            AZ = repmat(sweep_az,1,numel(sweep_range))';
            id=~isnan(dens);
            % radius/distance from center of an ellipse
            R = @(a,b,theta,x) a*b./sqrt(a^2*sind(x-theta).^2+b^2*cosd(x-theta).^2);
            % keep these as fixed
            a = 0.5;
            b = 1.5; 
            % fit only the angle
            Rfit = @(theta,x) R(a,b,theta,x);
            mdl = fit(AZ(id),dens(id),fittype(Rfit),"StartPoint",120);
            dens = dens ./ Rfit(mdl.theta,AZ);
        end

        if parms.smoothing
            % circular smoothing along azimuth but not range
            dens(end-1:end,:)=nan;
            dd = repmat(dens,1,1,5*5); u=1;
            for i=-2:2
                for j=-2:2
                    dd(:,:,u) = circshift(dens,[i j]);
                    u=u+1;
                end
            end

            W=fspecial('gaussian',5,1);
            W=reshape(W(:),1,1,[]);

            dens2 = sum(dd.*W,3,'omitnan') ./ sum(~isnan(dd).*W,3);
            dens2(isnan(dens))=nan;

            if false
                figure;
                nexttile; imagesc(dens,'AlphaData',~isnan(dens));
                nexttile; imagesc(dens2,'AlphaData',~isnan(dens2));
            end

            dens_n = sum(~isnan(dd),3);
            dens2(dens_n<5)=nan;

            dens = dens2;
        end

        % Remove datapoint too close or too high altitude
        dens(sweep_range<parms.minRange,:)=nan;
        [ ~, sweep_height ] = slant2ground( sweep_range, rdt.dz.sweeps(i_elev).elev );
        dens(sweep_height>parms.maxAlt,:)=nan;

        if false
            dens0 = z_to_refl(idb(data0))/11;
            figure; tiledlayout(2,3,'TileSpacing','tight','Padding','tight')
            nexttile([1,2]); hold on;
            imagesc(dens0,'AlphaData',~isnan(dens0)); colorbar; axis tight; ylim([0 1000])
            nexttile; histogram(dens0); x_lim=xlim(); set(gca,"YScale",'log')
            nexttile([1,2]); hold on; imagesc(dens,'AlphaData',~isnan(dens));
            id = find(isnan(dens)&~isnan(dens0));
            [tmp_x, tmp_y] = ind2sub(size(dens0),id);
            plot(tmp_y,tmp_x,'.r','MarkerSize',2)
            colorbar; axis tight; ylim([0 1000])
            yline(max(find(sweep_range<parms.limitRange)),'-k')
            nexttile; histogram(dens);% xlim(x_lim);
            % set(gca,"YScale",'log')
        end

        % Downscale
        dens_dns = nan(size(dens));
        dens_dns(1:parms.data_res:end, 1:parms.data_res:end) = dens(1:parms.data_res:end, 1:parms.data_res:end);

        rdt.dens.sweeps(i_elev).data = dens_dns;
        rdt.dens.sweeps(i_elev).az = sweep_az;
        rdt.dens.sweeps(i_elev).range = sweep_range;
        rdt.dens.sweeps(i_elev).elev = rdt.dz.sweeps(i_elev).elev;

        % some additional data to keep useful for e.g. sweep2mat
        rdt.dens.sweeps(i_elev).range_bin1 = rdt.dz.sweeps(i_elev).range_bin1;
        rdt.dens.sweeps(i_elev).gate_size = rdt.dz.sweeps(i_elev).gate_size;
        rdt.dens.sweeps(i_elev).nbins = rdt.dz.sweeps(i_elev).nbins;
        rdt.dens.sweeps(i_elev).azim_v = rdt.dz.sweeps(i_elev).azim_v;
        rdt.dens.sweeps(i_elev).prf = rdt.dz.sweeps(i_elev).prf;
    end
    rdt = rmfield(rdt,{'dz'});
    rd{i_s} = rdt;
    send( q, i_s );
end
function UpdateWaitbar(~)
    waitbar(i_s_wb/N,wb,"Clean and load weather radar data: "+ num2str(i_s_wb) + "/"+num2str(N));
    i_s_wb = i_s_wb + 1;
end
if isvalid(wb)
    close(wb);
end
end

