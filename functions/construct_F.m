function [F, data, info] = construct_F(rd, df, g)%, varargin)


%% Extract radar scan 
% A radar scan contains all different scan (elevation) possible for all
% timestep for each radar
rs = cell(height(rd),1);
for i_r = 1:height(rd)
    tmp = cellfun(@(x) [x.i_df], rd(i_r,:), 'UniformOutput', false);
    rs{i_r}.i_df = unique(vertcat(tmp{:}));
    rs{i_r}.lon = rd{i_r,1}.lon;
    rs{i_r}.lat = rd{i_r,1}.lat;
    rs{i_r}.station = rd{i_r,1}.station;

    % Find location on grid
    [~,rs{i_r}.lon_id]=min(abs(g.lon-rs{i_r}.lon));
    [~,rs{i_r}.lat_id]=min(abs(g.lat-rs{i_r}.lat));

    % Check that the radar is within the grid
%     id_lon = rs{i_r}.lon_id - ceil(df.sz_coarse_grid(1)/2);
%     id_lat = rs{i_r}.lat_id - ceil(df.sz_coarse_grid(2)/2);
%     assert(id_lon>=0 & (id_lon+df.sz_coarse_grid(1))<=numel(g.lon))
%     assert(id_lat>=0 & (id_lat+df.sz_coarse_grid(2))<=numel(g.lat))
end

%% 
% Take the larger size of az and range to have everything on a regular grid
tmp = cellfun(@(x) x.sz_raz,df,'UniformOutput',false);
sz_rd = [max(vertcat(tmp{:})) numel(rs) numel(df)];

% Compute a index of possible value
id_start = reshape(cumsum(repmat(prod(sz_rd(1:2)), prod(sz_rd(3:4)),1)), sz_rd(3), sz_rd(4))-prod(sz_rd(1:2));

% Check if id_start is correct [id_range, id_az, id_r, id_df]=ind2sub(sz_rd,id_start+1);

% We also want to convert the index of F_raz_id into the new grid
% (e.g. 900x360 -> 1800 x 720)
df_i_f_F_raz_id = cell(size(df));
for i_f=1:numel(df)  
    [row,col] = ind2sub(df{i_f}.sz_raz, df{i_f}.F_raz_id);
    df_i_f_F_raz_id{i_f} = sub2ind(sz_rd(1:2), row, col);
end

%% Compute data
data_C = cell(numel(rs), numel(df));
data_pres = cell(numel(rs), numel(df));
for i_r = 1:numel(rs) % through the radars
    for i_f = 1:numel(df) % through the unique scan
        % build the data matrix for this elevation scan for all time. 
        data_tmp=nan([sz_rd(1:2) width(rd)]);
        for i_t = 1:width(rd) % time 
            % find the index in the radar data of the particular scan
            i_elev = find(rd{i_r,i_t}.i_df==i_f);

            if isempty(i_elev) % Does not exist in this radar/time
                continue
            end
            
            % Check that the first few range azimuth are similar (i.e., similar
            % start and same resolution) -> assess as the difference to the
            % reoslution.
            assert(abs(df{i_f}.elev-rd{i_r,i_t}.dens.sweeps(i_elev).elev)<0.1)
            assert(all(abs(df{i_f}.range(1:10)-rd{i_r,i_t}.dens.sweeps(i_elev).range(1:10))<diff(df{i_f}.range(1:2))))
            assert(all(abs(df{i_f}.az(1:10)-rd{i_r,i_t}.dens.sweeps(i_elev).az(1:10))<diff(df{i_f}.az(1:2))))
            
            % take the data on the full data grid (1820x720) instead of the
            % original size tailored to each df.
            min_raz = size(rd{i_r,i_t}.dens.sweeps(i_elev).data);
            data_tmp(1:min_raz(1),1:min_raz(2),i_t) = rd{i_r,i_t}.dens.sweeps(i_elev).data+eps; % trick to use sparse for nan. 

        end
  
        % Find all the nan (no data available)
        tmp2 = find(~isnan(data_tmp));
        % and convert their position to their index 
        [i_range,i_az,i_time] = ind2sub(size(data_tmp),tmp2);
        % convert range and azimuth in raz index
        i_raz=sub2ind(sz_rd(1:2),i_range,i_az);

        % check
        [id_range, id_az, id_r, id_df]=ind2sub(sz_rd,id_start(i_r,i_f) + i_raz);
        assert(all(id_r==i_r))
        assert(all(id_df==i_f))
        assert(all(max(id_range)<=numel(df{i_f}.range)))
        assert(all(max(id_az)<=numel(df{i_f}.az)))

        data_C{i_r,i_f} = [id_start(i_r,i_f)+i_raz i_time data_tmp(tmp2)];

        % For computing F, we need to know if there is at least one data
        % present in any timestep
        data_pres{i_r,i_f} = find(any(~isnan(data_tmp),3));
    end
end

tmp = vertcat(data_C{:});
data = sparse(tmp(:,1), tmp(:,2), tmp(:,3), prod(sz_rd), width(rd));


%% Compute F

F_build = cell(numel(rs), numel(df));

for i_r = 1:numel(rs) % through the radars

    % Convert the index of position in x,y,z on the radar volume to the
    % position (index) of the multi-radar grid. This is done by offest
    % the lat lon sub by the position of the radar. 
     % get the position of the 0,0 of the coarse grid
    id_lon = rs{i_r}.lon_id - ceil(df{1}.sz_coarse_grid(1)/2);
    id_lat = rs{i_r}.lat_id - ceil(df{1}.sz_coarse_grid(2)/2);
    
    for i_f = 1:numel(df) % through the unique scan

        % 
        id_pres_data = ismember(df_i_f_F_raz_id{i_f},data_pres{i_r,i_f});
  

        % Assemble all the df of this scan/radar into a table to be able to
        % filter
        dfe = table( ...
            df{i_f}.F_gr_x(id_pres_data), ...
            df{i_f}.F_gr_y(id_pres_data), ...
            df{i_f}.F_gr_z(id_pres_data), ...
            ... %df.F_raz_id{i_f}(id_pres_data), ... % index in raz coordinate
            id_start(i_r,i_f) + df_i_f_F_raz_id{i_f}(id_pres_data),... % indes in full data coordinate (sz_d)
            df{i_f}.Fval(id_pres_data), ...
            'VariableNames',{'x','y','z','id_d','Fval'});

        % Find location of the grid 
        id_g_xyz = [dfe.x+id_lon dfe.y+id_lat dfe.z];
        % Filter for data outside the grid
        id_F_keep = id_g_xyz(:,1)>=1 & id_g_xyz(:,1)<=g.sz(1) & id_g_xyz(:,2)>=1 & id_g_xyz(:,2)<=g.sz(2);
        if mean(id_F_keep)<1
            warning("Some part of radar is not covered by the grid")
        end
        dfe = dfe(id_F_keep,:);
        % compute the index of the grid
        dfe.id_g = sub2ind(g.sz, id_g_xyz(id_F_keep,1), id_g_xyz(id_F_keep,2), id_g_xyz(id_F_keep,3));
        
        F_build{i_r,i_f} = removevars(dfe,{'x', 'y', 'z'});

%         F = sparse(dfe.id_d, dfe.id_g, double(dfe.Fval), prod(sz_rd), prod(g.sz));
%         tmp = reshape(full(sum(F,1)), g.sz);
%         figure;imagesc(log(sum(tmp,3)))
    end
end

%% 
tmp = vertcat(F_build{:});
F = sparse(tmp.id_d, tmp.id_g, double(tmp.Fval), prod(sz_rd), prod(g.sz));

%%
% norm = sum(F,2);
% % tmp = find(Fsp>0);
% % [row,~] = ind2sub(size(Fsp),tmp);
% % Fsp(tmp) = Fsp(tmp)./ norm(row);
% [row, col] = find(F>0);
% tmp = F(F>0) ./ norm(row);
% F = sparse(row,col,tmp,size(F,1),size(F,2));

%% 
info.sz_rd = sz_rd;
info.id_start = id_start;
end