function [df, rd] = load_F_elev(filepath, rd)
% load_F_elev load the F matrix for each elevation angle from files

%% Load filename of all possible df scan
df_file = struct2table(dir(filepath));
df_file.name = string(df_file.name);
df_file=df_file(startsWith(df_file.name, "F_"),:);
M = split(df_file.name,'_');

%% check the scan needed
% Store the 
rd_i_df = cell(size(rd));

dfnew=struct(); u=0;
for i_s=1:numel(rd)
    rd_i_df{i_s} = nan(numel(rd{i_s}.dens.sweeps),1);
    for i_elev=1:numel(rd{i_s}.dens.sweeps)
        % Find the match in the files

         tmp = find("el"+num2str(round(rd{i_s}.dens.sweeps(i_elev).elev*10)) == M(:,2) & "az" + num2str(numel(rd{i_s}.dens.sweeps(i_elev).az)) == M(:,4));
        
        if isempty(tmp)
            tmp_name = "el"+num2str(round(rd{i_s}.dens.sweeps(i_elev).elev*10)) + "_r" + num2str(numel(rd{i_s}.dens.sweeps(i_elev).range))+ "_az" + num2str(numel(rd{i_s}.dens.sweeps(i_elev).az));
            warning("scan missing: "+tmp_name)
            u=u+1;
            dfnew(u).elev = rd{i_s}.dens.sweeps(i_elev).elev;
            [~, dfnew(u).range, dfnew(u).az] = sweep2mat(rd{i_s}.dens.sweeps(i_elev));
            dfnew(u).az = round(dfnew(u).az / 0.25)*.25;
            dfnew(u).name = tmp_name;
        else
            rd_i_df{i_s}(i_elev) = tmp;
        end
    end
end

if u>0
    save('data/scan_pattern_new','dfnew')
    warning("New df pattern neew to be computed. Saved in 'data/scan_pattern_new'")
    return
end

%% Load scan

rd_i_df_unique = unique(vertcat(rd_i_df{:}));

df = cell(size(rd_i_df_unique));
for i_f=1:numel(df)
    f =  df_file(rd_i_df_unique(i_f),:);
    tmp=load([f.folder{:} '/' f.name{:}]);

    df{i_f} = tmp.dfe;
    % df{i_f}.gridsize = gridsize;
    
    df{i_f}.sz_raz = [numel(tmp.dfe.range) numel(tmp.dfe.az)];
    df{i_f}.F_raz_id=sub2ind( df{i_f}.sz_raz , tmp.Felevidc(:,1), tmp.Felevidc(:,2));
    
    df{i_f}.x = tmp.x;
    df{i_f}.y = tmp.y;
    df{i_f}.z = tmp.z;
    df{i_f}.sz_coarse_grid = [numel(tmp.x) numel(tmp.y) numel(tmp.z)];
    
    [df{i_f}.F_gr_y, df{i_f}.F_gr_x, df{i_f}.F_gr_z] = ind2sub(df{i_f}.sz_coarse_grid, tmp.Felevidc(:,3));

    df{i_f}.Fval = tmp.Felevfc;

%     F = sparse(double(df{i_f}.F_raz_id), double(tmp.Felevidc(:,3)), double(tmp.Felevfc), prod(df{i_f}.sz_raz), prod(df{i_f}.sz_coarse_grid));
%     tmp2 = reshape(full(sum(F,1)), df{i_f}.sz_coarse_grid);

end

%% Order

[~,I] = sort(cellfun(@(x) x.elev, df));
df=df(I);
rd_i_df_unique_s=rd_i_df_unique(I);

for i_s=1:numel(rd_i_df)
    [~, rd{i_s}.i_df] = ismember(rd_i_df{i_s}, rd_i_df_unique_s);
end

assert(numel(unique(cellfun(@(x) numel(x.x), df)))==1)
assert(numel(unique(cellfun(@(x) numel(x.y), df)))==1)
assert(numel(unique(cellfun(@(x) numel(x.z), df)))==1)

%% Compute windspeed

RMIN_M  = 5000;
RMAX_M  = 37500;
ZSTEP_M = diff(df{1}.z(1:2));
ZMAX_M  = ZSTEP_M*numel(df{1}.z);

parfor i = 1:numel(rd)
    [ ~, ~, u, v, rmse, nll, cov] = epvvp(rd{i}, ZSTEP_M, RMIN_M, RMAX_M, ZMAX_M);
    rd{i}.uv = [u v]; % write speed (m/s)
end

end