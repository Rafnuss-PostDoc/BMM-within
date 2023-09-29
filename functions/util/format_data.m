function [D, id_nnz, A, E, R, H] = format_data(data, df, info)

    % retrieve id of non-zeros
    id_nnz = find(data(:,1));
    
    % keep only non-zeros
    D = full(data(id_nnz,1));

    if nargout>2
    
        % find the index of the data in the dimension of radar, df, range and
        % azimuth for the non-zeron data
        [id_range, id_az, id_r, id_df] = ind2sub(info.sz_rd, id_nnz);
        
        df_range = cellfun(@(x) x.range, df,'UniformOutput',false);
        tmp = df_range(id_df);
        R = cellfun(@(x,y) x(y),tmp, num2cell(id_range));
        
        df_az = cellfun(@(x) x.az, df,'UniformOutput',false);
        tmp = df_az(id_df);
        A = cellfun(@(x,y) x(y), tmp, num2cell(id_az));
        
        df_e = cellfun(@(x) x.elev, df);
        E = df_e(id_df);

    end

    if nargout>5
        [ ~, H ] = slant2ground( R(:), E(:) );
    end
end