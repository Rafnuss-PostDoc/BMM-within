function [preds, classes, probs] = mistnet_polar_radar(radar)
%MISTNET_POLAR Classify sample volumes in polar coordinates using mistnet


% Run mistnet to get predictions for Cartesian volume
[PREDS, PROBS, classes, x, y, elevs] = mistnet( radar );

% Create interpolating function
F = griddedInterpolant({y, x, elevs}, PREDS, 'nearest');

preds = cell(numel(radar.dz.sweeps),1);
probs = cell(numel(radar.dz.sweeps),1);
for i_elev = 1:numel(radar.dz.sweeps)
    % Convert sample volume coordinate matrices to XYZ
    [~, range, az] = sweep2mat(radar.dz.sweeps(i_elev));
    [RANGE, AZ, ELEV] = ndgrid(range, az, radar.dz.sweeps(i_elev).elev);
    [X, Y, ~] = radar2xyz(RANGE, AZ, ELEV);

    % Interpolate predictions onto polar volume
    preds{i_elev} = F(Y, X, ELEV);

    % If class probabilities are requested, also interpolate those, one at a time
    if nargout > 2
        n_classes = size(PROBS, 4);
        probs{i_elev} = zeros([sz n_classes]);
        for c = 1:n_classes
            F = griddedInterpolant({y, x, elevs}, PROBS(:,:,:,c), 'nearest');
            probs{i_elev}(:,:,:,c) = F(Y, X, ELEV);
        end    
    end
end
end

