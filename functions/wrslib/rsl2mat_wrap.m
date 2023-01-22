function radar = rsl2mat_wrap(filepath, params)
% Wrapper of rsl2mat adding a few improvements

% Get the radar name (callid from the filepath
[~,filename]=fileparts(filepath);
callid = char(filename);

% Ensure that filepath is a char
radar = rsl2mat(char(filepath), callid(1:4), params);

% Convert year, month, day etc.. to datetime and remove them
radar.datetime = datetime(radar.year,radar.month,radar.day,radar.hour,radar.minute,radar.second);
radar = rmfield(radar,{'year','month','day','hour','minute','second'});
end