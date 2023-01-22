function [ fileinfo ] = aws_list_datetime( station, time )
%AWS_LIST_DATETIME Get a list of archive files
%
% [ files ] = aws_list_datetime( station, time )
%

y = year(time);
m = month(time);
d = day(time);
h = hour(time);


% If hour is specified, append timestamp to narrow selection to selected
% hour
if h > 0
    [ fileinfo ] = aws_list( station, y, m, d, h );
else
    [ fileinfo ] = aws_list( station, y, m, d );
end

end
