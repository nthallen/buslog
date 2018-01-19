function output_txt = buslog_cursor_text_func(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
output_txt = {['Y: ',num2str(pos(2),8)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end

output_txt{end+1} = ['X: ',num2str(pos(1),8)];

[Hi,Mi,Si] = time2hms(pos(1)*3600);
Tsec = floor(Si/10);
Osec = Si - Tsec*10;
output_txt{end+1} = sprintf('T: %02d:%02d:%d%g\n', Hi, Mi, Tsec, Osec);
