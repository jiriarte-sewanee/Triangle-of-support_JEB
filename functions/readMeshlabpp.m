% read meshlab .pp file
% return the x y z data as n x 3 array
function output = readMeshlabpp(fileName)
    % feature_num x 3 array
    data = [];
    names = {} ;
    fid = fopen(fileName); 
    line = fgetl(fid);
    while ischar(line)
        indices = strfind(line, '<point');
        if length(indices) ~= 1 || indices ~= 2
%             fprintf("%d", indices);
            line = fgetl(fid);
            continue;
        end
        indices = strfind(line, '"'); 
        assert(mod(length(indices), 2) == 0); % check that " comes in pairs
        for i = 1 : length(indices) / 2
            % quote positions
            pos1 = indices(2 * i - 1);
            pos2 = indices(2 * i);
            
            if strcmp(line(pos1-3:pos1-2), ' x')
                xstr = line(pos1 + 1 : pos2 - 1);
            elseif strcmp(line(pos1-3:pos1-2), ' y')
                ystr = line(pos1 + 1 : pos2 - 1);
            elseif strcmp(line(pos1-3:pos1-2), ' z')
                zstr = line(pos1 + 1 : pos2 - 1);
            elseif strcmp(line(pos1-3:pos1-2), 'me')
                nstr = line(pos1+1 : pos2-1) ;
            end
        end
%         fprintf("%s %s %s\n", xstr, ystr, zstr);
%         fprintf('%s\n', line);
        x = str2double(xstr);
        y = str2double(ystr);
        z = str2double(zstr);
        data = [data; [x y z]];
        names = [names;{nstr}] ;
        line = fgetl(fid);
    end
    output = [names num2cell(data)] ;
    fclose(fid);
end