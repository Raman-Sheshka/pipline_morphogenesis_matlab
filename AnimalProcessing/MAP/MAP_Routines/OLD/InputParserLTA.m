function [res] = InputParserLTA(Q)

if isempty(Q)
    % GUI
    res = 'gui';
else
    switch class(Q)
        case {'char';'string'}
            % KeyWord
            res = 'key';
        case 'cell'
            % multiple boxes or areas
            if isequal(size(Q{1}), [1 2])
                % multiple boxes
                res = 'mb';
            elseif isequal(size(Q{1}), [1 4])
                % multiple areas
                res = 'ma';
            else
                % ?
                res = 'wtf';
            end
        otherwise
            % ?
            res = 'wtf';
    end
end

if strcmp(res,'wtf')
    msg = sprintf(['Wrong input area format, please use folowing :\n' ...
        'Qarea = [ ]                     : Start GUI\n' ...
        'Qarea = keyword                 : Analyse a predifined area\n' ...
        'Qarea = {[x y];[x y];...}       : Analyse corresponding grid box\n' ...
        'Qarea = {[xmin ymin xmax ymax]} : Analyse corresponding area\n']);
    f = errordlg(msg,'LTA Input Error', 'modal');
end

end

