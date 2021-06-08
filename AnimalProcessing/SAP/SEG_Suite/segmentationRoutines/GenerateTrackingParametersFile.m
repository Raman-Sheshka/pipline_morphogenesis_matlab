function generateTrackingParametersFile(filefullpath,parametersSet)
% GENERATETRACKINGPARAMETERSFILE generates the parameters file used by
% the C++ tracking program (v3.4- and higher)
%
% USAGE: GENERATETRACKINGPARAMETERSFILE(FILEFULLPATH,PARAMETERSSET) with
% × FILEFULLPATH the fullpath to the file to generate.
% × PARAMETERSSET the structure containing the parameters value.
%
% @author GOYA Yûki
% @version 0.1 (2012-06-25)

    % open a file with write permission (discard previous content if any)
    fileID = fopen(filefullpath,'w');
    % write SPATH
    fprintf(fileID,'SPATH %s\n',parametersSet.SPATH);
    % write OPATH
    fprintf(fileID,'OPATH %s\n',parametersSet.OPATH);
    % write FIRST
    fprintf(fileID,'FIRST %d\n',parametersSet.FIRST);
    % write LAST
    fprintf(fileID,'LAST %d\n',parametersSet.LAST);
    % write PPATH
    fprintf(fileID,'PPATH %s\n',parametersSet.PPATH);
    % write COPT if specified
    if isfield(parametersSet,'COPT') && ~isempty(parametersSet.COPT)
        fprintf(fileID,'COPT %d\n',parametersSet.COPT);
    end
    % write OOPT if specified
    if isfield(parametersSet,'OOPT') && ~isempty(parametersSet.OOPT)
        fprintf(fileID,'OOPT %s\n',parametersSet.OOPT);
    end
    % write MDR2 if specified
    if isfield(parametersSet,'MDR2') && ~isempty(parametersSet.MDR2)
        fprintf(fileID,'MDR2 %0.2f\n',parametersSet.MDR2);
    end
    % write MDR3 if specified
    if isfield(parametersSet,'MDR3') && ~isempty(parametersSet.MDR3)
        fprintf(fileID,'MDR3 %0.2f\n',parametersSet.MDR3);
    end
    % write DDR2 if specified
    if isfield(parametersSet,'DDR2') && ~isempty(parametersSet.DDR2)
        fprintf(fileID,'DDR2 %0.2f\n',parametersSet.DDR2);
    end
    % write DDR3 if specified
    if isfield(parametersSet,'DDR3') && ~isempty(parametersSet.DDR3)
        fprintf(fileID,'DDR3 %0.2f\n',parametersSet.DDR3);
    end
    % write ANISO if specified
    if isfield(parametersSet,'ANISO') && ~isempty(parametersSet.ANISO)
        fprintf(fileID,'ANISO %0.2f\n',parametersSet.ANISO);
    end
    % write APOP if specified
    if isfield(parametersSet,'APOP') && ~isempty(parametersSet.APOP)
        fprintf(fileID,'APOP %0.2f\n',parametersSet.APOP);
    end
    % write NAPO if specified
    if isfield(parametersSet,'NAPO') && ~isempty(parametersSet.NAPO)
        fprintf(fileID,'NAPO %d\n',parametersSet.NAPO);
    end
    % close file
    fclose(fileID);
end

% --- HISTORY ---
% 0.1 (2012-06-25)
% * creation
% ---------------