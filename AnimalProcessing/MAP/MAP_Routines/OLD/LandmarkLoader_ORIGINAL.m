% LandmarkLoader
%
% Function that go search for the macrocaete coordinate file in the folder structure
% It uses regular expression in order to work on multiple version of the path
% The path still need to contain specific keywords, see the regular expression,
% and must have the structure of raw / registration / landmarks / macro.txt
% 
% Stephane Rigaud
% version 1.0
%

function [macroCoord] = LandmarkLoader(pathFolderRaw, animal, sideShort)

if strcmp(sideShort,'r'),     sideLong = 'RIGHT';
elseif strcmp(sideShort,'b'), sideLong = 'BOTH' ;
else                          sideLong = 'LEFT' ;
end

%%% Regular Expression ----------------------------------------------------
regexpRegistrationFolder = ['\w*Spacetime\w*Registration'];
regexpLandmarksFolder    = [animal '\w*Landmarks\w*' sideLong];
regexpMacroFile          = ['\w*Macro\w*[Cc]oordinate\w*' animal];
%--------------------------------------------------------------------------

%%% Registration Folder ---------------------------------------------------
dirRawFolder = dir(pathFolderRaw);
dirRawFolder(strncmp({dirRawFolder.name}, '.', 1)) = [];
ind = cellfun(@(x)( ~isempty(x) ), regexpi({dirRawFolder.name}, regexpRegistrationFolder));
dirRawFolder(~ind) = [];
if length(dirRawFolder) ~= 1
    disp('Error LandmarkLoader: could not find the Registration Folder!');
    return
else
    RegistrationFolder = dirRawFolder.name;
end
%--------------------------------------------------------------------------

%%% Landmarks Folder ------------------------------------------------------
dirRegistrationFolder = dir([pathFolderRaw filesep RegistrationFolder]);
dirRegistrationFolder(strncmp({dirRegistrationFolder.name}, '.', 1)) = [];
ind = cellfun(@(x)( ~isempty(x) ), regexpi({dirRegistrationFolder.name}, regexpLandmarksFolder));
dirRegistrationFolder(~ind) = [];
if length(dirRegistrationFolder) ~= 1
    disp('Error LandmarkLoader: could not find the Landmarks Folder!');
    return
else
    LandmarksFolder = dirRegistrationFolder.name;
end
%--------------------------------------------------------------------------

%%% Macro File ------------------------------------------------------------
dirLandmarksFolder = dir([pathFolderRaw filesep RegistrationFolder filesep LandmarksFolder]);
dirLandmarksFolder(strncmp({dirLandmarksFolder.name}, '.', 1)) = [];
ind = cellfun(@(x)( ~isempty(x) ), regexpi({dirLandmarksFolder.name}, regexpMacroFile));
dirLandmarksFolder(~ind) = [];
if length(dirLandmarksFolder) ~= 1
    disp('Error LandmarkLoader: could not find the Coordinate File!');
    return
else
    MacroFile = dirLandmarksFolder.name;
end
%--------------------------------------------------------------------------

%%% lock and load ---------------------------------------------------------
macroCoordPath = [pathFolderRaw filesep RegistrationFolder filesep LandmarksFolder filesep MacroFile];
macroCoord = dlmread(macroCoordPath);
%--------------------------------------------------------------------------
end

%% Historic

% 09/05/2016: v1.0


