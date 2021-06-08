function filename_out = FilenameFormatter(filename_in)

% Version 1.1
% Boris Guirao
%
% filename_out = FilenameFormatter(filename_in)
%
% first removes last '_' of "filename_in", then crops it to its first
% nmax_char = 15 characters.
%
% NB: "nmax_char" SHOULD NOT BE CHANGED SO AS TO BE ABLE TO USE PREVIOUS BACKUP FILES!!


%% Code %%

%%% maximal number of characters kept:
nmax_char = 15;        

%%% Removes last character, namely "_"
filename_out = substr(filename_in,  0, -1);

%%% Crops name at its first "nmax_char" characters:
filename_out = substr(filename_out, 0, nmax_char);  

%%% Removes last character if it is '_' (1.1):
if max(strfind(filename_out, '_')) == nmax_char                             % yields maximum index where '_' appears and compares it to nmax_char: if equal, it's cut it
    filename_out = substr(filename_out,  0, -1);
end

%% History %%

% 07/09/2011: 1.1
% - prevents name cut ending with '_'. It so, removes it.

% 12/08/2011: creation