function UnionsegMaker()

%% Initialisation
S1params = Stage0_Initial_info();

% digit format
digits= ['%0' num2str(S1params.DIGITNUMBER) 'd'];
% create directory if necessary
if ~isdir([S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_backup']), mkdir(S1params.PATHFOLDER,[S1params.OUTPUTNAME '_backup']), end;

% clean folders
% remove old suspects files
delete([S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_cellout' filesep 'suspects-*']);
% delete keyframes memory
delete([S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_backup' filesep 'keyframes.*']); 
% delete old modification masks
delete([S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_backup' filesep 'modif_mask_' S1params.ROOTFILENAME '*']);
% delete old correction masks
delete([S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_backup' filesep 'correction_mask_' S1params.ROOTFILENAME '*']);  
% delete old segmentation boxes
if exist([S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_backup' filesep S1params.ROOTFILENAME 'boxes.mat'],'file') > 0;
    delete([S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_backup' filesep S1params.ROOTFILENAME 'boxes.mat']);
end

%% PREPROCESSING %%
disp('Creating Unionseg ...');
nFrames = S1params.LASTIMAGE - S1params.FIRSTIMAGE;
parfor_progress(nFrames);
parfor i = S1params.FIRSTIMAGE:S1params.LASTIMAGE
    
    % Loading of raw image
    image = imread([S1params.PATHFOLDER filesep S1params.ROOTFILENAME num2str(i,digits) '.' S1params.FILEEXTENSION]);
    % union of selected segmentations
    union = zeros(size(image));
    union = union + ~imread([S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_results' filesep 'directskel_' S1params.ROOTFILENAME num2str(i,digits) '.png']);
    % binarize
    union = im2bw(union);   
    % remove small cells
    union = ~SmallCellRemover(~union, S1params.FILLSMALLCELLS);
    % Removes incomplete sides
    union = ~watershed(union,4);
    % save result    
    imwrite(~union,[S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_results' filesep 'Unionseg_' S1params.ROOTFILENAME num2str(i,digits) '.png'],'png');    
    parfor_progress;
end
parfor_progress(0);


end % end of Stage

%% History %%
 
% 2015-10-07: 3
% - Code cleaning and variable renaming

% 2012-09-18: 2.5
% - code cleaning done.

% 2012-07-20: 2.4.2
% - nodisplay option added.

% 2012-07-17: 2.4.1
% - proprer unionseg cleaning process added.

% 2011-09-15: 2.4
% - previously made correction mask & modification deletion added.
% - minor change on starting warning message.

% 2011-09-14: 2.3
% - previously made suspects file and keyframes memory deletion added.

% 2011-08-25: 2.2
% - input parameters loading changed.

% 2011-07-06: 2.1
% - version number removed from filename.
% - Initial input adjusted.

% 22/09/2010: 2.0
% - removed saving of Unionseg_ backups
% - added warning message before overwriting Unionseg_ images

% 18-20/09/2010: 2.0GM
% - Removed all Benoit files and folder
% - use of Small_Cell_Remover
% - added display of Unionseg_ images
% - added questdlg at the beginning checking that >1 segmentation was
% selected.
% - added series of  "if threedseg == 1...end" at the end so as not to save
% meaningless images


