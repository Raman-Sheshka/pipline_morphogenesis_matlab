%% Initialisation
S1params = Stage0_Initial_info();

digits = ['%0' num2str(S1params.DIGITNUMBER) 'd'];
% creates directory if necessary
if ~isdir([S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_results_alt']), mkdir(S1params.PATHFOLDER,[S1params.OUTPUTNAME '_results_alt']),end;
if ~isdir([S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_results_final']), mkdir(S1params.PATHFOLDER,[S1params.OUTPUTNAME '_results_final']),end;

% process parameters
preprocessedTag = '';
if S1params.PREPROCESSED
    preprocessedTag = [S1params.OUTPUTNAME '_preprocessed' filesep 'preprocessed_'];
end
connectivityExtendedMinima = 8;
connectivityWatershed = 8;

prompt = {'Enter gaussian sigma (between [1 3]):','Enter minimum intensity (between [5 30] for 8b):'};
dlg_title = 'Input';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);

S1params.PARA_GAUSSIAN = str2num(answer{1});
S1params.PARA_EXTMIN = str2num(answer{2});

disp('Calculating alternative segmentation ...');

%% PROCESSING
filterGauss = fspecial('gaussian',S1params.PARA_GAUSSIAN,S1params.PARA_GAUSSIAN/4);
% image loop
nFrames = S1params.LASTIMAGE-S1params.FIRSTIMAGE+1;
parfor_progress(nFrames);
parfor i = S1params.FIRSTIMAGE:S1params.LASTIMAGE
    % read image
    image = imread([S1params.PATHFOLDER filesep preprocessedTag S1params.ROOTFILENAME num2str(i,digits) '.' S1params.FILEEXTENSION]);
    % apply gaussian filter
    bluredImage = imfilter(double(image),filterGauss,'replicate');
    % search the region of extended minima
    extendedMinima = imextendedmin(bluredImage,S1params.PARA_EXTMIN,connectivityExtendedMinima);
    % creates an image that forces the previous regions to be minima
    minimaI = imimposemin(bluredImage,extendedMinima);
    % simple watershed on the image with seeded minima
    segmentedImage = ~watershed(minimaI,connectivityWatershed);
    % Skeletonize
    segmentedImage = ~bwmorph(segmentedImage,'skel',Inf);
    % Suppress small and big cells
    segmentedImage = SmallCellRemover(segmentedImage,S1params.FILLSMALLCELLS);
	% Suppress border junction
    if S1params.BORDERJUNCTIONS
        segmentedImage= Border_Junction_Remover(segmentedImage);
    end
    % Save result
    imwrite(segmentedImage,[S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_results_alt' filesep 'directskel_' S1params.ROOTFILENAME num2str(i,digits) '.png'],'png');
    % display
    parfor_progress;
    % disp(['Segmentation of frame #' num2str(i,digits) ' ... Done!']);
end
parfor_progress(0);

disp('Creating Unionseg ...');
parfor_progress(nFrames);
parfor i = S1params.FIRSTIMAGE:S1params.LASTIMAGE
    
    % Loading of raw image
    image = imread([S1params.PATHFOLDER filesep S1params.ROOTFILENAME num2str(i,digits) '.' S1params.FILEEXTENSION]);
    % union of selected segmentations
    union = zeros(size(image));
    union = union + ~imread([S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_results_alt' filesep 'directskel_' S1params.ROOTFILENAME num2str(i,digits) '.png']);
    % binarize
    union = im2bw(union);   
    % remove small cells
    union = ~SmallCellRemover(~union, S1params.FILLSMALLCELLS);
    % Removes incomplete sides
    union = ~watershed(union,4);
    % save result    
    imwrite(~union,[S1params.PATHFOLDER filesep S1params.OUTPUTNAME '_results_alt' filesep 'Unionseg_' S1params.ROOTFILENAME num2str(i,digits) '.png'],'png');    
    parfor_progress;
end
parfor_progress(0);

%% History %%

% 2015-10-07: 3
% - Code cleaning and variable renaming
% - Add preprocessing option

% 2014-06-12: 2.5
% - Suppress border junction option

% 2012-09-17: 2.4
% - input arguments changed.

% 2012-09-05: 2.3
% - use of S1.display option added.

% 2011-09-26: 2.2.3
% - adaptative histogram equalization removed.

% 2011-09-15: 2.2.2
% - adaptative histogram equalization corrected.
% - figure closing added.

% 2011-09-09: 2.2.1
% - input parameters loading adjusted.
% - adaptative histogram equalization added.

% 2011-08-25: 2.2
% - input parameters loading adjusted.

% 2011-07-06: 2.1
% - Histogram equalizaton removed.
% - version number removed from filename.
% - Initial input adjusted.

% 2011-06-15: 2.01
% - Histogram equalizaton added.

% 22/09/2010: 2.0

% 16/09/2010: 2.0GM
% - Removed all Benoit files and folder
% - use of function "imoverlay" to created 8bits RGB image "directskel"
% - added color of skeleton as parameter
% - made image display in full screen

% 01/06/09 : change to be like other segmentations functions

% 10/05/09 : add benoit and normal output

% 20/04/09 : first version
% simple iteration with min_watershed without using previous skeleton but
% only with parameters para_gaussian and para_extmin