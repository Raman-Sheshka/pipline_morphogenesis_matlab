%% PIV %%

% Version 2.2
% Boris Guirao
% Yûki Goya
% Anaïs Bailles
% last update: 2014-02-26



%% Additional Information, Notations and creation of directories %%

Version_name = '2.2';
Program = 'PIV';                                                             % This program name. Used for saving files
Movie = ['_' Animal];
image_extension_input_Raw=['.' image_format_input_Raw];
image_extension_output_PIV=['.' image_format_output_PIV];
print_format=['-d' image_format_output_PIV];
print_resolution_PIV=['-r' num2str(Resolution_PIV)];

% CHECKS THAT AN APPROPRIATE GRIDSIZE WAS ENTERED (changed 2.0GMc):
if strcmp(PIV_Grid,'XS')
    grid_size = 1;
elseif strcmp(PIV_Grid,'S')
    grid_size = 2;
elseif strcmp(PIV_Grid,'M')
    grid_size = 4;
elseif strcmp(PIV_Grid,'L')
    grid_size = 8;
elseif strcmp(PIV_Grid,'XL')
    grid_size = 16;
elseif strcmp(PIV_Grid,'XXL')
    grid_size = 32;
else
    disp('Error: Please choose a proper grid size: XS, S, M, L, XL, XXL')
    return
end


%% Creation of directories %%

%%% Creates PIV SUBfolder:
Path_save_SUBfolder = [Path_save_folder filesep Program  '_' Animal];
mkdir(Path_save_SUBfolder);

%%% Grid size folder:
Grid_size = ['Grid_' PIV_Grid];
Gridsize_folder = [Path_save_SUBfolder filesep Grid_size];
mkdir(Gridsize_folder);

%%% Backup folder:
Backup_folder = [Path_save_SUBfolder filesep Grid_size filesep 'Backups'];
mkdir(Backup_folder)

%%% Frame folder:
if Velocity_display==1 || Speed_display==1 || Divergence_display==1
    Frame_folder = [Path_save_SUBfolder filesep Grid_size filesep 'Frames'];
    mkdir(Frame_folder);
end

%%% Saving txt file to specify program version:
Today = datestr(now,29);                   % format 29 displays date yyyy-mm-dd style. Look up date for details
Today_Full_Start = datestr(now);
Info = ['Program "' Program '" version ' Version_name ' with grid size ' PIV_Grid  '  Run on animal "' Animal '". Started on ' Today_Full_Start ];
dlmwrite([Gridsize_folder filesep Today '_' Program '_' Version_name '_' PIV_Grid  '_' Animal '.txt'], Info, 'delimiter', '', 'newline','pc')


%% Iteration over frames

%%% the grid size is specified here
if grid_size<=16
    size_grid_1 = 16*grid_size;
    size_grid_2 = 8*grid_size;
else                                                                       % for XXL grid_size
    size_grid_1 = 8*grid_size;
    size_grid_2 = 8*grid_size;
end
arrow_color = 'm';                                                          % changed arrow color to magenta (2.2)

%%% Loads raw images. All raw images assumed to be in same folder (1.3):
raw_seg_image = imread([Path_folder_Raw filesep Filename_Raw num2str(Start_frame,digits_format) image_extension_input_Raw]);               % selects 1st raw image (2.0GMa)
image_size = size(raw_seg_image);
% 2011-10-13: set WIDTH & HEIGHT
HEIGHT = image_size(1);
WIDTH = image_size(2);

if Figure_Replot_PIV == 0                                                  % checks whether it is not just a graphics plot (2.0a)

    tic
    %%% Parfor Progress bar initialization:
    N_Vfields = Final_frame-Start_frame;                                                                                    % number of velocity fields calculated
%     ppm = ParforProgMon(['Processing PIV ' Version_name ' with Grid ' PIV_Grid ' for "' Animal '" images. Please wait...'], N_Vfields);           % 2.0GMb
    
    parfor n = Start_frame:Final_frame-1
        
        %% PIV processing %%
        
        % All raw images assumed to be in same folder (1.3):
        disp(' ');
        disp(['Applying "matpiv" on image # ' num2str(n) '...']); % 2.2
        im1 = imread([Path_folder_Raw filesep Filename_Raw num2str(n,digits_format) image_extension_input_Raw]);
        im2 = imread([Path_folder_Raw filesep Filename_Raw num2str(n+1,digits_format) image_extension_input_Raw]);
        
        [x,y,u,v,snr,pkh] = matpiv(im1,im2,...      % selects 1st raw image
                                   [size_grid_1 size_grid_1; size_grid_2 size_grid_2],1,0.5,'multin'); %#ok<PFBNS,PFTUS>
        
%         [x,y,u,v,snr,pkh] = matpiv([Path_folder_Raw,filesep,Filename_Raw{1} num2str(n,digits_format) image_extension_input_Raw],...        % selects 1st raw image
%                                    [Path_folder_Raw,filesep,Filename_Raw{1} num2str(n+1,digits_format) image_extension_input_Raw],...      % selects 1st raw image
%                                    [size_grid_1 size_grid_1; size_grid_2 size_grid_2],1,0.5,'multin'); %#ok<PFBNS,PFTUS>
        % ... [s1 s1 ; s2 s2] =gridsize; 1 =timescale; 0.5 =overlap; 'multin' =several iterations       
        
        
        %% Applying filters (edit 2.2)%%

        [u v] = InfFilter(u,v,[WIDTH HEIGHT]);        % removes "Inf" values from displacement field (2011-10-14)
        
        disp(' ');
        disp('Applying "snrfilt" filter based on signal to noise ratio: putting NaNs at locations corresponding to areas of image deemed irrelevant (snr<1.3)...');
        disp(' ');
        [u,v] = snrfilt(x,y,u,v,snr,1.3);
        % NB: this filter is pretty good at finding the borders of the animal!
        
        disp('Applying "peakfilt" filter based on normalized height of correlation peaks: putting NaNs at locations where a weak correlations was found (pkh<0.3)...');
        disp(' ');
        [u,v] = peakfilt(x,y,u,v,pkh,0.3);          % peak height filter
        disp(' ');
        
        %[u,v] = globfilt(x,y,u,v,3);                % global filter (based on std)(commented 2.2)
        % NB: commented in 2.2 because in fullt thorax, velocity fields can be very heterogeneous. Applying not-so-local filter instead
        
        disp(' ');
        disp('Applying local filter "localfilt" replacing displacement values outside [median (over 5²-1 neighbors) +/- 4*sigma] by NaNs...');
        disp(' ');
        [u,v] = localfilt(x,y,u,v,4,'median',5);    % local filter (based on 'median' or 'mean') (changed parameters 2.2)
        % OLD:
        %[u,v] = localfilt(x,y,u,v,2,'median',3);    % local filter (based on 'median' or 'mean')

        % Loading image to determine regions with no signal (0 intensity)(2011-10-14)
%         im = imread([Path_folder_Raw filesep Filename_Raw num2str(n,digits_format) image_extension_input_Raw]);
        
        % SPLITTING VELOCITY FIELDS to be used for segmentation/tracking (u,v) and for quantitative analysis (uq,vq) (2.2):
        % quantitative analysis: (uq,vq)
        %-------------------------------------------------------------------------------------------
        % Put NaNs at xy gridpoints whose surounding boxes (of size "size_grid_2") only contain pixels with 0 gray level.
        [uq vq] = SignalFilterNan(x,y,u,v,im1,size_grid_2); % changed size of the box used to cut signal closer to black regions (to size_grid_2 from 2*size_grid_2)
        %-------------------------------------------------------------------------------------------
        
        % segmentation/tracking: (u,v)
        %-------------------------------------------------------------------------------------------
        % Put 0s at xy gridpoints whose surounding boxes (of size "size_grid_2") only contain pixels with 0 gray level.
        [u v] = SignalFilter(x,y,u,v,im1,size_grid_2); % changed size of the box used to cut signal closer to black regions (to size_grid_2 from 2*size_grid_2)
        % OLD:
        % [u v] = SignalFilter(x,y,u,v,im,2*size_grid_2);
        [u,v] = naninterp(u,v,'linear');                    % interpolate NaN from the field
        %-------------------------------------------------------------------------------------------
        
        %%% Saving x and y for start frame ONLY 
        if n == Start_frame
            display('Saving "XY_Grid"');
            ParSavePIVxy(Backup_folder, [Program  '_' PIV_Grid  '_' Animal '_Backup_XY_Grid'], x, y)        
        end
        %%% Saving "u" and "v" for each interframe WITHIN THE PARFOR LOOP (1.4a):
        display(['Saving Backup file #' num2str(n,digits_format)]);
        ParSavePIVuv(Backup_folder, [Program  '_' PIV_Grid  '_' Animal '_Backup_'], digits_format, n, u, v, uq, vq); % 2.2
        % parSave_PIV_uv(Backup_folder, [Program  '_' PIV_Grid  '_' Animal '_Backup_'], digits_format, n, u, v)
        
        %%% Progress display:
%         ppm.increment();                                                       %#ok<PFBNS>
    end
    loop_one = toc;
    disp(['Loop #1 duration: ' num2str(loop_one) ' seconds'])
end


%% Display Vector field, Speed map, Divergence map %%

% Moved out of parfor loop in 1.4a beacause don't think it's possible to save images in parallel computing!!

if Velocity_display == 1 || Speed_display == 1 || Divergence_display == 1
    
    %%% progressbar initialization:
    progressbar('PIV display and saving...')                               % use progressbar (2.0GMa)

    tic
    for n = Start_frame:Final_frame-1
        
        When = datestr(now);
        
        %% Loading x,y,u,v from backups %%
        
        %%% Loads grid xy only for 1st frame:
        if n == Start_frame
            XY_Grid = load([Backup_folder filesep [Program  '_' PIV_Grid  '_' Animal '_Backup_'] 'XY_Grid.mat']);
            x = XY_Grid.x;
            y = XY_Grid.y;
        end
        
        %%% EITHER loads u,v backups OR uq,vq (contain NaNs)(edit 2.2):
        PIV_backup_n = load([Backup_folder filesep [Program  '_' PIV_Grid  '_' Animal '_Backup_'] num2str(n,digits_format) '.mat']);
        if exist('PIV_type','var') && strcmp(PIV_type,'q')
            u = PIV_backup_n.uq;
            v = PIV_backup_n.vq;
            tag = 'q';
        else
            u = PIV_backup_n.u;
            v = PIV_backup_n.v;
            tag = '';
        end
        
        %%% Vector field:
        if Velocity_display == 1
            figure('PaperPositionMode','auto')     % REQUIRED SO THAT "SAVEAS" DOES THE SAME AS "FILE->SAVE AS" FROM THE FIGURE WINDOW. ALSO REQUIRED TO SAVE IMAGES WITHOUT BORDERS
            %All raw images assumed to be in same folder (1.3):
            imshow([Path_folder_Raw filesep Filename_Raw num2str(n,digits_format)  image_extension_input_Raw],'Border', 'tight') % REQUIRED TO REMOVE BORDERS IN THE FIGURE (.fig) DIPLAYED BY MATLAB
            hold all 
            quiver(x,y,u,v,1.,arrow_color);                                % put back to 1 (2.0GMa)
            %quiver(x,y,u,v,0,arrow_color);                                % 0 for no autoscale, y for yellow
            % title([ 'Velocity Field      '  'Grid size: ' PIV_Grid ' (' num2str(size_grid_2) ' pixels)      ' Program,'   ', Version_name, '     ', When,'   ','Frame: ',  num2str(n,digits_format) '-' num2str(n+1,digits_format)],'FontSize',7)
            disp(['Saving Velocity map for interframe #' num2str(n,digits_format) '-' num2str(n+1,digits_format)]);
            print (print_format,print_resolution_PIV,[Frame_folder filesep  Program  '_' PIV_Grid  '_' Animal '_' tag 'Velocity_'  num2str(n,digits_format)  image_extension_output_PIV]); % added tag (2.2)
            hold off
            close
        end
        
        %%% Rebuilds x_image, y_image, u_image, v_image, ONLY if u,v and NOT uq,vq were loaded (2.2):
        if Speed_display == 1 || Divergence_display == 1
            if isempty(tag)
                [x_image, y_image, u_image, v_image] = PIV_Interpolate(image_size,x,y,u,v, 'cubic', PIV_Grid);
            else
                disp('WARNING: speed map and velocity divergence map cannot be computed with (uq,vq) that contain NaNs!! They will be skipped.')
            end
        end
        
        %%% Speed display 
        if Speed_display == 1 && isempty(tag)
            
            % moved inside if 2.0GMa:
            u_image_square = u_image.^2;
            v_image_square = v_image.^2;
            speed_field = sqrt(u_image_square + v_image_square);
            
            figure('PaperPositionMode','auto')     % REQUIRED SO THAT "SAVEAS" DOES THE SAME AS "FILE->SAVE AS" FROM THE FIGURE WINDOW. ALSO REQUIRED TO SAVE IMAGES WITHOUT BORDERS
            imshow(speed_field,'Border', 'tight')  % REQUIRED TO REMOVE BORDERS IN THE FIGURE (.fig) DIPLAYED BY MATLAB
            colormap(jet);
            caxis([0 max(max(speed_field))]); % very important so as not to saturate the colorscale to value 1
            %caxis([0 max(max(norm_speed_field))]);
            %colorbar
            % title([ 'Speed Field      ' '    Grid size: ' PIV_Grid  ' (' num2str(size_grid_2) ' pixels)      '  Program,' ', Version_name, '   ', When,'   ','Frame: ', num2str(n,digits_format) '-' num2str(n+1,digits_format)],'FontSize',7)
            disp(['Saving Speed map for interframe #' num2str(n,digits_format) '-' num2str(n+1,digits_format)]);
            print (print_format,print_resolution_PIV,[Frame_folder filesep  Program  '_' PIV_Grid  '_' Animal '_' tag 'Speed_'  num2str(n,digits_format)  image_extension_output_PIV]); % added tag (2.2)
            close
        end
        

        %%% Velocity gradient display:
        if Divergence_display == 1 && isempty(tag)
            
            % moved inside if 2.0GMa:
            [ux_image, uy_image] = gradient(u_image,1);
            [vx_image, vy_image] = gradient(v_image,1);
            div_velocity = ux_image + vy_image;
            
            figure('PaperPositionMode','auto')     % REQUIRED SO THAT "SAVEAS" DOES THE SAME AS "FILE->SAVE AS" FROM THE FIGURE WINDOW. ALSO REQUIRED TO SAVE IMAGES WITHOUT BORDERS
            imshow(div_velocity,'Border', 'tight') % REQUIRED TO REMOVE BORDERS IN THE FIGURE (.fig) DIPLAYED BY MATLAB
            colormap(jet);
            caxis([min(min(div_velocity)) max(max(div_velocity))]); % THE VELOCITY DIVERGENCE CAN BE NEGATIVE!!
            % title([ 'Velocity Divergence      ' '    Grid size: ' PIV_Grid  ' (' num2str(size_grid_2) ' pixels)      '  Program,' ', Version_name, '   ', When,'   ','Frame: ', num2str(n,digits_format) '-' num2str(n+1,digits_format)],'FontSize',7)
            disp(['Saving Divergence map for interframe #' num2str(n,digits_format) '-' num2str(n+1,digits_format)]);
            print (print_format,print_resolution_PIV,[Frame_folder filesep  Program  '_' PIV_Grid  '_' Animal '_' tag 'Divergence_'  num2str(n,digits_format)  image_extension_output_PIV]); % added tag (2.2)
            close
        end
        
        %%% Progress display (2.0GMa):
        progressbar(n / (Final_frame-1));
    end
    loop_two = toc;
    disp(['Loop #2 duration: ' num2str(loop_two) ' seconds'])
end

% Update of txt file:
Today_Full_End = datestr(now);
Info = ['Program "' Program '" version ' Version_name ' with grid size ' PIV_Grid  '  Run on animal "' Animal '". Started on ' Today_Full_Start ', ended on ' Today_Full_End];
dlmwrite([Gridsize_folder filesep Today '_' Program '_' Version_name '_' PIV_Grid  '_' Animal '.txt'], Info, 'delimiter', '', 'newline','pc')




%% History %%

% 24/02/2014: 2.2
% - now saves 2 versions of velocity fields in PIV backups:
%   * (u,v) for segmentation/tracking use, that were interpolated and are "NaN free"
%   * (uq,vq) for quantitative analysis, that were NOT interpolated and still contain all the NaNs added by the filters, including "SignalFilterNan"
% - created "SignalFilterNan" from "SignalFilter" that puts NaNs instead of 0s where gray level = 0
% - removed global filter and changed parameters of filters according to discussions with Anaïs Bailles
% - applying SignalFilter on a box twice as small as before to limit overflowing beyond the animal
% - removed Border_parameter_PIV to always use "tight".
% - fixed minor bugs to display speed and divergence maps
% - added messages specifying which filter is being applied
% - changed names of figure saved to Velocity, Speed and Divergence.

% 2011-10-14: 2.1.1
% - absurd displacement values (beyond image, Inf...) & non-image data filters replaced by functions.
% - increased boxe size for "SignalFilter" to 2*size_grid_2

% 2011-10-13: 2.1
% - use of polymask.mat for PIV filtering removed.
% - filter of absurd displacement values (beyond image, Inf...) added.
% - non-image data filter added (regions of raw images with 0 gray level).

% 15/02/2011: 2.0GMc
% - replaced specification of grid size by number "grid_size", by string "PIV_Grid" (that replaced "grid_size_letter")

% 08/02/2011: 2.0GMb CHANGED NAME FROM "PIV" TO "Particle_Image_Velocimetry"
% - removed version from program name and changed name
% - added name program version, animal name and grid used in progressbar display
% - only one number of digits specified now: "digits" + "digit_format" now defined in AIA.

% 06/12/2010: 2.0GMa
% - minor adjustments to make it compatible with n_raw_images >1
% - use of function "progressbar" instead of "waitbar"
% - changed  "Image_..."  to  "image_..." 

% 01/09/2010: 2.0b
% - fixed bug: "Path_folder_PIV" to "Backup_folder" when loading backups to
% plot images

% 29/06/2010: 2.0a WORKS
% - removed opening/closing of matlabpool that is done in AIA
% - skips 1st loop if Figure_Replot_PIV == 1
% - new name (1.4 to 2.0) to acknowledge implementation of parallel
% computing.

% 28/06/2010:1.4c WORKS
% - execute display (2nd loop) ONLY if >=1 display has been selected
% - update the txt file with starting time AND ending time

% 21/06/2010: 1.4b WORKS
% - now only saves u,v (Backup) for each frame and x,y (Backup_XY_Grid) when processing "Start_frame".
% - use of homemade function "PIV_Interpolate" to get x,y_image and u,v_image
% - plot and save all graphics AFTER the main (parfor) loop.
% - progress display for parfor AND classic loop
% - Now uses "naninterp(mu,mv,'linear')" instead of "naninterp(mu,mv,'linear','polymask.mat',x,y)"
% which REMOVES EVERY SINGLE NaNs in the field and allows better subsequent
% interpolation using "interp2" within "PIV_Interpolate". Checked the
% relevance and the locality of modified values.

% 17/06/2010: 1.4a NOT WORKING (imshow issue)
% - now use a parfor loop
% - removed U_image and V_image that could require a lot of memory and were
% not used anymore.
% - remove N_Vfields, checking of parameter Multiple_Backups

% 18/05/2010:
% - Removed case where each raw image is in a different folder (Benoit)

% 06/04/10: 1.2
% - only saves variables called by CT to lighten backups

% 31/03/2009: start modifying Fanny's program

% 07/04/2009: switch to original_image_PIV 01
% - interpolation of u v to u_image, v_image where each image pixel has
% been assigned a velocity using 3 types of interpolation: (bi)linear,
% first neighbors and cubic (look up "interp2" in Matlab help)
% - graphic display of the speed in the whole image with colormap(jet)

% 08-09/04/2009: switch to 0.1a, 0.1b, 0.1c:
% - Selection of the cubic interpolation function.
% - using the velocity gradients to extrapolates u,v values to borders (yielding u_extended, v_extended)
% before performing the interpolation to get u_image, v_image.
% - graphic display of the velocity dirvergence and use to test relevant
% grid sive to use for PIV

% 14/04/2009: switch to 0.1d
% - creates subfolders for each grid size now
% - change of file naming to make quantity names appear first

% 27/04/2009
% - introduced "Image_extension_input_Raw" to specify the source image format

% 06/07/2009: switch to "PIV 1.0a"
% - possibility to choose which images will be saved
% - choose input and output number of digits
% - choose input/output image format
% - choose the resolution of output image
% - choose between saving individual backup files for each interframe OR an
% overall backup file (as in older version)

% 21/07/2009: switch to "PIV 1.0c" occured earlier
% - Doesn't specify the date when creating the folder because, in contrast
% to SIA, CT, DTaS, PIV doesn't depend on the segmentation.

% 21-22/07/2009: switch to 1.0d
% - Adjustments for execution from Advanced_Image_Inalysis and use for
% Fanny's segmentation set of programs.

% 23/07/2009: switch to 1.0e
% - Removed the date from the name of the folder created for PIV since it
% has to be executed only once.

% 18/08/2009: switch to 1.0f then 11
% - removed the program version from all filenames given + created a txt
% file containing all info on date and program version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


