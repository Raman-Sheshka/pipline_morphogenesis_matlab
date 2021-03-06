function ALL_FRAMES = BuildTrackingStructure(trackingPath, firstTrackingFrame, lastTrackingFrame)
% Export into function of the ALL_FRAMES building code from the Dynamic
% Processing Stage of the segmentation pipeline
% by Stephane Rigaud
%
% Building ALL_FRAMES,FRAME_n and cell_tags_n (0.8c) %%
% NB: See Boris' notebook #4 to see cell tags and how lists were built up

% Loads last frame to get max number of divisions and built AN Vector style: 

correspondence_last = dlmread([trackingPath filesep 'correspondence_' num2str(lastTrackingFrame) '.txt']);

width_last = size(correspondence_last, 2);

% Initialization of frist frame:
NEXT_correspondence = dlmread([trackingPath filesep 'correspondence_' num2str(firstTrackingFrame) '.txt']);
size_NEXT_correspondance = size(NEXT_correspondence);
% pad if necessary
if size_NEXT_correspondance(2) < width_last
    NEXT_correspondence = [NEXT_correspondence zeros(size_NEXT_correspondance(1),width_last-size_NEXT_correspondance(2))];
elseif size_NEXT_correspondance(2) > width_last % update width_last
    width_last = size_NEXT_correspondance(2);
end

% "unique" to get all cell numbers only once
NEXT_all_cells_RN = unique(NEXT_correspondence(:,1));                  
NEXT_n_cells = length(NEXT_all_cells_RN);
NEXT_cell_tags_column = zeros(NEXT_n_cells, 1);

progressbar('Building tracking color patches');
for n = firstTrackingFrame:lastTrackingFrame
    % CURRENT frame quantities:
    % Update with already loaded next frame quantities:
    correspondence = NEXT_correspondence;
    all_cells_RN = NEXT_all_cells_RN;
    cell_tags_column = NEXT_cell_tags_column; % already partially filled
    
    % Reads txt files generated by "tracking.exe":
    % load lists of Border cells ("border_cell_RN_1.txt, _2...) (0.8b)
    Border_cells = dlmread([trackingPath filesep 'border_cells_RN_' num2str(n) '.txt']);
    % load coalesced cells list
    coalesced_cells_RN = dlmread([trackingPath filesep 'coalesced_cells_RN_' num2str(n) '.txt']);
    % load just divided cells list
    just_divided_cells_RN = dlmread([trackingPath filesep 'just_divided_cells_RN_' num2str(n) '.txt']);
    % load new cells list 
    new_cells_RN = dlmread([trackingPath filesep 'new_cells_RN_' num2str(n) '.txt']);
    
    % replacing -1 by empty sets:
    if coalesced_cells_RN(1) == -1; 
        coalesced_cells_RN = []; 
    end                                               
    if just_divided_cells_RN(1) == -1
        just_divided_cells_RN= []; 
    else
        just_divided_cells_RN = reshape(just_divided_cells_RN,numel(just_divided_cells_RN),1); %2011-11-30: reshape because of new format
    end                                    
    if new_cells_RN(1) == -1; 
        new_cells_RN= []; 
    end
    
    % restriction for list not existing in last frame:
    if n < lastTrackingFrame
        dividing_cells_RN = dlmread([trackingPath filesep 'dividing_cells_RN_' num2str(n) '.txt']);
        if dividing_cells_RN(1) == -1; 
            dividing_cells_RN = []; 
        end
    else
        dividing_cells_RN = [];
    end
    
    % Buildging  "no_pb_just_divided_cells_RN":
    no_pb_just_divided_cells_RN = setdiff(just_divided_cells_RN, coalesced_cells_RN); % removes coalesced cells from the list of divided cells
    
    % Building AN lists:
    coalesced_cells_rows = ismember(correspondence(:,1), coalesced_cells_RN);
    coalesced_cells_AN = correspondence(coalesced_cells_rows, 2:end);

    % NEXT frame quantities:
    if n < lastTrackingFrame
    
        % Re-formatting correspondence to largest possible AN:
        % 2012-09-06
        % load correspondence file
        NEXT_correspondence = dlmread([trackingPath filesep 'correspondence_' num2str(n+1) '.txt']);
        size_NEXT_correspondance = size(NEXT_correspondence);
        % pad if necessary
        if size_NEXT_correspondance(2) < width_last
            NEXT_correspondence = [NEXT_correspondence zeros(size_NEXT_correspondance(1),width_last-size_NEXT_correspondance(2))]; %#ok<AGROW>
        elseif size_NEXT_correspondance(2) > width_last
            width_last = size_NEXT_correspondance(2); % update width_last
            % resize coalesced_cells_AN
            coalesced_cells_AN = [coalesced_cells_AN zeros(size(coalesced_cells_AN,1),1)]; %#ok<AGROW>
        end               
        
        % "unique" to get all cell numbers only once
        NEXT_all_cells_RN = unique(NEXT_correspondence(:,1));
        NEXT_n_cells = length(NEXT_all_cells_RN);
        % intialization of the NEXT_cell_tags
        NEXT_cell_tags_column = zeros(NEXT_n_cells,1);
        
        % reading txt files generated by "tracking":
        NEXT_coalesced_cells_RN = dlmread([trackingPath filesep 'coalesced_cells_RN_' num2str(n+1) '.txt']);
        NEXT_just_divided_cells_RN = dlmread([trackingPath filesep 'just_divided_cells_RN_' num2str(n+1) '.txt']);
        if NEXT_coalesced_cells_RN(1) == -1; 
            NEXT_coalesced_cells_RN = []; 
        end 
        if NEXT_just_divided_cells_RN(1) == -1; 
            NEXT_just_divided_cells_RN = []; 
        end
        
        % Building AN lists:
        NEXT_coalesced_cells_rows = ismember(NEXT_correspondence(:,1), NEXT_coalesced_cells_RN);
        NEXT_coalesced_cells_AN = NEXT_correspondence(NEXT_coalesced_cells_rows,2:end);
        NEXT_just_divided_cells_rows = ismember(NEXT_correspondence(:,1), NEXT_just_divided_cells_RN);
        NEXT_just_divided_cells_AN = NEXT_correspondence(NEXT_just_divided_cells_rows,2:end);        
        
        % Building "just_coalesced_cells_AN", "old_coalesced_cells_AN", "coalescing_cell_AN":
        NEXT_just_coalesced_cells_AN = setdiff(NEXT_coalesced_cells_AN, coalesced_cells_AN, 'rows');     % NEXT_ because refers to next frame
        NEXT_old_coalesced_cells_AN = intersect(NEXT_coalesced_cells_AN, coalesced_cells_AN, 'rows');    % NEXT_ because refers to next frame
        coalescing_cells_AN = setdiff(NEXT_just_coalesced_cells_AN, NEXT_just_divided_cells_AN, 'rows'); % REMOVES AN OF DAUGHTERS FROM WRONG DIVISIONS
                                                                                                         % NB: not NEXT_ because AN numbers are the same in both frames since
                                                                                                         % dividing cells are not involved and can thus refer to current frame
        % Building "just_coalesced_cells_RN", "old_coalesced_cells_RN", "coalescing_cell_RN":                                                                                                 
        NEXT_just_coalesced_cells_rows = ismember(NEXT_correspondence(:,2:end), NEXT_just_coalesced_cells_AN, 'rows');
        NEXT_just_coalesced_cells_RN = unique(NEXT_correspondence(NEXT_just_coalesced_cells_rows,1));
        
        NEXT_old_coalesced_cells_rows = ismember(NEXT_correspondence(:,2:end), NEXT_old_coalesced_cells_AN, 'rows');
        NEXT_old_coalesced_cells_RN = unique(NEXT_correspondence(NEXT_old_coalesced_cells_rows,1));
        
        % NB: SOME RN CAN BELONG TO BOTH "NEXT_just_coalesced_cells_RN" and
        % "NEXT_old_coalesced_cells_RN" when a new cell joins a region already coalesced
      
        % Only in CURRENT FRAME
        coalescing_cells_rows = ismember(correspondence(:,2:end), coalescing_cells_AN(:,1:size(correspondence,2)-1), 'rows');
        coalescing_cells_RN = unique(correspondence(coalescing_cells_rows,1));   
        
        % Filling cell_tags_column AND NEXT_cell_tags_column:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NB: ONLY FILLED WHEN THERE IS A NEXT FRAME!!
        cell_tags_column(coalescing_cells_RN) = -2;                         
        NEXT_cell_tags_column(NEXT_old_coalesced_cells_RN) = 3; 
        NEXT_cell_tags_column(NEXT_just_coalesced_cells_RN) = 2; % NB: just coalesced cells have priority since they are tagged at the end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end % end if not last frame
    
    % Complete cell_tags_column:
free_indices = find(cell_tags_column == 0); % indices not previously taken by higher priority indices

% Filter indices to avoid overwriting:
filtered_dividing_cells_RN = intersect(dividing_cells_RN, free_indices); % NB: dividing_cells empty if n = lastframe
filtered_no_pb_just_divided_cells_RN = intersect(no_pb_just_divided_cells_RN, free_indices);
filtered_new_cells_RN = intersect(new_cells_RN, free_indices);

% fills up cell_tags_column:
cell_tags_column(filtered_dividing_cells_RN) = -1;
cell_tags_column(filtered_no_pb_just_divided_cells_RN) = 1;
cell_tags_column(filtered_new_cells_RN) = 4;

%% Full cell_tags:
cell_tags = [all_cells_RN cell_tags_column];

%%% saving cell_tags in a structure:
THIS_FRAME.Border_cells = Border_cells;
THIS_FRAME.cell_tags = cell_tags;
this_name = ['FRAME_' num2str(n)];
ALL_FRAMES.(this_name) = THIS_FRAME;
    
progressbar(n/lastTrackingFrame);
end % end loop over frames

end

%% Historique

% 19/01/2017 - v1