function [ind_side_dilated, side_I_dilated, side_J_dilated] = SideDilatator(image_size, ind_side, value)

% Version 1.0
% Boris Guirao
%
% will dilate initial side made up by indices "ind_side" of "value" layers of pixels.
%
% INPUTS: image size (1x2 vector), list of indices (COLUMN VECTOR!) making up the side to
% dilate(ind_side), number of pixels to add around initial side
%
% OUTPUTS: list of indices making up the dilated side, list of subscripts
% such that ind_side_dilated = sub2ind(image_size, side_I_dilated, side_J_dilated);

%% Code

ind_side_dilated = ind_side;                                               % initialize

for dil = 1:value                                                          % dilatation of pixels made one layer after another
    n_pix = length(ind_side_dilated);
    [side_I side_J] = ind2sub(image_size,ind_side_dilated);
    side_IJ = [side_I side_J];
    
    % dilatation:
    dilate_one_pix = ones(n_pix,1);
    side_I_plus = side_I + dilate_one_pix;
    side_I_minus = side_I - dilate_one_pix;
    side_J_plus = side_J + dilate_one_pix;
    side_J_minus = side_J - dilate_one_pix;
    
    % new pixels (starting right, finishing upper right, rotating clockwise = 8 position)
    one = [side_I side_J_plus];
    two = [side_I_plus side_J_plus];
    three = [side_I_plus side_J];
    four = [side_I_plus side_J_minus];
    five = [side_I side_J_minus];
    six = [side_I_minus side_J_minus];
    seven = [side_I_minus side_J];
    eight = [side_I_minus side_J_plus];
    
    % merging all pixels together making up the dilated side:
    side_IJ_dilated = unique([side_IJ ; one ; two ; three ; four ; five ; six ; seven ; eight], 'rows');
    side_I_dilated = side_IJ_dilated(:,1);
    side_J_dilated = side_IJ_dilated(:,2);
    
    % removes pixels out of image:
    row_I_out = find(side_I_dilated > image_size(1) | side_I_dilated <= 0);
    row_J_out = find(side_J_dilated > image_size(2) | side_J_dilated <= 0);
    row_out = unique([row_I_out ; row_J_out]);
    all_rows = (1:length(side_I_dilated))';
    row_keep = setdiff(all_rows, row_out);
    
    % crops side_I/J_dilated to indices kept:
    side_I_dilated = side_I_dilated(row_keep);
    side_J_dilated = side_J_dilated(row_keep);
    
    % gets linear indices:
    ind_side_dilated = sub2ind(image_size, side_I_dilated, side_J_dilated);
end

%% History

% 09-10/06/2010: creation