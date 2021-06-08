function [x_image, y_image, u_image, v_image] = InterpolatePIV(image_size,x,y,u,v, interp_type, PIV_Grid)
%
% [x_image, y_image, u_image, v_image] = InterpolatePIV(image_size,x,y,u,v, interp_type, PIV_Grid)
%
%
% INPUTS:
% - x,y,u,v = matpiv outputs (x,y gives location coordinates where displacement field u,v was caculated)
% - interp_type = type of interpolation. 'cubic' recommended (look up "interp2" function to see accepted types)
% - image_size = size of image (1x2 matrix)
% - PIV_Grid = 'XS'(1),'S'(2),'M'(4),'L'(8),'XL'(16),'XXL'(32)
%
% OUTPUTS:
% - x_image, y_image, u_image, v_image = x,y,u,v values interpolated (between points) and extrapolated (up to borders) to all image pixels
%
% Version 1.2
% Boris Guirao


%% Conversion of "PIV_Grid" into "grid_size" (1.1) %%

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
    [x_image, y_image, u_image, v_image] = deal([]); % 1.2
    return
end


%% Creates x_image, y_image %%

%%% defining "size_grid_1" "size_grid_2":
if grid_size<=16
    size_grid_2 = 8*grid_size;
else                                                                       % for XXL grid_size
    size_grid_2 = 8*grid_size;
end

%%% get x,y bounds:
x_image_max = image_size(2);                                               % THE FIRST NUMBER CORRESPONDS TO THE NUMBER OF LINES AND THUS Y!!
y_image_max = image_size(1);

%%% Create meshgrid for interpolation:
[x_image, y_image] = meshgrid(1:x_image_max,1:y_image_max);                % Creates a "x_image" matrix of 1s in column 1, 2s in column
                                                                           % 2...x_image_maxs in the last column WITH y_image_max NUMBER OF LINES AND
                                                                           % a "y_image" matrix with 1s in line 1,...y_image_maxs in last line WITH
                                                                           % x_image_max NUMBER OF COLUMNS.
%%% number of lines and columns of the native PIV matrices x,y,u,v:
n_lines_piv = size(x,1);                                                   
n_columns_piv = size(x,2);

%%  Extrapolation of x,y,u,v on image borders: x_extended,y_extended,u_extended,v_extended %%

% NB: values on borders are extrapolated from last known values using u and v gradients in x and y directions

%%% x,y:
x_max = max(max(x));
y_max = max(max(y));
[x_extended_temp y_extended_temp] = meshgrid(0:size_grid_2/2:x_image_max, 0:size_grid_2/2:y_image_max);
n_lines = size(x_extended_temp,1);
n_columns = size(x_extended_temp,2);
x_extended_temp(:,1) = ones(n_lines,1);
x_extended_temp(:,n_columns) = x_image_max*ones(n_lines,1);
x_extended = x_extended_temp;
y_extended_temp(1,:) = ones(1,n_columns);
y_extended_temp(n_lines,:) = y_image_max*ones(1,n_columns);
y_extended = y_extended_temp;

%%% u:
[ux,uy] = gradient(u,size_grid_2/2);                % calculates u gradients along x and y directions. The distance between 2 points (2 x values for example) is size_grid_2/2

u_extended = zeros(n_lines,n_columns);              % u, v have same size as x,y, so do u,v_extended
u_extended(2:n_lines-1,2:n_columns-1) = u;          % puts u right in the middle of u_extended, leaving one line and columns of 0s on the borders

% First: extending the u matrix in the x direction using (filling u_extended
% first and last columns) ux values at its first and last columns:
u_extended(:,1) = [0;u(:,1);0]-[0;ux(:,1);0].*(size_grid_2/2-1);                                    % fills up the first column of u_extended using ux value in the first column
                                                                                                    % of ux multiplied by the difference of number of pixels
u_extended(:,n_columns) = [0;u(:,n_columns_piv);0]+[0;ux(:,n_columns_piv);0].*(x_image_max-x_max);  % fills up the last column of u_extended

% Second: extending the u_EXTENDED matrix in y direction using uy on the first and last lines:
% NB: now that the first and last columns of u_extended are filled up, uy is also used to extrapolate these values. For this purpose, the first and last values
% of uy on each lines (uy(1,1), uy(1,n_columns_piv) and uy(n_lines_piv,1), uy(n_lines_piv,n_columns_piv)) are used.
u_extended(1,:) = u_extended(2,:)-[uy(1,1) uy(1,:) uy(1,n_columns_piv)].*(size_grid_2/2-1);                                               % fills up the first line of u_extended
u_extended(n_lines,:) = u_extended(n_lines-1,:)+[uy(n_lines_piv,1) uy(n_lines_piv,:) uy(n_lines_piv,n_columns_piv)].*(y_image_max-y_max); % last line


%%% v:                 
[vx,vy] = gradient(v,size_grid_2/2);

v_extended = zeros(n_lines,n_columns);              % u, v have same size as x,y, so do u,v_extended
v_extended(2:n_lines-1,2:n_columns-1) = v;          % puts u right in the middle of u_extended, leaving one line and columns of 0s on the borders

% First extending the v matrix in the x direction using vx:
v_extended(:,1) = [0;v(:,1);0]-[0;vx(:,1);0].*(size_grid_2/2-1);                                    % fills up the first column of v_extended using vx value in the first column of vx multiplied by the difference of number of pixels
v_extended(:,n_columns) = [0;v(:,n_columns_piv);0]+[0;vx(:,n_columns_piv);0].*(x_image_max-x_max);  % fills up the last column

% Second extending the v matrix in y direction using vy:
v_extended(1,:) = v_extended(2,:)-[vy(1,1) vy(1,:) vy(1,n_columns_piv)].*(size_grid_2/2-1);
v_extended(n_lines,:) = v_extended(n_lines-1,:)+[vy(n_lines_piv,1) vy(n_lines_piv,:) vy(n_lines_piv,n_columns_piv)].*(y_image_max-y_max);


%% Interpolation of u,v_extended to each image pixel: u_image, v_image %%

u_image = interp2(x_extended,y_extended,u_extended,x_image,y_image, interp_type);
v_image = interp2(x_extended,y_extended,v_extended,x_image,y_image, interp_type);


%% History

% 26/02/2014: 1.2
% - added [x_image, y_image, u_image, v_image] = deal([]) in case PIV_Grid doesn't fit.

% 17/02/2011: 1.1
% - Replaced "grid_size" (number) by "PIV_Grid" (string)

% 21/06/2010: creation










