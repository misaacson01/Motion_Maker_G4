function arena_projection(Pats, plot_type, arena_phi, arena_theta, p_rad, frame, param)
% FUNCTION arena_projection(Pats, gs_val, x, y, axes, dot_size, frame, rot180)
% 
% Plots a single frame of the input pattern as a scatter plot of dots.
%
% inputs:
% Pats: 2-D or 3-D array of brightness values for each pixel of pattern
% gs_val: # of brightness bits (1 or 4)/whether Pats values are 0-1 or 0-15
% x/y: matrices of x/y-coordinates corresponding to each pixel in the arena
% axes: [min-x max-x min-y max-y] sets the axes limits for the figure
% dot size: sets the size of dot of each pixel
% frame: for a 3-D Pats variable, frame sets the 3rd dimension index to
%        specify the 2-D pattern to be displayed
% rot180: 1 if arena is mounted upside-down (0 otherwise)

%figure size settings
max_w = 1.15;
max_h = 0.5;
center_h = 0.33;
center_w = 0.5;

%determine drawing parameters based on plot type
switch plot_type
    case 1 %mercator projection
        x = rad2deg(arena_phi);
        y = rad2deg(arena_theta)-90; %convert to lattitude
        dot_size = 140*p_rad;
        axes = [-180 180 -90 90];
    case 2 %grid projection
        rows = length(arena_phi(:,1));
        cols = length(arena_phi(1,:));
        x = repmat(1:cols,rows,1);
        y = repmat((1:rows)',1,cols);
        dot_size = 320*p_rad;
        axes = [0 cols+1 0 rows+1];
end

%turn values into vectors for easy plotting
num_pixels = numel(x);
x_vec = reshape(x,[num_pixels 1]);
y_vec = reshape(y,[num_pixels 1]);

%convert Pats to vector of values between 0 and 1
if param.rot180==1
    Pats = rot90(Pats,2);
end
Pats_vec = reshape(Pats(:,:,frame)/(2^param.gs_val - 1),[num_pixels 1]);

%draw pixels in green
Color = [zeros(num_pixels,1) Pats_vec zeros(num_pixels,1)];

%plot pattern
scatter(x_vec,y_vec,dot_size,Color,'filled')
grid on
axis(axes)

%scale plot size to match pattern size
ax = gca;
size_ratio = diff(axes(3:4))/diff(axes(1:2));
width = max_w;
height = width*size_ratio;
if height>max_h
    width = width*(max_h/height);
    height = max_h;
end
ax.OuterPosition = [center_w-width/2 center_h-height/2 width height];

end