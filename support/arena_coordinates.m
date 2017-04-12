function [x, y, z, p_rad, rot180] = arena_coordinates
% FUNCTION [x, y, z, p_rad, rot180] = arena_coordinates
%
% Calculates the cartesian coordinates of every pixel in a cylindrical LED
% arena. Also calculates the distance between pixels in radians (where the
% arena's circumference is 2*pi.) Check the user-defined arena parameters 
% within this script to see if it fits your specific LED arena.
% 
% This script makes assumptions about how an arena is typically mounted 
% (e.g. what part of the arena is considered "straight ahead" of the fly.)
% For differently rotated arenas, you can use the rotate_coordinates script 
% at the end of this one to adjust the arena's coordinates as desired.
%
% outputs:
% x/y/z: matrices of cartesian coordinates for all pixels in the arena
% p_rad: distance between pixels (along rows/column directions) in radians
% rot180: 1 if physical arena is flipped upside-down (0 otherwise)
 

%% user-defined arena parameters
Psize = 16; %# of pixels per row (or column) of a single LED panel
Pcols = 12; %# of columns of LED panels 
Prows = 3; %# of rows of LED panels
Pcircle = 18; %# of panel columns that would fully enclose the arena
rot180 = 1; %if physical arena is flipped upside-down
model = 'poly'; %'smooth' or 'poly' cylinder 
%%

rows = Prows*Psize; %# of rows of pixels in arena
cols = Pcols*Psize; %# of columns of pixels in arena
Pan_rad = 2*pi/Pcircle; %radians between panel columns
p_rad = Pan_rad/Psize; %radians between pixels

%calculate height (z) of each pixel
z = p_rad*(1-rows)/2:p_rad:p_rad*(rows-1)/2; 
z = repmat(flipud(z'),1,cols);

%calculate the angle of each panel column's center from straight ahead
cphi = -Pan_rad*(Pcols-1)/2:Pan_rad:Pan_rad*(Pcols-1)/2;

%for a fully-enclosed arena with an even # of col panels, assume the panel
%center is straight ahead (rather than a vertex)
if mod(Pcols,2) && Pcols==Pcircle
    cphi = cphi - Pan_rad/2; %shift panel center to straight ahead
end
cphi = repmat(cphi,[Psize 1]);
cphi = reshape(cphi,[1 cols]);
        
%calculate each pixel's distance relative to its panel center
points = (p_rad-Pan_rad)/2:p_rad:(Pan_rad-p_rad)/2;
points = repmat(points,[1 Pcols]);
        
switch model
    case 'smooth' %model as smooth-surfaced cylinder (radius = 1)
        %calculate angle from straight ahead of each pixel column
        col_phi = cphi + points;
        
        %calculate x,y coordinates from angle
        x = repmat(sin(col_phi), rows, 1);
        y = repmat(cos(col_phi), rows, 1);
        
    case 'poly' %model as polygonal cylinder
        apothem = Pan_rad/(2*tan(pi/Pcircle));
        
        %calculate x,y coordinates of each pixel's panel column center
        x = apothem*sin(cphi);
        y = apothem*cos(cphi);
        
        %adjust each pixel's coordinate from center of panel
        x = x + points.*cos(-cphi);
        y = y + points.*sin(-cphi);
        x = repmat(x,[rows 1]);
        y = repmat(y,[rows 1]);
end

%make any desired adjustments to arena position/orientation here

end