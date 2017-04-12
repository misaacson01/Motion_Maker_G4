function [Pats, arena_phi, arena_theta, p_rad, true_step_size, rot180] = Motion_Maker_G4(param)
% FUNCTION [Pats, arena_phi, arena_theta, true_step_size, rot180] = Motion_Maker_G4(param)
%
% inputs:
% see Motion_Maker_G4_script for list of inputs needed
%
% outputs:
% Pats: array of brightness values for each pixel in arena (multiple frames
%      in 3rd dimension)
% arena_phi/theta: spherical coordinates of pixels in arena
% p_rad: distance between pixels (along rows/column directions) in radians
% true_step_size: value (in radians) of corrected step_size (may have been
%      altered in order to divide evenly into spatial frequency of pattern)
% rot180: if physical arena is flipped upside-down


%% check for parameter errors
if param.gs_val~=1 && param.gs_val~=4
    error('gs_val must be either 1 or 4');
end
if any(mod(param.levels,1)~=0) || any(param.levels<0) || any(param.levels>(2^param.gs_val)-1)
    error(['for selected gs_val, all levels must be positive integers ' ...
        'no greater than ' num2str((2^param.gs_val)-1)]);
end
if any(strcmpi(param.motion_type(1), {'r' 't' 'e'}))==0
    error('invalid choice of motion type')
end
if any(strcmpi(param.pattern_type(1:2), {'sq' 'si' 'ed' 'st' 'of'}))==0
    error('invalid choice of pattern type')
end
if any(strcmpi(param.pattern_fov(1), {'l' 'f'}))==0
    error('invalid choice of pattern fov')
end
if any(strcmpi(param.dot_size(1), {'s' 'd'}))==0
    error('invalid choice of dot size')
end
if any(strcmpi(param.dot_occ(1), {'c' 's' 'm'}))==0
    error('invalid choice for dot occlusion')
end
if param.aa_samples<1 || mod(param.aa_samples,1)~=0
    error('aa_samples must be a positive integer')
end
if any([param.arena_pitch, param.pole_coord, param.motion_angle, param.spat_freq, ...
        param.step_size, param.dot_radius, param.sa_mask(1:3), param.long_lat_mask(1:4)]>2*pi)
    error('some units might be in degrees rather than radians (error catches units > 2*pi)');
end
if param.duty_cycle<0 || param.duty_cycle > 100
    error('duty_cycle must be value between 0-100');
end

%% calculate arena coordinates
%get starting arena parameters
[arena_x, arena_y, arena_z, p_rad, rot180] = arena_coordinates;
param.p_rad = p_rad;
[param.rows, param.cols] = size(arena_x);

%shift arena coordinates for arena pitch
[arena_x, arena_y, arena_z] = rotate_coordinates(arena_x, arena_y, arena_z, [0 param.arena_pitch 0]);
[arena_phi, arena_theta, ~] = cart2sphere(arena_x, arena_y, arena_z);

%% create pattern
switch lower(param.pattern_type(1:2))
    case 'st' %starfield
        [Pats, num_frames, true_step_size] = make_starfield(param, arena_x, arena_y, arena_z);
    case 'of' %off_on
        true_step_size = param.step_size;
        [Pats, num_frames] = make_off_on(param);
    otherwise %grating or edge
        [Pats, num_frames, true_step_size] = make_grating_edge(param, arena_x, arena_y, arena_z);
end

%% apply masks to pattern
%solid angle mask
if param.sa_mask(3)<pi
    mask = sa_mask(arena_x, arena_y, arena_z, param.sa_mask, param.aa_samples);
    mask = repmat(mask,[1, 1, num_frames]);
    Pats = Pats.*mask + param.levels(3)*ones(size(Pats)).*(1-mask);
end
%lattitude-longitude mask
if abs(diff(param.long_lat_mask(1:2)))<2*pi || abs(diff(param.long_lat_mask(3:4)))<pi
    mask = long_lat_mask(arena_phi, arena_theta, param.long_lat_mask, param.aa_samples);
    mask = repmat(mask,[1, 1, num_frames]);
    Pats = Pats.*mask + param.levels(3)*ones(size(Pats)).*(1-mask);
end
    
%round pattern values
Pats = round(Pats);

%% Make final adjustments to pattern
%left-right flip the right half of the pattern (if desired)
if param.flip_right
    Pats(:,1+param.cols/2:end,:) = fliplr(Pats(:,1+param.cols/2:end,:));
end

%flip pattern according to physical arena placement
if rot180
    Pats = rot90(Pats,2);
end

%shift order of pattern frames (if desired)
if param.phase_shift~=0
    Pats = circshift(Pats,-param.phase_shift,3);
end

%add background frame as frame 1 (if desired)
if param.back_frame
    num_frames = num_frames+1;
    Pats(:,:,2:num_frames) = Pats;
    Pats(:,:,1) = param.levels(3);
end

end