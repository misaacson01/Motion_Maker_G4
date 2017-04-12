function [Pats, num_frames, true_step_size] = make_starfield(param, arena_x, arena_y, arena_z)
% function [Pats, num_frames, true_step_size] = make_starfield(param, arena_x, arena_y, arena_z)
% 
% Creates moving pattern of starfield of dots
% 
% inputs: (all angles/distances/sizes in units of radians)
% arena_x/y/z: cartesian coordinates of pixels in arena
% param.motion_type: 'rotation', 'translation', or 'expansion-contraction'
% param.pattern_fov: 'local' or 'full-field'
% param.gs_val: bits of intensity value (1 or 4)
% param.levels: brightness level of [dots, background, mask]
% param.dot_level: 0 = dot brightness set to 1st level; 1 and 2 = random brightness (0-1st; 0 or 1st)
% param.arena_pitch: pitch of arena 
% param.step_size: distance dot travels per frame
% param.pole_coord: location of pattern pole [longitude, lattitude]
% param.sa_mask: [mask-center-longitude, mask-center-lattitude]
% param.motion_angle: direction of motion (0=left->right, positive values rotate clockwise)
% param.num_dots: number of random dots generated in arena space
% param.dot_radius: radius of each dot
% param.dot_size: 'static' or 'distance-relative' (where closer dots appear larger)
% param.view_radius: distance dots can be seen
% param.dot_occ: how occluding dots are handled ('closest', 'sum', or 'mean')
% param.aa_samples: # of samples taken to calculate a single pixel's brightness
% 
% outputs:
% Pats: array of brightness values for each pixel in the arena
% num_frames: # of frames in the Pats variable (3rd dimension in Pats)
% true_step_size: recalculated value of step size (in order to divide evenly into motion range)


%% calculate arena coordinates needed to produce for desired pattern
if strncmpi(param.pattern_fov,'l',1)==1
    if strncmpi(param.motion_type,'r',1)
        rotations = [-param.sa_mask(1:2) -param.motion_angle];
    else
        rotations = [-param.sa_mask(1:2) -param.motion_angle-pi/2];
    end
else
    rotations = [-param.pole_coord(1) -param.pole_coord(2)-pi/2 0];
end
[pat_x, pat_y, pat_z] = rotate_coordinates(arena_x, arena_y, arena_z, rotations);

% map arena coordinates to surface of unit sphere
[pat_phi, pat_theta, ~] = cart2sphere(pat_x, pat_y, pat_z);
[pat_x, pat_y, pat_z] = sphere2cart(pat_phi, pat_theta);


%% calculate random starting coordinates of dots
num_dots = param.num_dots;
view_rad = param.view_radius;
dots_x = (rand(num_dots, 1)*2 - 1)*view_rad;
dots_y = (rand(num_dots, 1)*2 - 1)*view_rad;
dots_z = (rand(num_dots, 1)*2 - 1)*view_rad;
[~, ~, dots_rho] = cart2sphere(dots_x, dots_y, dots_z);

% recalculate coordinates of dots that overlap the arena center   
if strncmpi(param.motion_type,'t',1)
    dots_pole_dist = hypot(dots_x, dots_y);
    idx = dots_pole_dist<=param.dot_radius;
    while any(idx)
        dots_x(idx) = (rand(sum(idx),1)*2 - 1)*view_rad;
        dots_y(idx) = (rand(sum(idx),1)*2 - 1)*view_rad;
        dots_pole_dist(idx) = hypot(dots_x(idx), dots_y(idx));
        idx = dots_pole_dist<=param.dot_radius;
    end
    [dots_phi, dots_theta, dots_rho] = cart2sphere(dots_x, dots_y, dots_z);
    
else
    idx = dots_rho<param.dot_radius;
    while any(idx)
        dots_x(idx) = (rand(sum(idx),1)*2 - 1)*view_rad;
        dots_y(idx) = (rand(sum(idx),1)*2 - 1)*view_rad;
        dots_z(idx) = (rand(sum(idx),1)*2 - 1)*view_rad;
        [~, ~, dots_rho(idx)] = cart2sphere(dots_x(idx), dots_y(idx), dots_z(idx));
        idx = dots_rho<param.dot_radius;
    end
    idx = dots_rho>=view_rad;
    dots_x(idx) = [];
    dots_y(idx) = [];
    dots_z(idx) = [];
    [dots_phi, dots_theta, dots_rho] = cart2sphere(dots_x, dots_y, dots_z);
    num_dots = length(dots_x);
end

%Method for maintaining even pixel density as dots travel to/from poles in 
%expansion-contraction motion:
%have some dots disappear before reaching poles and traverse equatorial
%region again (on a new randomised phi trajectory, but keeping same rho)
%before returning to their starting trajectory and completing a loop
if strncmpi(param.motion_type,'e',1)
    %array of 2 phi trajectories for each dot
    dots_2_phi = [dots_phi'; 2*pi*(rand(1,num_dots)*2 - 1)];
  
    %array of 1-D distance from pole (randomly chosen from x and y)
    xy_idx = logical(randi(2,[num_dots 1])-1);
    dots_pole_dist = nan(num_dots,1);
    dots_pole_dist(xy_idx) = dots_x(xy_idx);
    dots_pole_dist(~xy_idx) = dots_y(~xy_idx);
    
    %maximum allowed theta for each trajectory before the dot will switch
    max_theta = [pi-abs(asin(dots_pole_dist')); pi-(pi/2-abs(asin(dots_pole_dist')))];
    
    %index of trajectories - all dots start with original trajectories (row 1)
    traj_idx = logical([ones(1,num_dots); zeros(1,num_dots)]); 
    
    %randomize distances, otherwise poles will favor closer (i.e. larger) dots
    dots_rho = dots_rho(randperm(num_dots)); 
end


%% calculate brightness level of dots (difference from 2nd level)
if param.dot_level == 1  %for random dot brightness between 0 and 1st level
    dot_level = param.levels(1)*rand(num_dots, 1) - param.levels(2);
elseif param.dot_level == 2 %for random dot brightness (either 0 or 1st level)
    dot_level = param.levels(1)*(randi(2,[num_dots 1]) - 1) - param.levels(2);
else %for dot brightness = 1st level
    dot_level = param.levels(1)*ones(num_dots, 1) - param.levels(2);
end


%% calculate number of frames needed for complete loop of dot motion
if strncmpi(param.motion_type,'r',1)
    num_frames = round(2*pi/param.step_size);
    true_step_size = 2*pi/num_frames;
elseif strncmpi(param.motion_type,'t',1)
    num_frames = round(2*view_rad/param.step_size); 
    true_step_size = 2*view_rad/num_frames;
else
    num_frames = round(pi/param.step_size);
    true_step_size = pi/num_frames;
end


%% use standard (unpitched) arena coordinates to check if dots will fall within arena fov
[cpat_x, cpat_y, cpat_z] = rotate_coordinates(arena_x, arena_y, arena_z, [0 -param.arena_pitch 0]);
[cpat_phi, cpat_theta, ~] = cart2sphere(cpat_x, cpat_y, cpat_z);


%% pre-allocate space for pattern
Pats = zeros([size(arena_x) num_frames]);


%% loop for each frame
for f=1:num_frames
    

%% dot checking: calculate dot coordinates for standard arena coordinates
check_rotations = [rotations(1) rotations(2)+param.arena_pitch rotations(3)];
[cdots_x, cdots_y, cdots_z] = rotate_coordinates(dots_x, dots_y, dots_z, -check_rotations, 1);
[cdots_phi, cdots_theta, cdots_rho] = cart2sphere(cdots_x, cdots_y, cdots_z);


%% calculate apparent radius (solid angle) of each dot
if strncmpi(param.dot_size,'d',1) %closer dots appear larger
    dots_rad = atan(param.dot_radius./dots_rho); 
else %all dots appear the same size (as if all dots are on surface of unit sphere)
    dots_rad = param.dot_radius*ones(size(dots_rho));
end
    
    
%% find dots that will fall within arena coordinates
%valid dots fall within arena fov and within view radius
vdot_idx = cdots_theta<max(max(cpat_theta))+dots_rad & ...
                cdots_theta>min(min(cpat_theta))-dots_rad & ...
                cdots_phi<max(max(cpat_phi))+dots_rad & ...
                cdots_phi>min(min(cpat_phi))-dots_rad & ...
                cdots_rho<view_rad+dots_rad;
num_vdots = sum(vdot_idx);


%% map valid dots to surface of unit sphere (rho=1)
[sdots_x, sdots_y, sdots_z] = sphere2cart(dots_phi(vdot_idx), dots_theta(vdot_idx)); 


%% find pixels that are close to each valid dot
%convert pattern and dot matrices to arrays (valid dots in 3rd dim) and
%calculate distance between pixels/dots in each dimension individually
dx = abs(permute(repmat(sdots_x,[1 size(pat_x)]),[2 3 1]) - repmat(pat_x,[1 1 num_vdots]));
dy = abs(permute(repmat(sdots_y,[1 size(pat_x)]),[2 3 1]) - repmat(pat_y,[1 1 num_vdots]));
dz = abs(permute(repmat(sdots_z,[1 size(pat_x)]),[2 3 1]) - repmat(pat_z,[1 1 num_vdots]));

%valid pixels: close to a dot in all 3 dimensions individually
p_rad = param.p_rad;
dots_rad = permute(repmat(dots_rad(vdot_idx),[1 size(pat_x)]),[2 3 1]);
vpix_inds = find(dx<dots_rad+p_rad & dy<dots_rad+p_rad & dz<dots_rad+p_rad);


%% calculate exact distance from each valid dot center to all valid pixels
v_rads = acos((2-(dx(vpix_inds).^2 + dy(vpix_inds).^2 + dz(vpix_inds).^2))/2);


%% draw pixels
dot_level_rep = permute(repmat(dot_level(vdot_idx),[1 size(pat_x)]),[2 3 1]);

fullPats = nan(size(dx)); %3rd dim is for each valid dot
frac_array = nan(size(dx)); %fraction of brightness of corresponding elements in fullPats

if param.aa_samples==1 %no edge smoothing (frac_array will only have ones)
    face_idx = v_rads<=dots_rad(vpix_inds);
    frac_array(vpix_inds(face_idx)) = 1;
    fullPats(vpix_inds(face_idx)) = dot_level_rep(vpix_inds(face_idx));
else
    face_idx = v_rads<=dots_rad(vpix_inds)-p_rad/2; %index of pixels which display the dot's face
    edge_idx = v_rads>dots_rad(vpix_inds)-p_rad/2 & v_rads<dots_rad(vpix_inds)+p_rad/2; % " " dot's edge
    
    %for frac array: face=full brightness, edge=fraction (depending on amount dot overlaps with pixel fov)
    frac_array(vpix_inds(face_idx)) = 1; 
    frac_array(vpix_inds(edge_idx)) = 1-((v_rads(edge_idx)-(dots_rad(vpix_inds(edge_idx))-p_rad/2))/p_rad);
    fullPats(vpix_inds(face_idx|edge_idx)) = dot_level_rep(vpix_inds(face_idx|edge_idx));
    
    %multiply fullPats with frac_array to get actual brightness (edges scaled down)
    fullPats = fullPats.*frac_array; 
end


%% combine dots (depending on occlusion type)
if strncmpi(param.dot_occ,'s',1) %summing dots together (brightness clipped if exceeds gs_val limit)
    Pats(:,:,f) = nansum(fullPats,3);
    
elseif strncmpi(param.dot_occ,'m',1) %average of all dots overlapping each pixel
    fullPats = nansum(fullPats,3);
    frac_sum = nansum(frac_array,3);
    frac_sum(frac_sum<1) = 1; %to prevent purely edge pixels from getting scaled back up
    Pats(:,:,f) = fullPats./frac_sum;
    
elseif strncmpi(param.dot_occ,'c',1) %only shows closest dot (edge pixels will be partially occluded)
    %sort by dot closeness
    [~, order] = sort(dots_rho(vdot_idx)); 
    %order = flipud(order); %to sort by farthest
    fullPats = fullPats(:,:,order);
    frac_array = frac_array(:,:,order);
    
    %find indices of closest dots for each pixel (can see past edge pixels, to an extent)
    frac_array(isnan(frac_array)) = 0;
    cum_frac_array = cumsum(frac_array,3);
    closest_idx = circshift(cum_frac_array<1,1,3);
    closest_idx(:,:,1) = 1;
   
    %take average of only the closest dots
    frac_array(~closest_idx) = nan;
    fullPats(~closest_idx) = nan;
    fullPats = nansum(fullPats,3);
    frac_sum = nansum(frac_array,3);
    frac_sum(frac_sum<1) = 1;
    Pats(:,:,f) = fullPats./frac_sum;
end


%% move dots for next frame
if f<num_frames
    if strncmpi(param.motion_type,'r',1)
        dots_phi = dots_phi + true_step_size;%for rotation, motion takes place through the phi coordinate
        [dots_x, dots_y, dots_z] = sphere2cart(dots_phi,dots_theta,dots_rho);
        
    elseif strncmpi(param.motion_type,'t',1)
        dots_z = dots_z + true_step_size; %for translation, motion takes place through z
        dots_z(dots_z>view_rad) = dots_z(dots_z>view_rad) - 2*view_rad;
        [dots_phi, dots_theta, dots_rho] = cart2sphere(dots_x, dots_y, dots_z);
        
    else
        dots_theta = dots_theta + true_step_size; %for exp-con, motion takes place through theta
        dots_theta(dots_theta>pi) = dots_theta(dots_theta>pi) - pi;
        
        %switch dot trajectories for dots that have completed their current path
        e_idx_change = dots_theta>max_theta(traj_idx); %indices of dots done with their current trajectory
        dots_theta(e_idx_change) = dots_theta(e_idx_change) - pi/2; %start theta of other trajectory
        traj_idx([e_idx_change'; e_idx_change']) = ~traj_idx([e_idx_change'; e_idx_change']);
        traj_idx = logical(traj_idx);   
        dots_phi = dots_2_phi(traj_idx); %switch to phi of other trajectory
        
        [dots_x, dots_y, dots_z] = sphere2cart(dots_phi,dots_theta,dots_rho);
    end
end


%% end for loop
end

 
%% map pattern to desired brightness range
Pats(isnan(Pats)) = 0;
Pats = Pats + param.levels(2);
Pats(Pats>(2^param.gs_val - 1)) = 2^param.gs_val - 1;


end