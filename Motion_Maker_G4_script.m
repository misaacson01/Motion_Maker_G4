% Motion_Maker_G4_script
%
% Generates patterns in script format (rather than using the GUI)

%% user-defined pattern parameters
% all angles/distances/sizes in units of radians
% some parameters are only needed in certain cirumstances {specified by curly brace}

param.pattern_type = 'square grating'; %'square grating', 'sine grating', 'edge', 'starfield', or 'Off/On'
param.motion_type = 'rotation'; %'rotation', 'translation', or 'expansion-contraction'
param.pattern_fov = 'full-field'; %'full-field' or 'local'
param.arena_pitch = deg2rad(0); %angle of arena pitch (0 = straight ahead, positive values = pitched up)
param.gs_val = 4; %bits of intensity value (1 or 4)
param.levels = [15 0 7]; %brightness level of [1st bar (in grating) or advancing edge, 2nd bar or receding edge, background (mask)]
param.pole_coord = deg2rad([0 0]); %location of pattern pole [longitude, lattitude] {for pattern_fov=full-field}
param.motion_angle = deg2rad(0); %angle of rotation (0=rightward motion, positive values rotate the direction clockwise) {fov=local}
param.spat_freq = deg2rad(30); %spatial angle (in radians) before pattern repeats {for gratings and edge}
param.step_size = deg2rad(1.25); %amount of motion per frame (in radians) {for type~=off/on}
param.duty_cycle = 50; %percent of spat_freq taken up by first bar {for square gratings}
param.num_dots = 500; %number of dots in star-field {for type=starfield}
param.dot_radius = deg2rad(1.25); %radius of dots (in radians) {for starfield}
param.dot_size = 'static'; %'static' or 'distance-relative' {for starfield}
param.dot_occ = 'closest'; %how occluding dots are drawn ('closest', 'sum', or 'mean') {for starfield}
param.dot_level = 0; %0 = dot brightness set to 1st level; 1 and 2 = random brightness (0-1st; 0 or 1st) {for starfield}
param.view_radius = 1; %distance from center that dots can be seen {for starfield}
param.sa_mask = [0 0 pi 0]; %location, size, and direction of solid angle mask [longitude, lattitude, solid_angle, out/in]
param.long_lat_mask = [-pi pi -pi/2 pi/2 0]; %coordinates of lattitude/longitude mask [min-long, max-long, min-lat, max-lattitude, out/in]
param.aa_samples = 15; %# of samples taken to calculate the brightness of each pixel (1 or 15 suggested)
param.aa_poles = 1; %1=anti-aliases the poles of rotation/translation grating/edge stimuli by matching them to the duty cycle
param.back_frame = 0; %1=adds a frame (frame 1) uniformly at background (mask) level
param.flip_right = 0; %1=left-right flips the right half of the pattern
param.phase_shift = 0; %shifts the starting frame of pattern


%% generate pattern
[Pats, arena_phi, arena_theta, p_rad, param.true_step_size, param.rot180] = Motion_Maker_G4(param);
param.stretch = zeros(size(Pats,3),1); %stretch increases (within limits) the per-frame brightness -- zeros add no brightness


%% save pattern
%set save directory
save_dir = 'C:\matlabroot\patterns_G4';
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end
cd(save_dir);

%set pattern ID
ID = 1;
patfiles = ls('Pattern*.mat');
if ~isempty(patfiles)
    ID = size(patfiles,1)+1;
end
patName = ['Pattern_' num2str(ID, '%04d') '_G4.mat'];

%store pattern and all parameters
pattern.Pats = Pats;
pattern.gs_val = param.gs_val;
pattern.param = param;

%save the mat file
matFileName = fullfile(save_dir, patName);
if exist(matFileName,'file')
    error('pattern already exists in save folder with that name')
end
save(matFileName, 'pattern');
disp(['Pattern saved as "' matFileName '"'])