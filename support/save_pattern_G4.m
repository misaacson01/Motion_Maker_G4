function save_pattern_G4(Pats, param, save_loc, filename)
% FUNCTION save_pattern_G4(Pats, param, save_loc, filename)
% 
% Saves the Pats variable to both .mat and .pat files, the former of which
% can be easily read back into Matlab and the latter a binary file which is
% used to display the pattern on the G4 LED arena.
%
% inputs:
% Pats: array of brightness values for each pixel in the arena
% param: full parameters of the input pattern to be stored in .mat file
% save_loc: directory to store the pattern files
% filename: desired name of the .mat pattern file

pattern.Pats = Pats;
pattern.x_num = length(Pats(1,1,:,1));
pattern.y_num = length(Pats(1,1,1,:));
pattern.gs_val = param.gs_val; 
pattern.stretch = param.stretch;
pattern.param = param; %store full pattern parameters

%get the vector data for each pattern through function make_pattern_vector_g4
if exist('make_pattern_vector_g4','file')
    pattern.data = make_pattern_vector_g4(pattern);
else
    disp('could not save binary .pat file; missing script from PControl');
end

%save the mat file
if ~exist(save_loc,'dir')
    mkdir(save_loc)
end
matFileName = fullfile(save_loc, filename);
if exist(matFileName,'file')
    error('pattern already exists in save folder with that name')
end
save(matFileName, 'pattern');

%save the corresponding binary pat file
if exist('make_pattern_vector_g4','file')
    patFileName = fullfile(save_loc, [num2str(param.ID,'%04d') '.pat']);
    fileID = fopen(patFileName,'w');
    fwrite(fileID, pattern.data);
    fclose(fileID);
end

end