% Obtain absolute path to this script
this_path = fileparts(mfilename('fullpath'));

% Add appropriate paths to matlab environment
addpath(genpath(fullfile(this_path, 'lib')));
addpath(genpath(fullfile(this_path, 'data_pipeline')));
addpath(genpath(fullfile(this_path, 'demos_and_examples')));

% Global variable with absolute path to the l1ktools/matlab folder
% Referenced within the library
MORTARPATH = this_path;