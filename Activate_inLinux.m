
% First, you need specify the location of the CRESLIDE 1.0Beta executable
addpath('./Hydro_Geo_Src'); 

% Second, you need provide a project file, which provides all necessary
% paths, files, initial conditions, parameter values and other information.
prediction=main_function('ControlFile.ProjectLinux');