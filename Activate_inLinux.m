
% First, you need specify the location of the CRESLIDE 1.0Beta executable
addpath('./iHydroSlide3D_src'); 

% Second, you need provide a project file, which provides all necessary
% paths, files, initial conditions, parameter values and other information.
tic
prediction = Main_program('ControlFile.ProjectLinux');
toc