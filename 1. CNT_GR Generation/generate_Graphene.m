% MATLAB script for Graphene generation and save to LAMMPS format
% MATLAB (R2023b)
% By Subrata Barman  
% The DOI will be added upon publication

close all; 
clear all; 
clc;

% Input parameters for Graphene generation
lx = input('Enter graphene sheet width in nm: ');
ly = input('Enter graphene sheet length in nm: ');
type = input('Enter edge type (armchair/zigzag): ', 's');
nlayers = 1; %input('Enter number of graphene layers: ');

% Validate inputs
if lx <= 0 || ly <= 0
    error('Graphene edge lengths must be positive values.');
end
if ~strcmp(type, 'armchair') && ~strcmp(type, 'zigzag')
    error('Edge type must be either "armchair" or "zigzag".');
end
if nlayers < 1 || floor(nlayers) ~= nlayers
    error('Number of layers must be a positive integer.');
end

% Define constants
a = 1.418;  % Carbon-carbon bond length in Ångstroms
pi_val = pi;

% Calculate unit cell parameters
if strcmp(type, 'armchair')
    Lx_cell = 2 * a * sin(60 * pi_val / 180);
    Ly_cell = 3 * a;
else
    Lx_cell = 3 * a;
    Ly_cell = 2 * a * sin(60 * pi_val / 180);
end

Nx_cell = ceil(lx * 10 / Lx_cell); % Convert nm to Å
Ny_cell = ceil(ly * 10 / Ly_cell);

% Initialize atom coordinates
coords = [];

% Generate graphene unit cell coordinates
if strcmp(type, 'armchair')
    r1 = [0, 0, 0];
    r2 = [-a * sin(60 * pi_val / 180), a * cos(60 * pi_val / 180), 0];
    r3 = r2 + [0, a, 0];
    r4 = [0, 2 * a, 0];
    l_shift = [0, a, 0];
else
    r1 = [0, 0, 0];
    r2 = [-a * cos(60 * pi_val / 180), a * sin(60 * pi_val / 180), 0];
    r3 = [a, 0, 0];
    r4 = r2 + [2 * a, 0, 0];
    l_shift = [a / 2, a * sin(60 * pi_val / 180), 0];
end

% Generate full graphene sheet coordinates
for k = 0:(nlayers - 1)
    for j = 0:(Ny_cell - 1)
        for i = 0:(Nx_cell - 1)
            r_shift = [i * Lx_cell, j * Ly_cell, k * 3.35];
            if mod(k, 2) ~= 0
                r_shift = r_shift + l_shift;
            end
            coords = [coords; r1 + r_shift; r2 + r_shift; r3 + r_shift; r4 + r_shift];
        end
    end
end


% Choose plane of graphene generation
disp('Select plane for graphene generation:');
disp('1: XY plane');
disp('2: XZ plane');
disp('3: YZ plane');
plane_choice = input('Enter your choice (1/2/3): '); % Enter 2 for XZ plane

% Adjust coordinates to the selected plane
if plane_choice == 2  % XZ Plane
    coords = coords(:, [1, 3, 2]); % Swap Y and Z
elseif plane_choice == 3  % YZ Plane
    coords = coords(:, [3, 1, 2]); % Swap X and Z
end

% Optional rotation
apply_rotation = input('Do you want to rotate graphene? (1 = Yes, 0 = No): ');
if apply_rotation
    rotation_axis = input('Enter rotation axis (x/y/z): ', 's');
    rotation_angle = input('Enter rotation angle in degrees: ');
end

% Apply rotation if requested
if apply_rotation
    coords = rotate_graphene(coords, rotation_axis, rotation_angle);
end

% Writing header for LAMMPS data file
sigma = 3.57; 
%lattice constant of AlCoCrFeNi HEA can be changed to adjust the simulation box diementions 

% Get user-defined simulation box size
box_size_x = 25; %input('Enter the X dimension of the simulation box: ');
box_size_y = 25; %input('Enter the Y dimension of the simulation box: ');
box_size_z = 50; %input('Enter the Z dimension of the simulation box: ');
box_size = [box_size_x * sigma, box_size_y * sigma, box_size_z * sigma];

xlo = 0; xhi = box_size(1);
ylo = 0; yhi = box_size(2);
zlo = 0; zhi = box_size(3);

% Center graphene in the simulation box
box_center = [(xhi + xlo) / 2 - (max(coords(:,1)) - min(coords(:,1))) / 2,...
              (yhi + ylo) / 2 - (max(coords(:,2)) - min(coords(:,2))) / 2,...
              (zhi + zlo) / 2 - (max(coords(:,3)) - min(coords(:,3))) / 2];

% The variable `shifted_coords` now holds all the atom coordinates for the centered nanotube
shifted_coords = coords + box_center;

% Display graphene properties
disp(['Edge type: ', type]);
disp(['Graphene dimensions (nm): ', num2str(lx), ' x ', num2str(ly)]);
disp(['Number of layers: ', num2str(nlayers)]);
disp(['Total atoms: ', num2str(size(shifted_coords, 1))]);

% Write LAMMPS data file
filename = sprintf('Generated_Graphene_%s_%.1f,%.1f.dat', type, lx, ly);
fid = fopen(filename, 'wt');

fprintf(fid, 'LAMMPS data file\n');
fprintf(fid, '%d atoms\n', size(shifted_coords, 1));
fprintf(fid, '1 atom types\n\n');
fprintf(fid, '%.6f %.6f xlo xhi\n', xlo, xhi);
fprintf(fid, '%.6f %.6f ylo yhi\n', ylo, yhi);
fprintf(fid, '%.6f %.6f zlo zhi\n\n', zlo, zhi);
fprintf(fid, 'Masses\n\n');
fprintf(fid, '1 12.0107  # Carbon Atom (C)\n\n');

fprintf(fid, 'Atoms # atomic\n\n');
for i = 1:size(shifted_coords, 1)
    fprintf(fid, '%d 1 %f %f %f\n', i, shifted_coords(i, 1), shifted_coords(i, 2), shifted_coords(i, 3));
end

fclose(fid);
disp(['LAMMPS data file saved as ', filename]);
