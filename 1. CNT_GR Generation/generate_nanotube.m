% MATLAB script for CNT generation and save to LAMMPS format
% MATLAB (R2023b)
% By Subrata Barman  
% The DOI will be added upon publication

close all 
clear all
clc;

% Input parameters for CNT generation
n = input('Enter the Chiral Index n>0 & integer : ');
m = input('Enter the Chiral Index m>0 & integer : ');
length_nm = input('Enter the nanotube length in nm : ');

bonds_on = 0; % Turn bonds calculation on/off (1 = on, 0 = off)
angles_on = 0; % Turn angles calculation on/off (1 = on, 0 = off)
dihedrals_on = 0; % Turn dihedrals calculation on/off (1 = on, 0 = off)
impropers_on = 0; % Turn impropers calculation on/off (1 = on, 0 = off)

% Validate inputs
if length_nm <= 0
    error('Nanotube length must be a positive value');
end
if n <= 0 || m < 0 || floor(n) ~= n || floor(m) ~= m
    error('n and m must be non-negative integers');
end
if m == 0 && n == 0
    error('n and m cannot both be zero');
end

% Define constants
a = 1.418;  % Carbon-carbon bond length in Angstroms
pi_val = pi;

% Calculate the chiral vector and radius
C = a * sqrt(3 * (n^2 + n * m + m^2));
R = C / (2 * pi_val);

% Calculate the unit cell length along the tube axis
d_R = gcd(2 * m + n, 2 * n + m); % Calculate greatest common divisor
L_cell = sqrt(3) * C / d_R;

% Calculate number of unit cells required for desired tube length
N_cell = ceil(length_nm * 10 / L_cell);  % Convert nm to Ångstroms

% Initialize arrays for atom coordinates
coords = [];
i = 1;

% Generate unit cell coordinates
pmin = 0;
pmax = ceil(n + (n + 2 * m) / d_R);
qmin = floor(-(2 * n + m) / d_R);
qmax = m;

for q = qmin:qmax
    for p = pmin:pmax
        % Calculate coordinates for the two basis atoms
        xprime1 = 3 * a^2 * (p * (2 * n + m) + q * (n + 2 * m)) / (2 * C);
        yprime1 = 3 * sqrt(3) * a^2 * (p * m - q * n) / (2 * C);
        xprime2 = xprime1 + 3 * a^2 * (n + m) / (2 * C);
        yprime2 = yprime1 - a^2 * sqrt(3) * (n - m) / (2 * C);

        phi1 = xprime1 / R;
        phi2 = xprime2 / R;

        % Transform to cylindrical coordinates
        x1 = R * cos(phi1);
        x2 = R * cos(phi2);
        y1 = R * sin(phi1);
        y2 = R * sin(phi2);
        z1 = yprime1;
        z2 = yprime2;

        % Store coordinates for the unit cell
        if 0 <= xprime1 && p * (2 * n + m) + q * (n + 2 * m) < 2 * (n^2 + n * m + m^2) && ...
           0 <= yprime1 && d_R * (p * m - q * n) < 2 * (n^2 + n * m + m^2)
            coords(i, :) = [x1, y1, z1];
            coords(i + 1, :) = [x2, y2, z2];
            i = i + 2;
        end
    end
end

% Expand coordinates for the full nanotube length
nanotube_coords = [];
for j = 0:(N_cell - 1)
    nanotube_coords = [nanotube_coords; coords(:, 1), coords(:, 2), coords(:, 3) + j * L_cell];
end

% Display basic nanotube properties
disp(['Chirality: (' num2str(n) ',' num2str(m) ')']);
disp(['Radius (Å): ' num2str(R)]);
disp(['Number of atoms: ' num2str(size(nanotube_coords, 1))]);
disp(['Length (nm): ' num2str(length_nm)]);

% Optional: add bonds, angles, dihedrals, and impropers if requested
if bonds_on
    % Placeholder for bonds calculation (implement if needed)
    disp('Bonds calculation enabled (placeholder).');
end
if angles_on
    % Placeholder for angles calculation (implement if needed)
    disp('Angles calculation enabled (placeholder).');
end
if dihedrals_on
    % Placeholder for dihedrals calculation (implement if needed)
    disp('Dihedrals calculation enabled (placeholder).');
end
if impropers_on
    % Placeholder for impropers calculation (implement if needed)
    disp('Impropers calculation enabled (placeholder).');
end

disp(['CNT Atoms: ', num2str(size(nanotube_coords, 1))]);

% Writing header for LAMMPS data file
sigma = 3.57; 
%lattice constant of AlCoCrFeNi HEA can be changed to adjust the simulation box diementions 

% Variables to hold cell dimensions
box_size_x = 25; %input('Enter the X dimension of the simulation box : ');   %
box_size_y = 25; %input('Enter the Y dimension of the simulation box : ');   %
box_size_z = 50; %input('Enter the Z dimension of the simulation box : ');   %
box_size = [box_size_x*sigma box_size_y*sigma box_size_z*sigma];

xlo = 0; xhi = box_size(:,1);
ylo = 0; yhi = box_size(:,2);
zlo = 0; zhi = box_size(:,3);

% Center of the simulation box
box_center = [(xhi + xlo) / 2, (yhi + ylo) / 2, (zhi + zlo) / 2-(max(nanotube_coords(:,3))-min(nanotube_coords(:,3)))/2];

% The variable `shifted_coords` now holds all the atom coordinates for the centered nanotube
shifted_coords = zeros(size(nanotube_coords, 1), 3);
% Shift atoms to the center of the simulation box
for i = 1:size(nanotube_coords, 1)
    shifted_coords(i, :) = nanotube_coords(i, :) + box_center;
end

% Generating output in LAMMPS format
filename = sprintf('Generated_CNT(%d,%d,%.2f).dat', n, m, length_nm);  % Dynamic file name

fid = fopen(filename, 'wt');  % Open file to write in LAMMPS format

fprintf(fid, 'LAMMPS data file\n');
fprintf(fid, '%d atoms\n', size(nanotube_coords, 1));
fprintf(fid, '1 atom types\n\n');
fprintf(fid, '%.6f %.6f xlo xhi\n', xlo, xhi); % X boundaries
fprintf(fid, '%.6f %.6f ylo yhi\n', ylo, yhi); % Y boundaries
fprintf(fid, '%.6f %.6f zlo zhi\n', zlo, zhi); % Z boundaries
fprintf(fid, '\nMasses\n\n');
fprintf(fid, '1 12.0107  # CA (Carbon Atom)\n\n');

% Writing atoms section
fprintf(fid, 'Atoms # atomic\n\n');
for i = 1:size(shifted_coords, 1)
    fprintf(fid, '%d 1 %f %f %f\n', i, shifted_coords(i, 1), shifted_coords(i, 2), shifted_coords(i, 3));
end

fclose(fid);
disp(['LAMMPS data file saved as ', filename]);
