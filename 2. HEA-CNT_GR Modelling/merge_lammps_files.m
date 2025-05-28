function merge_lammps_files(file1, file2, outputFile)
    % Read first LAMMPS file (HEA)
    [num_atoms1, atom_data1, masses1, box_dim1, num_types1] = read_lammps_file(file1);
    
    % Read second LAMMPS file (HCS)
    [num_atoms2, atom_data2, masses2, box_dim2, num_types2] = read_lammps_file(file2);
    
    % Update atom types in the second file to avoid overlap
    atom_data2(:, 2) = atom_data2(:, 2) + num_types1; % Shift atom types of second file
    
    % Combine atoms data
    merged_atom_data = [atom_data1; atom_data2];
    
    % Sort atoms by atom type
    merged_atom_data = sortrows(merged_atom_data, 2);
    
    % Update atom serial numbers after sorting
    merged_atom_data(:, 1) = 1:size(merged_atom_data, 1);
    
    % Update the total number of atom types
    total_types = num_types1 + num_types2;
    
    % Write the merged data to a new LAMMPS file
    write_lammps_file(outputFile, merged_atom_data, [masses1; masses2], box_dim1, total_types);
    
    disp('Files merged successfully!');
end

function [num_atoms, atom_data, masses, box_dim, num_types] = read_lammps_file(filename)
    % Open LAMMPS file for reading
    fid = fopen(filename, 'r');
    
    % Check if the file was successfully opened
    if fid == -1
        error('Error: File %s could not be opened. Please check the file path or name.', filename);
    end
    
    atom_data = [];
    masses = [];
    box_dim = [];
    reading_atoms = false;
    reading_masses = false;
    num_atoms = 0;
    num_types = 0;
    
    % Read through the file line by line
    while ~feof(fid)
        line = fgetl(fid);
        
        % Get the number of atoms
        if contains(line, 'atoms')
            num_atoms = sscanf(line, '%d atoms');
        end
        
        % Get the number of atom types
        if contains(line, 'atom types')
            num_types = sscanf(line, '%d atom types');
        end
        
        % Get box dimensions
        if contains(line, 'xlo xhi')
            dims = sscanf(line, '%f %f xlo xhi');
            box_dim(1, :) = dims;
        elseif contains(line, 'ylo yhi')
            dims = sscanf(line, '%f %f ylo xhi');
            box_dim(2, :) = dims;
        elseif contains(line, 'zlo zhi')
            dims = sscanf(line, '%f %f zlo xhi');
            box_dim(3, :) = dims;
        end
        
        % Start reading masses section
        if contains(line, 'Masses')
            reading_masses = true;
            fgetl(fid); % Skip the header
            continue;
        end
        
        % Read mass values
        if reading_masses
            if isempty(line)
                reading_masses = false; % Stop reading masses
            else
                mass_data = sscanf(line, '%d %f');
                masses = [masses; mass_data'];
            end
        end
        
        % Start reading atom data
        if contains(line, 'Atoms')
            reading_atoms = true;
            fgetl(fid); % Skip header line
            continue;
        end
        
        % Extract atom data
        if reading_atoms
            atom = sscanf(line, '%d %d %f %f %f');
            if numel(atom) == 5 %~isempty(atom)
                atom_data = [atom_data; atom'];
            end
        end
    end
    
    fclose(fid);
end

function write_lammps_file(outputFile, atom_data, masses, box_dim, total_types)
    % Open output file for writing
    fid = fopen(outputFile, 'w');
    
    % Check if the output file is opened successfully
    if fid == -1
        error('Error: Output file %s could not be created.', outputFile);
    end
    
    % Write header information
    fprintf(fid, '# Merged LAMMPS data file\n');
    fprintf(fid, '%d atoms\n', size(atom_data, 1));
    fprintf(fid, '%d atom types\n\n', total_types);
    
    % Write box dimensions
    fprintf(fid, '%f %f xlo xhi\n', box_dim(1, 1), box_dim(1, 2));
    fprintf(fid, '%f %f ylo yhi\n', box_dim(2, 1), box_dim(2, 2));
    fprintf(fid, '%f %f zlo zhi\n\n', box_dim(3, 1), box_dim(3, 2));
    
    % Write masses section
    fprintf(fid, 'Masses\n\n');
    for i = 1:total_types
        fprintf(fid, '%d %f\n', i, masses(i, 2));  % Corrected to assign unique IDs to each mass
    end
    
    % Write atoms section
    fprintf(fid, '\nAtoms # atomic\n\n');
    for i = 1:size(atom_data, 1)
        fprintf(fid, '%d %d %f %f %f\n', atom_data(i, 1), atom_data(i, 2), atom_data(i, 3), atom_data(i, 4), atom_data(i, 5));
    end
    
    fclose(fid);
    disp(['LAMMPS data file saved as ', outputFile]);
end
