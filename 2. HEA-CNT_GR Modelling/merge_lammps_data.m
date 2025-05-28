clear all
clc
close all

file1 = sprintf('HEA.xyz');
file2 = sprintf('Generated_CNT(4,4,14.90).dat');
outputFile = sprintf('HEA_CNT(4,4,14.90).dat');
merge_lammps_files(file1, file2, outputFile);

