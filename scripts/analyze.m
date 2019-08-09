#!/usr/bin/octave -q

if length(argv()) == 0
	error("Usage: analyze.m folder")
end

folder=argv(){1};

Mavg = zeros(3);

fname=strcat(folder, '/', 'purity.txt');
fd = fopen(fname);
[C,k] = textscan(fd, '%f');
fclose(fd);
M1 = cell2mat(C);

fname=strcat(folder, '/', 'rand.txt');
fd = fopen(fname);
[C,k] = textscan(fd, '%f');
fclose(fd);
M2 = cell2mat(C);

fname=strcat(folder, '/', 'adjusted_rand.txt');
fd = fopen(fname);
[C,k] = textscan(fd, '%f');
fclose(fd);
M3 = cell2mat(C);

M=[M1 M2 M3];

Mavg = [mean(M)];

Mavg

Mmedian = [median(M)];

Mmedian
