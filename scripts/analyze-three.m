#!/bin/octave

folder='../output/six';

subfolder = cell(3);

subfolder(1) = 'jain_neal_split';
subfolder(2) = 'algorithm8';
subfolder(3) = 'triadic';

Mavg = zeros(3, 3);

for i=1:3
	fname=strcat(folder, '/', subfolder(i){}, '/', 'purity.txt');
	fd = fopen(fname);
	[C,k] = textscan(fd, '%f');
	fclose(fd);
	M1 = cell2mat(C);

	fname=strcat(folder, '/', subfolder(i){}, '/', 'rand.txt');
	fd = fopen(fname);
	[C,k] = textscan(fd, '%f');
	fclose(fd);
	M2 = cell2mat(C);

	fname=strcat(folder, '/', subfolder(i){}, '/', 'adjusted_rand.txt');
	fd = fopen(fname);
	[C,k] = textscan(fd, '%f');
	fclose(fd);
	M3 = cell2mat(C);

	M=[M1 M2 M3];

	Mavg(i,:) =[mean(M)];
end

Mavg

%ylim([0 1])

%data=[mean(M/3); mean(M/2); mean(M)]';

%title("Clustering performance");
%figure(1)
%w=0.8;
%bar(data, w)
%set(gca,'XTick', 1:3, 'XTickLabel',{'Purity', 'Rand Index', 'Adjusted Rand Index'})
%legend("Auxiliary variable sampler", "Jain-Neal split-merge sampler", "Triadic sampler");
%hold on
%errorbar(mean(M),var(M))

%hold off;
