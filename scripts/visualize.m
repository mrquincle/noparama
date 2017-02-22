%clear all
clf
cmap = colormap;

%dirname = ['../output/20170215_13:05'];
%dirname = ['../output/20170215_13:10'];
dirname = ['../output/LATEST'];

fname = [dirname '/results.txt']
load(fname)

mu=mu';

K=size(mu,2)

for i = 1:K
	j=i*2;
	color = cmap(mod(5*j,63)+1,:); 

	R = sigma(:,:,i);
	mu_i = mu(:,i);

	% calculate the ellipse
	A     = chol ( R, "lower");
	theta = linspace (0, 2*pi, 1000);
	x     = mu_i + 2.5 .* A * [cos(theta); sin(theta)];

	% plot the ellipse
	hold on;
	plot(x(1,:), x(2,:), "r", "LineWidth", 2, 'color', color);

	% plot the data
	fname = [dirname '/results' num2str(i-1) '.txt'];
	data = load(fname)';
	plot(data(1,:), data(2,:), '.', 'MarkerSize', 20, 'color', color);
end


hold off

