clear all
clf
cmap = colormap;

K=5

mu=zeros(2,K);
sigma=zeros(2,2,K);

%mu(:,i) = [-1.2636 1.02329 ]; sigma(:,:,i) = [ 0.1272 0; 0 0.1272]; i++;
%mu(:,i) = [ 3.1578 6.12903 ]; sigma(:,:,i) = [ 0.627089 0; 0 0.627089]; i++;
%mu(:,i) = [  -1.793 0.949772 ]; sigma(:,:,i) = [ 0.314666 0; 0 0.314666]; i++;
%mu(:,i) = [0.712294 0.244961 ]; sigma(:,:,i) = [ 0.597651 0; 0 0.597651]; i++;
%mu(:,i) = [-0.535471   -1.5184 ]; sigma(:,:,i) = [ 0.498098 0; 0 0.498098]; i++;
%mu(:,i) = [6.53592 5.09471 ]; sigma(:,:,i) = [ 2.01896 0; 0 2.01896]; i++;
%mu(:,i) = [4.7638 4.8292 ]; sigma(:,:,i) = [ 0.31972 0; 0 0.31972]; i++;

mu(:,1) = [0.0123209 -0.773208];
sigma(:,:,1) = [1.36157 0; 0 1.36157];
mu(:,2) = [-0.420987 1.15129];
sigma(:,:,2) = [0.401986 0; 0 0.401986];
mu(:,3) = [3.78765 4.69144];
sigma(:,:,3) = [0.414196 0; 0 0.414196];
mu(:,4) = [5.03733 6.07106];
sigma(:,:,4) = [0.776407 0; 0 0.776407];
mu(:,5) = [5.94044 4.36555];
sigma(:,:,5) = [0.444987 0; 0 0.444987];


for i = 1:K
	j=i*2;
	color = cmap(mod(5*j,63)+1,:); 

	R = sigma(:,:,i);
	mu_i = mu(:,i)';

	% calculate the ellipse
	A     = chol ( R, "lower");
	theta = linspace (0, 2*pi, 1000);
	x     = mu_i' + 2.5 .* A * [cos(theta); sin(theta)];

	% plot the ellipse
	hold on;
	plot(x(1,:), x(2,:), "r", "LineWidth", 2);

	% plot the data
	fname = ['../output' num2str(i-1) '.txt'];
	data = load(fname)';
	plot(data(1,:), data(2,:), '.', 'MarkerSize', 20, 'color', color);
end


hold off

