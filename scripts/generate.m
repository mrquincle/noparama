mu=[0 0];
sigma=eye(2);
N=100;

Y0=mvnrnd(mu, sigma, N);

mu=[5 5];

Y1=mvnrnd(mu, sigma, N);

Y=[Y0;Y1];

save('-ascii', 'output.data', 'Y')



