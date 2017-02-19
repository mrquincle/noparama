mu=[0 0];
sigma=eye(2);
N=100;

Y0=mvnrnd(mu, sigma, N);

Y0=[Y0 ones(size(Y0,1),1)*0];

mu=[5 5];

Y1=mvnrnd(mu, sigma, N);

Y1=[Y1 ones(size(Y1,1),1)*1];

Y=[Y0;Y1]

save('-ascii', '../datasets/twogaussians.data', 'Y')
