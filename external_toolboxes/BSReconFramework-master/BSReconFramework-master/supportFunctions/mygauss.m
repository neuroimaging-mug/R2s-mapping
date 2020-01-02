function out = mygauss(dim)
x=linspace(-3.5,3.5,dim);
% out = 1/sqrt(2*pi) * exp(-.5 * x.^2);
out = exp(-.5 * x.^2);
end