function [R] = normrnd(mu_desired, std_desired, N)
%
% USAGE [R] = normrnd(mu_desired, std_desired, N)
% 
% normrnd returns a vector of random numbers with a given mean
% and standard deviation, where
% mu_desired is the desired mean
% std_desired is the desired standard deviation
% N is the number of row elements in the vector

R = randn(N,1);

% Measure
mu_tmp = mean(R);
std_tmp = std(R);

% Normalise and denormalise
R = (R - mu_tmp) / std_tmp;
R = (R * std_desired) + mu_desired;