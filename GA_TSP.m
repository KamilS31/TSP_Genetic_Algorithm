close all; clc; clear;
load('usborder.mat','x','y','xx','yy');
nCities = 50;  % O(N^2)
citiesLat = zeros(nCities,1);
citiesLon = citiesLat;
order = 1:nCities;
nPopulation = 20; % population size

% generating an initial population
chromosome = zeros(nPopulation, nCities);
for i = 1:nPopulation
    perm = randperm(length(order));
    chromosome(i,:) = order(perm);
end

n = 1;
while (n <= nCities)
    xp = 299.*rand() + 1;
    yp = 299.*rand() + 1;
    if inpolygon(xp,yp,214*x,300*y) % Test if inside the border
        citiesLat(n) = xp;
        citiesLon(n) = yp;
        n = n+1;
    end
end
plot(214*x,300*y)
hold on
scatter(citiesLat,citiesLon,'Filled')

[dist] = objective_function(citiesLat, citiesLon, chromosome, nPopulation, nCities);
[probability] = selection_probability(dist, nPopulation, chromosome);

mate = [];
for i = 1:5
    [mate] = [mate; randsample(1:nPopulation, 2, true, probability(:,1))];
end


% Function definitions:
function [dist] = objective_function(citiesLat, citiesLon, chromosome, nPopulation, nCities)
    dist_mx = zeros(nPopulation, nCities);
    for i = 1:nPopulation
        dist_mx(i,1:end-1) = sqrt((citiesLat(chromosome(i,2:end)) - citiesLat(chromosome(i,1:end-1))).^2 + (citiesLon(chromosome(i,2:end)) - citiesLon(chromosome(i,1:end-1))).^2);
        dist_mx(i,end) = sqrt((citiesLat(chromosome(i,1)) - citiesLat(chromosome(i,end))).^2 + (citiesLon(chromosome(i,1)) - citiesLon(chromosome(i,end))).^2);
    end 
    dist = sum(dist_mx, 2);
end

function [probability] = selection_probability(dist, nPopulation, chromosome)
    list = [dist, chromosome];
    list = sortrows(list, 1, 'descend');
    rank_sum = nPopulation .* (nPopulation + 1) ./ 2;
    probability = (1:nPopulation) / rank_sum;
    probability = [probability', list];
end





