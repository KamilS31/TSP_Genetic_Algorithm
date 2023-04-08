close all; clc; clear;
load('usborder.mat','x','y','xx','yy');
nCities = 5;  % O(N^2)
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
for i = 1:nPopulation
    [mate] = [mate; randsample(1:nPopulation, 2, true, probability(:,1))];
end

[cros,y,z] = crossover(mate, probability, nCities);

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

function [cros,y,z] = crossover(mate, probability, nCities)
    [cros] = [probability(mate(1,1),3:end); probability(mate(1,2),3:end)];
    start = randi(nCities); % random start point
    dl = randi([0,nCities - start]); % random number of indexes to take, 0 or sth.
    for i = start:start+dl
        temp = cros(1,i); % value of first vector
        cros(1,i) = cros(2,i);
        cros(2,i) = temp;
    end
    disp(start)
    disp(dl)
    disp('----------')
    T1 = [cros(1,1:start-1), zeros(1,length(start:start+dl)), cros(1,start+dl+1:end)];
    y = [];
    for j = start:start+dl
        [y] = [y, find(T1(1,:) == cros(1,j))];
    end
    szukana = start:start+dl;
    temp = cros(2, szukana);
    cros(2, szukana) = 0;
    z = [];
    for k = temp
        [z] = [z, find(cros(2,:) == k)];
    end
    cros(2, szukana) = temp;

    temp = cros(1,y); 
    cros(1,y) = cros(2,z);
    cros(2,z) = temp;
end






