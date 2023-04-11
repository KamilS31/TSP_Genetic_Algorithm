close all; clc; clear;
load('usborder.mat','x','y','xx','yy');
nCities = 10;  % O(N^2)
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

[cros, tablica, tablica2] = crossover(mate, probability, nCities);

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

function [cros, tablica, tablica2] = crossover(mate, probability, nCities)

    [cros] = [probability(mate(1:end,1),3:end); probability(mate(1:end,2),3:end)];
    start = randi(nCities); % random start point
    dl = randi([0,nCities - start]); % random number of indexes to take, 0 or sth.

    for i = start:start+dl
        temp = cros(1:2:end,i); % value of first vector
        cros(1:2:end,i) = cros(2:2:end,i);
        cros(2:2:end,i) = temp;
    end

    disp(start)
    disp(dl)
    disp('----------')

    T1 = [cros(1:2:end,1:start-1), zeros(length(cros(1:2:end,1)),length(start:start+dl)), cros(1:2:end,start+dl+1:end)];
    tablica = cell(length(T1), 1);
    for k = 1:length(T1)
        idx = [];
        for j = start:start+dl
            [idx] = [idx, find(T1(k,:) == cros(2.*k-1,j))];
        end

        if ~isempty(idx)
            tablica{k} = idx;
        elseif isempty(idx) && isempty(tablica{k})
            tablica{k} = {};
        end
    end

    T2 = [cros(2:2:end,1:start-1), zeros(length(cros(2:2:end,1)),length(start:start+dl)), cros(2:2:end,start+dl+1:end)];
    tablica2 = cell(length(T2), 1);
    for k = 1:length(T2)
        idx2 = [];
        for j = start:start+dl
            [idx2] = [idx2, find(T2(k,:) == cros(2.*k,j))];
        end

        if ~isempty(idx2)
            tablica2{k} = idx2;
        elseif isempty(idx2) && isempty(tablica2{k})
            tablica2{k} = {};
        end
    end

    for i = 1:length(tablica)
        for j = 1:length(tablica{i}(:))
            temp = cros(2.*i-1, tablica{i}(j));
            cros(2.*i-1, tablica{i}(j)) = cros(2.*i, tablica2{i}(j));
            cros(2.*i, tablica2{i}(j)) = temp;
        end
    end
end
