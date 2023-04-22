close all; clc; clear;
load('usborder.mat','x','y','xx','yy');
nCities = 11;  % O(N^2)
citiesLat = zeros(nCities,1);
citiesLon = citiesLat;
order = 1:nCities;
nPopulation = 20; % population size
iterations = 200;

% Check if number of cities is correct
if nCities <= 2
    error("Number of Cities must be greater than 2");
end

% Randomly generated location of cities
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

% Predetermined location of cities
% best = 856.9473025397323
%citiesLat = [38.9690580717583; 84.2709674412475; 287.294543794855; 287.192917524641; 223.196607969350; 196.987889163089; 30.0424025895184; 208.753758269769; 132.184563537263; 56.8749087617592];
%citiesLon = [274.099380985567; 164.517574242290; 289.501672024584; 146.127318968130; 118.275878840716; 52.1848196556570; 247.213890669861; 95.8127445381973; 115.085978670810; 147.439554340681];

% generating an initial population
chromosome = zeros(nPopulation, nCities);
for i = 1:nPopulation
    perm = randperm(length(order));
    chromosome(i,:) = order(perm);
end

[dist] = objective_function(citiesLat, citiesLon, chromosome, nPopulation, nCities);
[probability] = selection_probability(dist, nPopulation, chromosome);
[bestOffsprings] = probability(end,2:end); % Best offspring from every iteration

for k=1:iterations
    
    mate = [];
    for i = 1:round(nPopulation/2) 
        [mate] = [mate; randsample(1:nPopulation, 2, true, probability(:,1))];
    end
    
    [cros, tablica, tablica2] = crossover(mate, probability, nCities);
    [chromosome_mut] = mutation(cros, nCities);
    [dist] = objective_function(citiesLat, citiesLon, chromosome_mut, nPopulation, nCities);
    TEMP = [dist, chromosome_mut];
    [min_value, min_index] = min(TEMP(:,1));
    [bestOffsprings] = [bestOffsprings; TEMP(min_index,:)];

end


[bestOffsprings] = sortrows(bestOffsprings, -1);

% Drawing lines between points
% [cities] = [citiesLat'; citiesLon'];
% for j = 1:length(bestOffsprings)
%     clf;
%     plot(214*x,300*y)
%     hold on
%     scatter(citiesLat,citiesLon,'Filled')
%     for k = 3:nCities+1
%         hold on
%         plot([cities(1,bestOffsprings(j,k-1)), cities(1,bestOffsprings(j,k))], ...
%             [cities(2,bestOffsprings(j,k-1)), cities(2,bestOffsprings(j,k))], 'b-');
%         drawnow;
%         pause(0.01);
%     end
%     plot([cities(1,bestOffsprings(j,end)), cities(1,bestOffsprings(j,2))], ...
%             [cities(2,bestOffsprings(j,end)), cities(2,bestOffsprings(j,2))], 'b-');
%     drawnow;
%     hold off;
% end

[cities] = [citiesLat'; citiesLon'];

clf;
plot(214*x,300*y)
hold on
scatter(citiesLat,citiesLon,'Filled')
for k = 3:nCities+1
    hold on
    plot([cities(1,bestOffsprings(end,k-1)), cities(1,bestOffsprings(end,k))], ...
        [cities(2,bestOffsprings(end,k-1)), cities(2,bestOffsprings(end,k))], 'b-');
    drawnow;
    pause(0.01);
end
plot([cities(1,bestOffsprings(end,end)), cities(1,bestOffsprings(end,2))], ...
    [cities(2,bestOffsprings(end,end)), cities(2,bestOffsprings(end,2))], 'b-');
    drawnow;
hold off;

figure;
plot(1:length(bestOffsprings), bestOffsprings(:,1), 'b-', 'LineWidth', 2);
hold on
plot(length(bestOffsprings), bestOffsprings(end,1), 'r.', 'MarkerSize', 10);

% Function definitions:
% I
function [dist] = objective_function(citiesLat, citiesLon, chromosome, nPopulation, nCities)
    dist_mx = zeros(nPopulation, nCities);
    for i = 1:nPopulation
        dist_mx(i,1:end-1) = sqrt((citiesLat(chromosome(i,2:end)) - citiesLat(chromosome(i,1:end-1))).^2 + (citiesLon(chromosome(i,2:end)) - citiesLon(chromosome(i,1:end-1))).^2);
        dist_mx(i,end) = sqrt((citiesLat(chromosome(i,1)) - citiesLat(chromosome(i,end))).^2 + (citiesLon(chromosome(i,1)) - citiesLon(chromosome(i,end))).^2);
    end 
    dist = sum(dist_mx, 2);
end

% II
function [probability] = selection_probability(dist, nPopulation, chromosome)
    list = [dist, chromosome];
    list = sortrows(list, 1, 'descend');
    rank_sum = nPopulation .* (nPopulation + 1) ./ 2;
    probability = (1:nPopulation) / rank_sum;
    probability = [probability', list];
end

% III
function [cros, tablica, tablica2] = crossover(mate, probability, nCities)

    cros = [];
    for i = 1:size(mate,1)
        [cros] = [cros; probability(mate(i,1),3:end); probability(mate(i,2),3:end)];
    end

    start = randi(nCities); % random start point
    dl = randi([0,nCities - start]); % random number of indexes to take, 0 or sth.
    for i = start:start+dl
        temp = cros(1:2:end,i);
        cros(1:2:end,i) = cros(2:2:end,i);
        cros(2:2:end,i) = temp;
    end
    
    %disp(start)
    %disp(dl)
    %disp('----------')

    T1 = [cros(1:2:end,1:start-1), zeros(length(cros(1:2:end,1)), length(start:start+dl)), cros(1:2:end,start+dl+1:end)];
    disp(size(T1))
    disp(T1)
    disp(length(T1))
    tablica = cell(size(T1,1), 1);
    % watch out: length(X) returns the length of the largest array dimension in X
    for k = 1:size(T1,1)
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
    tablica2 = cell(size(T2,1), 1);
    for k = 1:size(T2,1)
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
    %disp(cros)
    for i = 1:length(tablica)
        for j = 1:length(tablica{i}(:))
            temp = cros(2.*i-1, tablica{i}(j));
            cros(2.*i-1, tablica{i}(j)) = cros(2.*i, tablica2{i}(j));
            cros(2.*i, tablica2{i}(j)) = temp;
        end
    end
end

% IV
function [chromosome_mut] = mutation(cros, nCities)
    prob_m = 1/nCities;
    chromosome_mut = [];
    rng('shuffle');
    for i = 1:length(cros)
        if rand() <= prob_m 
            locus1 = randi(nCities);
            locus2 = randi(nCities);
            while locus2 == locus1
                locus2 = randi(nCities);
            end
            %disp([locus1, locus2])
            temp = cros(i,locus1);
            cros(i,locus1) = cros(i,locus2);
            cros(i,locus2) = temp;
            [chromosome_mut] = [chromosome_mut; cros(i,:)];
        else
            %disp('no mutation')
            [chromosome_mut] = [chromosome_mut; cros(i,:)];
            continue
        end
    end
end