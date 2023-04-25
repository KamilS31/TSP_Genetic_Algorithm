tic
close all; clc; clear;

load('usborder.mat','x','y');
nCities = 100;  % O(N^2)
nPopulation = 500; % population size

% Check if number of cities is correct
if nCities <= 2
    error("Number of Cities must be greater than 2");
end

% Randomly generated location of cities
citiesLat = zeros(nCities,1);
citiesLon = citiesLat;
n = 1;
while (n <= nCities)
    x_point = 299.*rand() + 1;
    y_point = 299.*rand() + 1;
    if inpolygon(x_point,y_point,214*x,300*y) % Test if inside the border
        citiesLat(n) = x_point;
        citiesLon(n) = y_point;
        n = n+1;
    end
end


% Predetermined location of cities
% best = 1173.6020807052705
% citiesLat = [147.097255037233; 163.946552996298; 68.1776068040153; 38.2220966073422; 14.1449515079153; 131.439763983590; 204.896352790044; ...
%             122.780804496441; 70.8742781216051; 69.4075876162935; 183.293773235291; 108.813428622198; 258.738114647446; 128.843646504298; ...
%             74.3793072084819; 143.706484353342; 37.2563276313789; 209.572743876007; 192.058278041047; 242.656587639574; 182.430991857692; ...
%             52.9424378203251; 131.323359065516; 175.837947846222; 244.832345837630; 194.777371154247; 271.632153331198; 129.124425212286; ...
%             49.1480539376800; 150.430028618058];
% citiesLon = [180.256744225672; 91.3534234940900; 197.343657214707; 153.058873492159; 215.591710331776; 236.090134189278; 125.018702954118; ...
%             112.255029074254; 192.691665933279; 202.268176211748; 190.169386617769; 124.915461188450; 86.7874944957495; 69.4477471852106; ...
%             243.017185149319; 210.168506582939; 246.438488501040; 86.7244641535839; 140.894722450234; 119.695129551657; 51.2782918074025; ...
%             138.988969738223; 125.743079836118; 92.4902549966836; 152.083468775689; 218.069099491905; 234.216911095138; 142.106874588722; ...
%             206.103330261466; 116.780288026215];

% generating an initial population
chromosome = zeros(nPopulation, nCities);
for i = 1:nPopulation
    perm = randperm(nCities);
    chromosome(i,:) = perm;  
end

[dist] = objective_function(citiesLat, citiesLon, chromosome, nPopulation, nCities);
[probability] = selection_probability(dist, nPopulation, chromosome);
% Best offspring from every iteration
[bestOffsprings] = probability(end,2:end); % Best offspring from initial generation

lack_imp = 0;
counter = 0;
score = probability(end,2);

while lack_imp < 500 % it not necessarily go N times
    mate = [];
    for i = 1:round(nPopulation/2) 
        [mate] = [mate; randsample(1:nPopulation, 2, true, probability(:,1))];
    end

    [cros, tablica, tablica2] = crossover(mate, probability, nCities, bestOffsprings);
    [chromosome_mut] = mutation(cros, nCities);
    [dist] = objective_function(citiesLat, citiesLon, chromosome_mut, nPopulation, nCities);

    TEMP = [dist, chromosome_mut];
    [min_value, min_index] = min(TEMP(:,1));
    [bestOffsprings] = [bestOffsprings; TEMP(min_index,:)];

    [probability] = selection_probability(dist, nPopulation, chromosome_mut);
    
    
    disp(score)
    disp(bestOffsprings(end,1))
    if score < bestOffsprings(end,1) % lack of improvement
        lack_imp = lack_imp + 1;
        format long;
        disp(lack_imp)
    else  
        if counter == 50
            break
        elseif score - bestOffsprings(end,1) < 0.001 
            counter = counter + 1;
        end
        score = bestOffsprings(end,1);
        lack_imp = 0;
    end
end

[bestOffsprings] = sortrows(bestOffsprings, -1);

% Drawing lines between points
[cities] = [citiesLat'; citiesLon'];
clf;
plot(214*x,300*y)
hold on
scatter(citiesLat,citiesLon,'Filled')
scatter(citiesLat(bestOffsprings(end,2),:), citiesLon(bestOffsprings(end,2)),'Filled', 'green')
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
plot(1:size(bestOffsprings,1), bestOffsprings(:,1), 'b-', 'LineWidth', 2);
xlabel('Number of generations');
ylabel('value of objective function for offspiring')
hold on
plot(size(bestOffsprings,1), bestOffsprings(end,1), 'rx', 'LineWidth', 1.5, 'MarkerSize', 8);

disp(['Percentage of searching the state space: ', num2str(length(bestOffsprings)*100*100/factorial(nCities)) '%'])
disp(['Shortest found route: ', num2str(bestOffsprings(end,1))])
disp(' ')

toc


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
function [cros, tablica, tablica2] = crossover(mate, probability, nCities, bestOffsprings)

    cros = zeros(size(mate,1)*2, nCities);
    for i = 1:size(mate,1)
        cros((i-1)*2+1:i*2, :) = [probability(mate(i,1),3:end); probability(mate(i,2),3:end)];
    end
    % The best one to next generation:
    rng('shuffle');
    cros(randi(size(cros,2)),:) = bestOffsprings(end,2:end);

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
        for j = 1:numel(tablica{i}) % length(tablica{i}(:))
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
    for i = 1:size(cros,1)
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







