close all; clc; clear;
load('usborder.mat','x','y','xx','yy');
nCities = 20;  % O(N^2)
citiesLat = zeros(nCities,1);
citiesLon = citiesLat;

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

order = 1:nCities;

dist = objective_function(citiesLat,citiesLon,order);

% Function definitions:
function dist = objective_function(citiesLat, citiesLon,order)
    dist = sqrt((citiesLat(order(2:end)) - citiesLat(order(1:end-1))).^2 + (citiesLon(order(2:end)) - citiesLon(order(1:end-1))).^2);
    A = sqrt((citiesLat(order(1)) - citiesLat(order(end))).^2 + (citiesLon(order(1)) - citiesLon(order(end))).^2);
    dist = sum([dist', A]);
end

