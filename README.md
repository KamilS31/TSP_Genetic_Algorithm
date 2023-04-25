# TSP_Genetic_Algorithm

## Issue to be solved

The travelling salesman problem asks the following question: "Given a list of cities and the distances between each pair of cities, what is the shortest possible route that visits each city exactly once and returns to the origin city?" In other words, the TSP seeks to find a Hamiltonian cycle in an undirected complete graph with weighted edges, where the goal is to minimize the sum of the edge weights along the cycle.

In this version of the problem, I will try to find the shortest route between points (cities) plotted on a simplified map of the United States.

## Probelm solving approach

To solve this problem I will use evolutionary algorithms and, more specifically, their subclass which is genetic algorithms.

The entire algorithm is based on several important steps:
1) Generation of the initial population.
```
chromosome = zeros(nPopulation, nCities);
for i = 1:nPopulation
    perm = randperm(nCities);
    chromosome(i,:) = perm;  
end
```
2) Create a condition for the end of evolution.
```
while lack_imp < N    % No improvement for N iterations
```
3) Evaluation of individuals in a population.
```
[dist] = objective_function(citiesLat, citiesLon, chromosome, nPopulation, nCities);
```
4) Calculating the probability of selection for crossover.
```
[probability] = selection_probability(dist, nPopulation, chromosome);
```
5) Crossovering of individuals.
```
[cros, tablica, tablica2] = crossover(mate, probability, nCities, bestOffsprings);
```
6) Mutation of individuals with a specific mutation probability equal to 1/number_points.
```
[chromosome_mut] = mutation(cros, nCities);
```
7) Evolution ends when the loop condition is not met (no improvement in the best solution after N new populations) or when the best solutions differ by less than 0.001 by M new populations.

## Technologies
* Matlab R2021b

## Implemented idea of crossover

(Note, within a single chromosome, genes must not repeat)

Initially, we draw pairs for crossbreeding. Once the pairs are drawn, the crossover proceeds as follows:
1) We draw any locus (number from 1 to the total number of cities) from the chromosome.
2) We draw the number of genes we will take for crossover (a number from 0 to the total number of cities minus the locus above)
3) We swap genes with the drawn locus between the two paired chromosomes.
4) When two identical genes appear in the chromosome after crossover, the repeated genes are swapped between the paired chromosomes. (Note, this swap occurs only with genes that occur outside the drawn locus).

The image below visualizes the whole process.

<div align="center">
  <img src="https://user-images.githubusercontent.com/127042515/234414638-519efafb-d36c-46e6-bbde-eae08b153ede.png" alt="crossover">
</div>
<p align="center"> Crossover <p> <br />

## Implemented idea of mutation



## Analysis and visualization



