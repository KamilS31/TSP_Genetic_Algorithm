# TSP_Genetic_Algorithm

## Issue to be solved

The travelling salesman problem asks the following question: "Given a list of cities and the distances between each pair of cities, what is the shortest possible route that visits each city exactly once and returns to the origin city?" In other words, the TSP seeks to find a Hamiltonian cycle in an undirected complete graph with weighted edges, where the goal is to minimize the sum of the edge weights along the cycle.

In this version of the problem, I will try to find the shortest route between points (cities) plotted on a simplified map of the United States.

*Note, we assume that we are moving the aircraft at a constant altitude and do not take into account the curvature of the earth.*

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

*_It is worth mentioning that the probability of crossing individuals is calculated using the rank selection._


## Implemented idea of mutation

Mutation of individuals with probability 1/chromosome_length proceeds as follows:
1) We draw two different locus (two different numbers from 1 to the total number of cities) from the chromosome.
2) Genes under the drawn locus are swapped with each other to form a mutation.

The image below visualizes the whole process.

<div align="center">
  <img src="https://user-images.githubusercontent.com/127042515/234563837-af99553b-c429-4c6a-97df-07159b2b9a79.png" alt="mutation">
</div>
<p align="center"> Mutation <p> <br />

## Important parameters

Variables affecting the speed of operation and the quality of the result:
* nCities - The number of cities to visit, or in other words, the length of a chromosome.
* nPopulation - Size of the population.
* lack_imp - Number of populations generated without improving the best individual.

## Analysis and visualization

**The result of the algorithm for the number of cities equal to 10:**

* Percentage of searching the state space: 0.30864%
* Shortest found route: 709.5766
* An 35.0135% improvement over the best and the worst
* Elapsed time is 3.643849 seconds
<div align="center">
  <img src="https://user-images.githubusercontent.com/127042515/234570630-4f9237d6-b2d9-47c5-813c-a67d3fccbe48.png" alt="10_points">
</div>
<p align="center"> Route (nCities = 10) <p> <br />

<div align="center">
  <img src="https://user-images.githubusercontent.com/127042515/234571114-2a8cb4ea-f85e-449a-823a-462e37abfd61.png" alt="10_points_plot">
</div>
<p align="center"> Sorted route length against number of generations (nCities = 10) <p> <br />

Settings:
* nPopulation = 30
* lack_imp = 100
* counter == 50


**The result of the algorithm for the number of cities equal to 25:**

* Percentage of searching the state space: 1.1347e-19%
* Shortest found route: 1163.5745
* An 51.9150% improvement over the best and the worst
* Elapsed time is 9.535783 seconds.

<div align="center">
  <img src="https://user-images.githubusercontent.com/127042515/234576565-511e2f2f-6be2-4a96-83c6-73e4ee7bd953.png" alt="25_points">
</div>
<p align="center"> Route (nCities = 25) <p> <br />

<div align="center">
  <img src="https://user-images.githubusercontent.com/127042515/234576869-f3df0c86-ce28-4392-9018-44f1e5736b70.png" alt="25_points_plot">
</div>
<p align="center"> Sorted route length against number of generations (nCities = 25) <p> <br />

Settings:
* nPopulation = 300
* lack_imp = 300
* counter == 100


**The result of the algorithm for the number of cities equal to 50:**

* Percentage of searching the state space: 1.1771e-58%
* Shortest found route: 1606.5968
* An 64.5187% improvement over the best and the worst
* Elapsed time is 1606.5968 seconds.

<div align="center">
  <img src="https://user-images.githubusercontent.com/127042515/234582487-72f4c278-f3dd-4f26-82f4-415f7215d751.png" alt="50_points">
</div>
<p align="center"> Route (nCities = 50) <p> <br />

<div align="center">
  <img src="https://user-images.githubusercontent.com/127042515/234582723-488481de-daed-44cd-a7da-9bedf1b1711b.png" alt="50_points_plot">
</div>
<p align="center"> Sorted route length against number of generations (nCities = 50) <p> <br />

Settings:
* nPopulation = 5000
* lack_imp = 700
* counter == 100

### Visualizing the evolution of solutions

The graph below shows the value of the objective function for the best individual from each population. We can see a clear downward trend on it, which proves the evolution of solutions in the direction of the most optimal route.

<div align="center">
  <img src="https://user-images.githubusercontent.com/127042515/234588426-f19f5193-fd7b-4da7-a50e-7d242713d50f.png" alt="plot">
</div>
<p align="center"> Route length against number of generations (nCities = 30) <p> <br />


