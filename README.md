# Running CultureGroup simulation

## File overview
The following files are found in the `methods/` folder
* `sbatch.sh` initates a job array on slurm
* `run.jl` initates julia simulations after loading methods from `methods.jl` and the job array specific parameters from `main.jl`
* `summariseFiles.R` opens all `*.Rdata` files from the output directory, summarises values, returns a summary file into the working directory, zips the output folder, deletes the output folder and finally deletes all `slurm*.out` files

## Requirements
Running the simulation code requires:
* `julia` (dependencies: `Distributions`, `StatsBase`, `RCall`)
* `r` (dependencies: igraph, dplyr, reshape2)

## Simulations reported in our article
In our article, *Cultural Selection Shapes Network Structure* (https://www.biorxiv.org/content/early/2018/11/08/464883), we report results for the following iterations of the model:

1. Social learning dynamics for **fixed simple graphs**
2. Social learning dynamics for dynamic complex networks with **fixed linking parameters**
3. Social learning dynamics for dynamic complex networks with **evolving linking parameters**
4. Social learning dynamics for complex networks with evolving but **coupled linking parameters**
5. Social learning dynamics for dynamic complex networks with **switching selection regimes**

Furthermore, in the ESM we report results for simulations with

6. Low mutation rate
7. Connection costs
8. Varying population size and trait number
9. Varying innovation and social learning success rate

To run the individual simulations follow the steps outlined below to adjust the simulation.

### Adjust simulation
#### 1. Social learning dynamics for fixed simple graphs
To use a simple graph such as a ring (as in the main text), find `# Setup for evolution` in the `main.jl` file and make sure `evolveNetwork` is set to `true`. This will which will create a ring graph of size `nod` with neighbourhood `neibhood`, and keep the network shape fixed. Note: this is currently _only implemented for neutral selection_, as parents for a newborn have to be adjacent to the newborn and are not selected based on their fitness.

Can be combined with: 8, and 9

#### 2. Social learning dynamics for dynamic complex networks with fixed linking parameters
To let linking parameters evolve, find `# Setup for evolution` in the `main.jl` file and make sure `evolvePN` and `evolvePR` are set to `false`. This will keep the linking parameters fixed throughout a simulation, while the network is still dynamically rewired.

Can be combined with: 5, 7, 8, and 9

#### 3. Social learning dynamics for dynamic complex networks with evolving linking parameters
To let linking parameters evolve, find `# Setup for evolution` in the `main.jl` file and make sure `evolvePN` and `evolvePR` are set to `true`. This will let the parameters evolve (and mutate) throughout a simulation.

Can be combined with: 6, 7, 8, and 9

#### 4. Social learning dynamics for complex networks with evolving but coupled linking parameters
Change `grid` parameter in line 18 (`pnprcoupled`) from `false`(not coupled), to `true` (couples `pr` to `pn` given an average degree `k` and population size `N`). In the main text we used average degrees 2, 6, 10, for a population of 100 individuals.

Can be combined with: 6, 7 (note, in this case mutation only affects PN), 8, and 9

#### 5. Social learning dynamics for dynamic complex networks with switching selection regimes
To enable switching payoff regimes find `## 3EXPLOIT` in `methodsJL` and uncomment the section starting with `# Switching payoff method twice throughout the simulation`. This will change the payoff regime twice throughout a simulation run (after 1/3 and 2/3 of the rounds).

Can be combined with: 6, 7, 8, and 9

#### 6. Low mutation rate
Change `grid` parameter in line 19 (mutation rate, `mutRate`). We used `1` for the main text and `0.01` for the appropriate simulations in the ESM.

#### 7. Connection costs
Change `grid` parameter in line 17 (connection cost, `cc`). We used `0` for the _no cost_ and `0.01` for the _cost condition_.

#### 8. Varying population size and trait number
Change `grid` parameter in line 1 (population size) and line 4 (number of seed traits). In the main text we used `100` for each parameter.

#### 9. Varying innovation and social learning success rate
Change `grid` parameters in line 2 (social learning success rate, `socialLearningSuc`) and line 14 (individual learning success rate, `indSuccessRate`). In the main text we used `0.75` and `0.01` respectively.
