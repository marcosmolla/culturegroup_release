### 1- Initialize parameters ...
#println("Here is your queue number: ", queue)
grid = collect(Base.product(
                # lines with value = 99 are currently depracted
[100],          # 1, number of individual nodes
[.75],          # 2, social learning success; #[.01 .25 .5 .75 .95]
99,             # 3, individual learning success; # [.01 .05 .1 .25 .5]
[100],          # 4, number of seed traits; # [2 5 10]
99,             # 5, number of tree levels; # [2 4 6]
[.1],           # 6, p_n; # [.5 .01 .05 .1 .2 .3 .4 .6 .7 .8]
[.01],          # 7, p_r; # [.001 0.005 .01 .03 .05 .07 .09 .2]
99,             # 8, trait out-degree; # [-2 1 2]
99,             # 9, number of rounds
99,             # 10, selection strength; # [.0 .2 .8]
0,              # 11 transmission, 0:full_horizontal, 1:full_vertical
1,              # 12 repetitions (local)
[4],# 1 0 6     # 13 payoff method [4 1 2 3 0]
[.01],# .1,     # 14 individual learning probability
[0],#1          # 15 density dependent payoffs (0: density independent, 1: density dependent)
[0], #-1 0 1    # 16 innovationSuccess, 0 contant, +1 increasing, -1 decreasing
0,#0.01         # 17 connection cost
false,#true     # 18 are pr and pn coupled
1, #0.01        # 19 mutation rate
2, # [2 6 10]   # 20 neighbourhood size, used for ring graphs only, requires 'evolveNetwork=false'
[1 2]))[queue]; # additional repetitions (global; to spread jobs over more cores)

# Record parameters for export
setup = Dict(
  # Setup for simulation
  "nRound"            => grid[1]*10, # The number here determines the number of generations the simulation will run for, i.e. one generation is equivalent to replacing all individuals (on average) once
  "repe"              => grid[12], # number of repetitions

  # Setup for graph and connection dynamics
  "nod"               => grid[1], # number of nodes
  "deg"               => 4, # average degree
  "pn"                => grid[6], #.8 # probaility to connect with neighbours of parents
  "pr"                => grid[7], #0.028 # probaility to connect with strangers

  # Setup for culture
  "nSeedTraitsTotal"  => grid[4], # total number of potential starting point for trait trees, no more than this can be innovated
  "startSeedTraits"   => grid[4], # number of seed traits an individual acquires at birth (assigned randomly at initialization)
  "socialLearningSuc" => grid[2],
  "indLearningSuc"    => grid[3], # overall success rate at changing cultural landscape (can be set for individual mechanisms later, see below)
  "selection"         => grid[10], # selection strength, with values >>0 strong selection, <<1 weak selection
  "transmission"      => grid[11], # [0...1], what's the proportion of learning turns that allow vertical transmission
  # "inverseTT"         => (grid[8]<0),
  # "traitoutDegree"    => abs(grid[8]),
  # "traitLevels"       => grid[5],
  # "maxMemory"         => 100, #25/200 low/unlimited

  # Setup for evolution
  "evolveSocLea"      => false, # Do values for social learning evolve?
  "evolveInnovation"  => false, # Do values for individual learning evolve?
  "evolvePN"          => true, # Do values for linking propensity for social inheritance evolve?
  "evolvePR"          => true, # Do values for linking propensity for random connections evolve?
  "evolveNetwork"     => false, # Does the social network evolve or is it fixed (e.g. a ring graph)?
  "evolveTransmission"=> false, # Do values for transmission evolve?
  "pnprcoupled"       => grid[18], # For simulations where PR and PN are coupled
  "payoffMethod"      => grid[13], # [0,1,4] for neutral, or repertoire and proficeincy based
  "indSuccessRate"    => grid[14], # Individual learning probability
  "densityDependence" => grid[15], # Density Dependence of Payoffs
  "innovationSuccess" => grid[16], # Success rate of Innocation (increasing, decreasing or constant)
  "cc"                => grid[17], # Connection cost
  "mutRate"           => grid[19], # Mutation rate for inherited parameters
  "neibhood"          => grid[20]  # Neighbourhood size
  );

# Make parameters globally accessable
nRound            = setup["nRound"];
repe              = setup["repe"];
nod               = setup["nod"];
deg               = setup["deg"];
pn                = setup["pn"];
pr                = setup["pr"];
nSeedTraitsTotal  = setup["nSeedTraitsTotal"];
startSeedTraits   = setup["startSeedTraits"];
socialLearningSuc = setup["socialLearningSuc"];
indLearningSuc    = setup["indLearningSuc"];
selection         = setup["selection"];
transmission      = setup["transmission"];
evolveSocLea      = setup["evolveSocLea"];
evolveInnovation  = setup["evolveInnovation"];
evolvePN          = setup["evolvePN"];
evolvePR          = setup["evolvePR"];
evolveNetwork     = setup["evolveNetwork"];
evolveTransmission= setup["evolveTransmission"];
payoffMethod      = setup["payoffMethod"];
indSuccessRate    = setup["indSuccessRate"];
densityDependence = setup["densityDependence"];
innovationSuccess = setup["innovationSuccess"];
cc                = setup["cc"];
pnprcoupled       = setup["pnprcoupled"];
mutRate           = setup["mutRate"];
neibhood          = setup["neibhood"];
# inverseTT         = setup["inverseTT"];
# traitoutDegree    = setup["traitoutDegree"];
# traitLevels       = setup["traitLevels"];
# maxMemory         = setup["maxMemory"];


## Initialize variables to store simulation states ...
# Individual information
  indPay          = zeros(nod); # recording an individual's payoff
  indSocLea       = setValues(evolveSocLea,nod,socialLearningSuc,.2); # recording an individual's SL propensity
  indInnovation   = setValues(evolveInnovation,nod,indLearningSuc,.2);# recording an individual's IL propensity
  indPN           = setValues(evolvePN,nod,pn,.2); # recording an individual's social inheritance parameter
  indPR           = setValues(evolvePR,nod,pr,.02); # recording an individual's random linking parameter
  indTransmission = setValues(evolveTransmission,nod,transmission,1); # recording an individual's transmission parameter
  # indTribe        = collect(1:nod); # recording tribe affiliation
  # indPrestige     = zeros(nod); # recording an individual's prestige

# Relationship matrix, to keep track of who is related to whom
  # relationship=zeros(nod, nod);

# Trait repertoir, matrix that holds each individual's proficiency/level along trait chains
  repertoir   = zeros(Int, nod, nSeedTraitsTotal);

# Initialize seed traits in individual repertoir
  for i in 1:nod
    repertoir[i, sample(collect(1:startSeedTraits))] = 1;
  end

# Tech-tree
 # Techtree is chains only in this version of the model, so there is no tech-tree initialised here

# # Population strucutre (as adjacency matrix)
#   adjm = initGraph(nod, deg, pn, pr);


# Initialise data structure for reporting summarised data for each repetition
  mmaxLevel = zeros(repe);
  mshared50 = zeros(repe);
  mshared90 = zeros(repe);
  mpay = zeros(repe);
  mntraits = zeros(repe);
  mbetadiv = zeros(repe);
  mpn = zeros(repe);
  mpr = zeros(repe);
  mtransmission = zeros(repe);
  mMedNTraits = zeros(repe);
  mMedMaxTraitLevel = zeros(repe);
  mTraitDiversity = zeros(repe);
  mdeg = zeros(repe);
  mpath = zeros(repe);
  mclust = zeros(repe);
  mclustLocal = zeros(repe);
  # mRepertoireSimiliarityA = zeros(repe);
  # mRepertoireSimiliarityR = zeros(repe);
  # mcomponentNumber = zeros(repe);
  # mcomponentMedSize = zeros(repe);
  # mcomponentMaxSize = zeros(repe);
  # mlearnedsoc = zeros(repe);
  # mlearnedino = zeros(repe);

# Initialise data structure for summarising results in R
  R"recPNl <- recPRl <- recDegl <- recClustl <- recPathl <- recTraitNl <- recTraitLl <- recNTraitsl <- recBetaDivl <- recLearnedSocl <- recLearnedInol <- recPayl <- recTraitDiversityl <- adjml <- list()"

# Initialise often used data objects
  allInds         = 1:nod; # range of all individuals
  availableSeeds  = nSeedTraitsTotal; # fixed, but could be growing in a future itteration, but then this needs to be inside the for loop
  allSeeds        = collect(1:availableSeeds); # as availableSeeds are fixed so far, I can also get the vector here
  pM              = payoffMethod;
  vT              = round(transmission*100); # number of rounds with vertical transmission

# Set Payoff Calculation Method
## In the main text we only report results for:
## neutral selection (0),
## selection based on repertoire size (1), and
## selection based on proficiency level (4)
  if payoffMethod == 0 # Neutral model (for densityDependence 0 AND 1)
    function payM(newRep,rep,n)
      return 0
    end
  elseif payoffMethod == 1 && densityDependence == 0 # Number of traits known
    function payM(newRep,rep,n)
      return sum(newRep.>0)
    end
  elseif payoffMethod == 1 && densityDependence == 1
    function payM(newRep,rep,n) #!!!!!
      if sum(newRep)==0
        return 0
      else
        return sum(1.-sum(rep.>0,1)[newRep.>0]/n) / sum(newRep.>0) # sum of average novelty
      end
    end
  # elseif payoffMethod == 2 # Total level count
  #   function payM(newRep,rep,n)
  #     return mean(1.-(freq[newRep.>0]./nod))*sum(newRep.>0) # novelty * sum of all levels #sum( newRep[newRep.>0] ./ freq[newRep.>0 ])
  #   end
  # elseif payoffMethod == 3 # Numer of traits * Total level count
  #   function payM(newRep,rep,n)
  #     return sum(1.-(freq[newRep.>0]./nod)) * (mean(1.-(freq[newRep.>0]./nod))*sum(newRep.>0)) # sum(1./(freq[newRep.>0])) * sum( newRep[newRep.>0] ./ freq[newRep.>0 ])
  #   end
elseif payoffMethod == 4 && densityDependence == 0 # Highest level
    function payM(newRep,rep,n)
      return maximum(newRep) # Only the ONE highest
      # return sum(sort(newRep,rev=true)[1:5]) # The highest FIVE
    end
  elseif payoffMethod == 4 && densityDependence == 1 #!!!!!
    function payM(newRep,rep,n)
      if sum(newRep)==0
        return 0
      else
        maxr=maximum(newRep);
        t=sample(find(newRep.==maxr));
        return maxr * ( 1.-sum(rep[:,t].>0)/n ) # highest level * trait novelty
      end
    end
  # elseif payoffMethod == 5 && densityDependence == 0
  #   function payM(newRep,rep,n)
  #     return sum(newRep) # average novelty  #sum(1./(freq[newRep.>0]))
  #   end
  # elseif payoffMethod == 5 && densityDependence == 1
  #   function payM(newRep,rep,n)
  #     if sum(newRep)==0
  #       return 0
  #     else
  #       t=find(newRep.>0);
  #       s=zeros(length(t));
  #       i=1;
  #       for T in t
  #         s[i] = newRep[T] / mean(rep[rep[:,T].>0,T]); # focal trait level relative to trait level of all with the trait
  #         i+=1;
  #       end
  #       return sum(s) / sum(newRep.>0)
  #     end
  #   end
  # elseif payoffMethod == 6 && densityDependence == 0
  #   function payM(newRep,rep,n)
  #     total=0;
  #     for i in newRep
  #         total+=sum(1:i)
  #     end
  #     return total # here, leveling up is worth more than just increasing number of traits
  #   end
  # elseif payoffMethod == 7 && densityDependence == 0 # Number of traits known with a logistic payoff function
  #   function payM(newRep,rep,n)
  #     #return 10/(1+exp(-.5*( sum(newRep.>0) -5))) # L/(1+exp(-k*(x-x0)))
  #     return (1+.33)^(sum(newRep.>0))
  #   end
  else
    error("PayoffMethod ",payoffMethod," is not defined.")
  end

# Set correct innovation success probability calculation (depending on proficeincy level and whether innovation success is increasing (1), decreasing(-1), or constant (0))
  if innovationSuccess == 0#"constant"
    function probI(ISR, CLI)
      ISR
    end
  elseif innovationSuccess == -1#"decreasing"
    function probI(ISR, CLI)
      ISR/(CLI+1)
    end
  elseif innovationSuccess == 1#"increasing"
    function probI(ISR, CLI)
      ISR*(CLI+1)
    end
  end
