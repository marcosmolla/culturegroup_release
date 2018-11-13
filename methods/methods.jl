### Declare functions
# A function to bound the upper and lower limits of a vector (or single value)
  function bound(x,lower=0,upper=1)
    return min(upper, max(lower, x))
  end

# A function to set initial values for ind arrays
  function setValues(evolve::Bool,times,mean,sd)
    if evolve
      return bound.(rand(Normal(mean,sd), times))
    else
      return bound.(fill(mean, times))
    end
  end

# Dynamic network updating, takes adj.matrix (ADJM), probabilities to connecto to parent(s), neighbours, and randoms (PB, PN, PR), the index of the parent (if NULL, default, NPARENT number of individuals are randomly choosen as parent), number of parents (NPARENT, default is 1)
  function pbpnprDyn(ADJM, PN, PR, NOD, DB=0)
      # one iteration of the basic model
      if DB==0
        DB = sample(1:NOD,2, replace=false); #sample a death and a birth, without replacement; 1 newborn, 2 parent
      end
      inherit = ADJM[:,DB[2]].*(1.-(rand(NOD)-PN .> 0.)); # socially inherited connections
      randconn = (1-ADJM[:,DB[2]]).*(1.-(rand(NOD)-PR .> 0.)); # random connections
      newconn = inherit+randconn #total connections
      newconn[DB[2]]=1; #connect to the parent
      ADJM[:,DB[1]]=newconn; #replace the dead individual with the newborn
      ADJM[DB[1],:]=newconn; #same
      ADJM[DB[1],DB[1]]=0; #set self-link equal to zero
      return ADJM
  end

# Function to calculate PR when coupled to PN in relatin to population size N and aimed for average degree K
  function PR(K,N,PN)
    (K*(N-1-PN*(N-2))-N+1) / ((N-2)*(N-1-K))
  end

# Select lower triangle of a matrix (excluding the diagonal)
  function lowerTri(x)
    res = zeros(Int,size(x)[1],size(x)[1]);
    for i in 2:size(x)[1]
      res[i,1:(i-1)] = 1;
    end
    return Bool.(res)
  end

# Random Graph Generation, takes number of nodes (NOD), average degree (DEG)
  function initGraph(NOD, DEG, PN, PR)
    P = (DEG*NOD)/(NOD^2); # connection probability
    ADJM = zeros(Int, NOD, NOD); # Set adjacency matrix
    diag = lowerTri(ADJM); # Select one half of the matrix
    npairs = sum(diag); # Number of pairs without diagonal (using only one matrix traingle)
    ADJM[diag] = rand(npairs).<= P;
    ADJM = ADJM + transpose(ADJM); # Populate other half of the matrix
    for i in 1:(10*NOD) # Let PBPNPR Dynamics run for burning in period
      ADJM = pbpnprDyn(ADJM, PN, PR, NOD);
    end
    return ADJM # return adjacency matrix
  end

# Ring Graph Generation, takes number of nodes (NOD), and neighbourhood size (NEIBHOOD)
  function initRING(NOD, NEIBHOOD)
    @rput NOD NEIBHOOD
    R"""
    ## RING
    net <- graph.ring(n=NOD)
    cat(paste("Ring with",nod,"vertices\n"))
    net <- connect.neighborhood(graph=net, order=NEIBHOOD)

    ## SMALL WORLD
    # net <- watts.strogatz.game(dim=1, size=setup$nod, nei=2, p=setup$neibs/setup$nod, loops=F)
    # adjm <- adjFromEdge(get.edgelist(as.directed(net)))

    ## Strong/Weak connectecd components (maximal connected subgraph), e.g. 5 components (á 5 nodes), could be connected 1 with another (circle, or random), or each component with two others, etc till each is connected with each other comonponent (5)
    # net <- graph_components(conn=setup$neibs)

    ## For Simulations with fixed degree (4) and varliable clustering coefficient
    # net <- graph_fixedDegVariableClust(goal_trans=setup$neibs, n=setup$nod, m=ifelse(setup$nod==25, 50,500), iterations=10^5)

    # Return adjacency matrix
    adjm<-as.matrix(get.adjacency(net))
    """
    @rget adjm
    adjm=round.(Int, adjm);
  end

# A function to inherit parental values, or reuse the standard value if this trait does not evovle
  function inherit(;evolve=evolve, parents=parents, standard=standard, sd=0.01, rmutation=0.01)
    if evolve
      if rand()<=rmutation
        return bound.(rand(Normal(parents, sd))) # bound.(rand(Normal(mean(parents), sd))) #(use this for a sexual reproduction version of the model)
      else
        return parents #mean(parents) #(use this for a sexual reproduction version of the model)
      end
    else
      return standard
    end
  end

# This function is similar to pbpnprDyn, howerer, it uses only one parameter. (currently not used)
# function paDyn(ADJM, PA, NOD, DB=0)
#     if DB==0
#       DB = sample(1:NOD,2, replace=false); #sample a death and a birth, without replacement; 1 newborn, 2 parent
#     end
#     newconn = 1.-(rand(NOD)-PA .>0); #total connections
#     newconn[DB[2]]=1; #connect to the parent
#     ADJM[:,DB[1]]=newconn; #replace the dead individual with the newborn
#     ADJM[DB[1],:]=newconn; #same
#     ADJM[DB[1],DB[1]]=0; #set self-link equal to zero
#     return ADJM
# end

# Function to determine which traits WHO can readily learn based on what they already know (this is required for branching tech-trees, however, here we only use chains and so this function is currently not used)
  # function selectNewTrait(;NNK=NNK, RP=RP, SEEDS=SEEDS, TTM=TTM, OK=OK)
  #   nt = Int[];
  #   # Prepare while loop
  #   h = 1;
  #   # Go through all the traits from the neighbours based on their payoff, if it is a seed trait learn it, if not learn it only if previous knowledge allows direct link
  #   maxItr = length(NNK);
  #   while h < maxItr
  #     # select highest paying
  #     highest = NNK[RP][h];
  #     # root node?
  #     if highest in SEEDS
  #       # learn root node/seed trait
  #       nt = highest;
  #       # and exit while loop
  #       h = maxItr+1;
  #     else
  #       # Check whether this trait can be connected to traits that are in the individuals repertoir
  #       # OR if the new trait requires two inputs, check that the individual already possesses both of those inputs
  #       if all(TTM[highest.==TTM[:,2],1] in OK) # this checks whether all parent nodes needed for the new trait are known to the individual
  #         # learn this trait
  #         nt = highest;
  #         # and exit while loop
  #         h = maxItr+1;
  #       end
  #     end
  #     h = h+1
  #   end #// END while loop
  #
  #   return nt
  # end

# Function to convert an edge list into an adjacency matrix and back
  # function edgeFromAdj(techMatrix=Array{Int64,2})
  #   from = Int[];
  #   for i in 1:size(techMatrix)[1]
  #     from = vcat( from, repeat([i], inner=sum(techMatrix[:,i])) );
  #   end
  #
  #   to = Int[];
  #   for i in 1:size(techMatrix)[1]
  #     to = vcat(to, find(techMatrix[i,:].==1));
  #   end
  #
  #   return hcat(from, to)
  # end

# Function to calculate the similarity between individual repertoires (currently not used)
  # function repSim(X, NEWBORN, REPBOO, REPSIM)
  #   for j in X
  #     denom = sum(REPBOO[NEWBORN,:].|REPBOO[j,:]); # in the rare case that two individuals have knowledge at all
  #     if denom==0
  #       tmp = 0;
  #     else
  #       tmp = sum(REPBOO[NEWBORN,:] .& REPBOO[j,:]) / denom # 1 means identical, 0 means no similarity at all
  #     end
  #     REPSIM[NEWBORN,j]=REPSIM[j, NEWBORN]=tmp;
  #   end
  #   return REPSIM
  # end

# Returns column and row number, given the 1D id for a position in a matrix, and the number of columns of that matrix (currently not used)
  # function matrixRowCol(POS,NSEEDS)
  #   remainder = rem(POS,NSEEDS);
  #   row = div(POS,NSEEDS) + (remainder>0);
  #   col = remainder + NSEEDS*(remainder==0);
  #   return row, col
  # end


#@@@  @@@@@@@@@@@   @@@@@@@   @@@@@@@ @@@@@@@@    @@@@@@@
#@@@   @@@@@@@@@@@   @@@@@     @@@@@   @@@@@@@@    @@@@@
#@@@   @@@@   @@@@@  @@@@@     @@@@@   @@@@@@@@@   @@@@@
#@@@   @@@@   @@@@   @@@@@     @@@@@   @@@@@@@@@@  @@@@@
#@@@   @@@@@@@@@@    @@@@@     @@@@@   @@@@@ @@@@@ @@@@@
#@@@   @@@@@@@@@@@   @@@@@     @@@@@   @@@@@ @@@@@ @@@@@
#@@@   @@@@   @@@@   @@@@@     @@@@@   @@@@@  @@@@@@@@@@
#@@@   @@@@    @@@@  @@@@@     @@@@@   @@@@@   @@@@@@@@@
#@@@   @@@@    @@@@@  @@@@@@ @@@@@@    @@@@@    @@@@@@@@
#@@@  @@@@@@   @@@@@    @@@@@@@@@     @@@@@@@    @@@@@@@@

### RUN MODEL
function runit(;adjm=nothing)
  print("Starting simulation updatedSL .")

  # Run several repertitions for the same set of parameters
  for reps in 1:repe
    # For each repetition, reset world:
      # Population strucutre (as adjacency matrix)
        if evolveNetwork
          adjm = initGraph(nod, deg, pn, pr);
        else
          adjm = initRING(nod, neibhood);
        end

      # Individual information
        indTribe        = collect(1:nod);
        indPay          = zeros(nod);
        indPrestige     = zeros(nod);
        indSocLea       = setValues(evolveSocLea,nod,socialLearningSuc,.2);
        indInnovation   = setValues(evolveInnovation,nod,indLearningSuc,.2);
        indPN           = setValues(evolvePN,nod,pn,.2);
        indPR           = setValues(evolvePR,nod,pr,.02);
        indTransmission = setValues(evolveTransmission,nod,transmission,1);
        #relationship = zeros(nod, nod); # Relationship matrix, to keep track of who is related to whom (currently not used)

      # Individual repertoires (plus adding informatino to initial repertoire)
        repertoir   = zeros(Int, nod, nSeedTraitsTotal); # Trait repertoir, matrix that holds each individual's proficiency/level along trait chains
        for i in 1:nod
          repertoir[i, sample(collect(1:startSeedTraits))] = 1; # Initialize seed traits in individual repertoir
        end

      # Objects to store repertoire similarity
        repSimA=zeros(nod,nod);
        repSimR=zeros(nod,nod);

      # Recording data, once per generation
        recgen=1; # counter for how far into the generation we are, at recGen==nod record and reset
        recid = 1; # index for recording data
        recLength = (nRound/nod); #
        recPay = zeros(recLength);
        recNTraits = zeros(recLength);
        recBetaDiv = zeros(recLength);
        recPN = zeros(recLength);
        recPR = zeros(recLength);
        recTransmission = zeros(recLength);
        recDegree = zeros(recLength);
        recPath = zeros(recLength);
        recClustLocal = zeros(recLength);
        recRepertoire = zeros(recLength, setup["nSeedTraitsTotal"]);
        recMedNTraits = zeros(recLength);
        recMedMaxTraitLevel = zeros(recLength);
        recTraitDiversity = zeros(recLength); # Trait diversity within the populatin
        # recRepertoireSimiliarityA = zeros(recLength); # Repertoire similiarity between each diad
        # recRepertoireSimiliarityR = zeros(recLength); # Repertoire similiarity between relatives
        #currently not used:
        #recIndLea = zeros(recLength);
        #recSocLea = zeros(recLength);
        #recLearnedSoc = zeros(recLength);
        #recLearnedIno = zeros(recLength);
        #recTraitLevels = zeros(recLength);
        #recTribe = zeros(recLength);

        recRep = zeros(recLength,nSeedTraitsTotal);

  # Start single repetition of a simulation
    for times in 1:nRound
      ## 1 ////////////////////////
      ## DEATH-BIRTH PROCESS, MORAN
      # We start be removing an individual relative to the inverse of its fitness.
      # A new individual inherits traits for inividual and social learning, as well as
      # pn and pr values, and seed traits and trait proficiency from its parent.

      ## Select dead/newBorn individual
        newBorn = sample(allInds); # ALWAYS randomly
        # newBorn = sample(allInds[indPay.==minimum(indPay)]); # Select based on lowest payoff

      ## Select parent, choose parent from alive, relative to indPay, if any indPay are != 0
        if evolveNetwork
          if any((indPay.>0)[1:end .!= newBorn])
            # If some individuals have non zero payoff sample with wieghts ...
            parents = wsample(allInds[1:end .!= newBorn], indPay[1:end .!= newBorn]);
          else
            # ... otherwise choose randomly
            parents = sample(allInds[1:end .!= newBorn]);
          end
        else
          # For the case were networks remain stable (note, this is currently only implemented for neutral selection)
          parents = sample(find(adjm[:,newBorn]));
        end

      ## Inherit attributes
        indPN[newBorn]        = inherit(evolve=evolvePN, parents=indPN[parents], standard=pn, sd=0.1, rmutation=mutRate); # sd=0.01, rmutation=.01
        if pnprcoupled # PR is either adjusted based on PN (when coupled) or directly inherited with mutation from parent (when not coupled)
          indPR[newBorn]        = bound(PR(6,nod,indPN[newBorn])); # couple PR to PN
        else
          indPR[newBorn]        = inherit(evolve=evolvePR, parents=indPR[parents], standard=pr, sd=0.01, rmutation=mutRate); # PR is PN-independent
        end
          #indSocLea[newBorn]    = round.(inherit(evolve=evolveSocLea, parents=indSocLea[parents], standard=socialLearningSuc, sd=0.1, rmutation=.05),2);
          #indInnovation[newBorn]= round.(inherit(evolve=evolveInnovation, parents=indInnovation[parents], standard=indLearningSuc, sd=0.2, rmutation=.05),2);
        # indTransmission[newBorn]= round(Int,inherit(evolve=evolveTransmission, parents=indTransmission[parents], standard=transmission, sd=10, rmutation=0.01));

      ## Reset values
        # Reset memory
          repertoir[newBorn,:] = 0;
        # Reset payoff
          indPay[newBorn] = 0;
        # Connect new individual (position of the newBorn) to the network based on PBPNPR dynamics (onle when networks are dynamic)
          if evolveNetwork
            adjm = pbpnprDyn(adjm, indPN[newBorn], indPR[newBorn], setup["nod"], [newBorn, parents]); # for pbpnpr prob
          end
        # Reset tribe (maybe also move to new tribe if random is larger than PN) (currently not used)
          #indTribe[newBorn] = indTribe[parents];
        # # Update relationship (currently not used)
          # relationship[parents,newBorn] = 0; # reset newborn parent first, otherwise there is a carry over problem in the next step
          # relationship[newBorn,:] = relationship[:,newBorn] = round.(relationship[parents,:]/2,3); # Inherit parent relationship, divided by half
          # relationship[newBorn,parents] = relationship[parents,newBorn] = 1; # Set parent offspring relation == 1
        # # Reset repertoire similarity for newBorn and its connected individuals (currently not used)
          # repSimA[newBorn,:] = repSimA[:,newBorn] = 0;
          # repSimR[newBorn,:] = repSimR[:,newBorn] = 0;


      ## 2 /////
      ## EXPLORE
      # Subsequently the individual is performing several learning turns, in which it has the opportunity
      # to innovate or to socially learn from its connections.

      ## PREPARE SOCIAL LEARNING
      # Determine repertoire contents of neighbours
        # Neigbhours
        neibn = find(adjm[:, newBorn]); # neighbour ID
        neibs = length(neibn); # number of neighbours
        neibRep = repertoir[neibn,:]; # repertoire in neighbourhood
        neibRepBool = neibRep.!=0; # boolean repertoire of neighbourhood
        # Parent(s), required for vertical transmission
        pW = repertoir[parents,:]; # parents weights for traits
        pW[pW.!=0]=1; # use this if it is not based on the parents proficiency, skip it if you want to keep proficiency in there
        # Calculate frequency of performance
        tmp=sum(neibRepBool,1)./sum(neibRepBool); tmp[tmp.==0]=0;
        weightsSoc = Weights(tmp[:]); # weights based on frequency of performance

      for LT in 1:100
        # Social learning
          socTrait = 0; # reset social trait
          if (neibs!=0)&(any(neibRepBool)) # If there are any traits known to the neighbours, proceed
          ## Biased sampling based on frequency of performance among neighbours
            if (vT!=0 )& (vT<=LT)
              socTrait = sample(allSeeds, pW); # vertical transmission
            else
              socTrait = sample(allSeeds, weightsSoc); # vertical and oblique transmission
            end

            if socTrait!=0 # if there are any traits that can be learned socially
              if any(repertoir[newBorn,socTrait].<repertoir[neibn,socTrait]) # if any neighbour is more proficient than newborn ...
                if rand() <= (indSocLea[newBorn]*tmp[socTrait]) # soc success is based on proportion of time observing someone with a higher level
                  repertoir[newBorn,socTrait] += 1; # if successful, increase repertoire level for the trait by one
                end
              end
            end
          end

        # Individual learning
          inoTrait = sample(allSeeds); # unbiased sampling; maybe add Weights() to mimic an guided/directed innovation process
          if rand() <= probI(indSuccessRate, repertoir[newBorn,inoTrait])
            repertoir[newBorn,inoTrait] += 1; # if innovation is successful, increase repertoire level for the trait by one
          end
      end


      ## 3EXPLOIT
      # Now the newBorns lifetime foraging success is determined based on what it has learned before and how proficient it is with its traits
      # Alt1: standard payoff
        # indPay[newBorn] = payM(repertoir[newBorn,:], repertoir, nod); # no connection costs
        indPay[newBorn] = max(0, payM(repertoir[newBorn,:], repertoir, nod) - (sum(adjm[:,newBorn])*cc)); # with connection costs
      # Alt2A: payoff with variable selection strenght
        # # d=0.01; # weak selection; looks like this is too weak
        # d= 0.2;#0.1; # strong selection
        # indPay[newBorn] = payM(repertoir[newBorn,:], repertoir, nod); #indPay[newBorn] = (1+d)^payM(repertoir[newBorn,:], repertoir, nod);
      # Alt2B: payoff with variable selection strength AND degree cost (both for pn and pr)
        # d= 0.2;
        # cc=0.1; # connection cost (penalty per degree)
        # indPay[newBorn] = payM(repertoir[newBorn,:], repertoir, nod) - (sum(adjm[:,newBorn])*CC);
      # Alt3: penalty for having connections (both for pn and pr)
        # indPay[newBorn] = max(0,payM(repertoir[newBorn,:], repertoir, nod)-sum(adjm[:,newBorn])/(nod-1));

      # Switching payoff method twice throughout the simulation
        # if (times == floor(Int, nRound*1/3)) | (times == floor(Int,nRound*2/3))
        #   pM==1 ? pM=4 : pM=1; # (1=repertoire, 4=proficiency)
        # end
        # if pM == 4
        #   indPay[newBorn] = max(0,(maximum(repertoir[newBorn,:])));# - (sum(adjm[:,newBorn])*cc));
        # else
        #   indPay[newBorn] = max(0,(sum(repertoir[newBorn,:].>0)));# - (sum(adjm[:,newBorn])*cc));
        # end

      # // END EXPLOIT

      ## 4RECORD
      ## Record Results
      if recgen==nod
        repBoo = repertoir.≠0; # turn repertoire into boolean
        recPay[recid]         = mean(indPay);
        recNTraits[recid]     = count(sum(repertoir,1).>0); # number of all known traits in the population
        recBetaDiv[recid]     = sum(sum(repBoo,1).!=0) / mean(sum(repBoo,2));
        recPN[recid]          = mean(indPN);
        recPR[recid]          = mean(indPR);
        recTransmission[recid]= mean(indTransmission);
        recDegree[recid]      = sum(adjm)/nod;
        recMedNTraits[recid]  = mean(sum(repBoo,2)); # median number of known traits per individual
        recMedMaxTraitLevel[recid] = mean(maximum(repertoir, 2)); # median of highest trait level for each individual
        # recIndLea[recid]      = mean(indInnovation);
        # recSocLea[recid]      = mean(indSocLea);
        # recTraitLevels[recid] = maximum(repertoir); # highest overall trait level
        # recTribe[recid]       = length(unique(indTribe));
        # recRepertoire[recid,:] = mean(repertoir, 1)[:];

        @rput adjm nod
        R"""
        net <- graph_from_adjacency_matrix(adjm)
        cl <- transitivity(net, type='local')
         clustLocal <- sum(cl[!is.nan(cl)])/nod # not mean, because that would disregard the nodes without neighbours
        path <- average.path.length(net)
        """
        @rget path clustLocal

        recPath[recid] = mean(path);
        recClustLocal[recid] = mean(clustLocal);

        # Caluclate repertoir similiarity using Jaccard index and update repSim matrices
        # repSimA = repSim(find(adjm[newBorn,:]),newBorn,repBoo,repSimA);
        # repSimR = repSim(find(relationship[newBorn,:]),newBorn,repBoo,repSimR);
        # recRepertoireSimiliarityA[recid] = mean( repSimA[find(adjm)] ); # ... amongst neighbours
        # recRepertoireSimiliarityR[recid] = mean( repSimR[find(relationship)] ); # ... amongst relatives

        #> Calculate diversity index for repertoir
        # # For each individual
        # relative = repertoir./ (sum(repertoir,2)); #(p_i / p_total)
        # # for Simpson index
        #   x=(relative).^2;
        #   Hsimp=1-sum(x,2);
        # # for Shannon
        #   x=-(relative).*log(relative);
        #   Hshan=sum(x,2);
        # For the entire population
        repBooSum=sum(repBoo, 1);
        recTraitDiversity[recid]=1-sum( (repBooSum/sum(repBooSum)).^2); # Simpson's index # 1 means highly diverse, 0 means not diverse at all
        # recTraitDiversity=sum(-(x/sum(x)).*log(x/sum(x))); # Shannon
        #< END Calculate diversity index for repertoir

        recRep[recid,:] = sum(repBoo,1); # returns the number of traits known to the population

        # Recording several network snapshots per run
          # if (times == floor(Int, nRound*.165)) | (times == floor(Int,nRound*.5)) | (times == floor(Int,nRound*.66))
          #   @rput adjm
          #   R"adjml <- c(adjml, adjm)"
          # end

        recid = recid+1;
        recgen = 1;

        # println(gen);
        # gen+=1;
      else
        recgen += 1;
      end
    end # END SIMULATION FOR-LOOP


  # SUMMARISING RESULTS OF A SINGLE SIMULATION RUN
    # At the end of a simulation run summarise values
    subs=round(Int, (recLength-recLength*.2)):Int(recLength); # keep final 20% of rounds
    repCols = sum(repertoir.>0,1); # ,1] is colsums ,2] is rowsums
    mmaxLevel[reps] = median(maximum(repertoir,2));
    mshared50[reps] = sum(repCols.>=(.5*nod));
    mshared90[reps] = sum(repCols.>=(.9*nod));
    mpay[reps]=mean(recPay[subs]);
    mntraits[reps]=mean(recNTraits[subs]);
    mbetadiv[reps]=mean(recBetaDiv[subs]);
    mpn[reps]=mean(recPN[subs]);
    mpr[reps]=mean(recPR[subs]);
    mtransmission[reps]=mean(recTransmission[subs]);
    mMedNTraits[reps]=mean(recMedNTraits[subs]);
    mMedMaxTraitLevel[reps]=mean(recMedMaxTraitLevel[subs]);
    mTraitDiversity[reps]=mean(recTraitDiversity[subs]);
    mdeg[reps]=mean(recDegree[subs]);
    mpath[reps]=mean(recPath[subs]);
    mclustLocal[reps]=mean(recClustLocal);
    # mlearnedino[reps]=mean(recLearnedIno[subs]);
    # mlearnedsoc[reps]=mean(recLearnedSoc[subs]);
    # mRepertoireSimiliarityA[reps]=mean(recRepertoireSimiliarityA[subs]);
    # mRepertoireSimiliarityR[reps]=mean(recRepertoireSimiliarityR[subs]);
    # mcomponentNumber[reps]=componentNumber;
    # mcomponentMedSize[reps]=componentMedSize;
    # mcomponentMaxSize[reps]=componentMaxSize;

  # Store data in an R object
    @rput recPN recPR recDegree recClustLocal recPath recMedNTraits recMedMaxTraitLevel recNTraits recBetaDiv recPay reps setup adjm repertoir recRep recTraitDiversity; #recLearnedSoc recLearnedIno relationship
    R"""
    recPNl[[reps]] <- recPN
    recPRl[[reps]] <- recPR
    recDegl[[reps]] <- recDegree
    recClustl[[reps]] <- recClustLocal
    recPathl[[reps]] <- recPath
    recTraitNl[[reps]] <- recMedNTraits
    recTraitLl[[reps]] <- recMedMaxTraitLevel
    recNTraitsl[[reps]] <- recNTraits
    recBetaDivl[[reps]] <- recBetaDiv
    recPayl[[reps]] <- recPay
    recTraitDiversityl[[reps]] <- recTraitDiversity
    adjml[[reps]] <- adjm
    # recLearnedSocl[[reps]] <- recLearnedSoc
    # recLearnedInol[[reps]] <- recLearnedIno
    """

  end # END ALL REPETITIONS

  # SUMMARISING ALL DATA FROM ALL REPETITIONS AND STORING THEM IN AN *.RDATA FILE
  # At the end of all repetitions export data as an RData file
  if cluster
    # Cluster
    pat=string(paths, "/output/");
    if !isdir(pat)
      mkdir(pat)
    end
  else
    # Local
    pat="/Users/marco/Desktop/output/trial_";
  end

  outpath = string(pat,"pnpr_free_pay_rout_",queue);

  @rput outpath mmaxLevel mshared50 mshared90 mpay mntraits mbetadiv mpn mpr mtransmission mMedNTraits mMedMaxTraitLevel mTraitDiversity mdeg mpath  mclust mclustLocal; # mcomponentNumber mcomponentMedSize mcomponentMaxSize mlearnedsoc mlearnedino mRepertoireSimiliarityA mRepertoireSimiliarityR

  R"""
  setup$maxLevel <- mmaxLevel
  setup$shared50 <- mshared50
  setup$shared90 <- mshared90
  setup$recPay <- mpay
  setup$recNTraits <- mntraits
  setup$recBetaDiv <- mbetadiv
  setup$recPN <- mpn
  setup$recPR <- mpr
  setup$recTransmission <- mtransmission
  setup$recMedNTraits <- mMedNTraits
  setup$recMedMaxTraitLevel <- mMedMaxTraitLevel
  setup$recTraitDiversity <- mTraitDiversity
  setup$recDegree <- mdeg
  setup$path <- mpath
  #setup$clust <- mclust
  setup$clustLocal <- mclustLocal
  setup$generation <- round(setup[["nRound"]]/setup[["nod"]])

  # setup$componentNumber <- mcomponentNumber
  # setup$componentMedSize <- mcomponentMedSize
  # setup$componentMaxSize <- mcomponentMaxSize
  # setup$recLearnedSoc <- mlearnedsoc
  # setup$recLearnedIno <- mlearnedino
  # setup$recRepertoireSimiliarityA <- mRepertoireSimiliarityA
  # setup$recRepertoireSimiliarityR <- mRepertoireSimiliarityR

  setup=do.call(cbind,setup)
  save(setup, adjm, repertoir, recPNl, recPRl, recDegl, recClustl, recPathl, recTraitNl, recTraitLl, recNTraitsl, recBetaDivl, recPayl, recRep, recTraitDiversityl, adjml, file=outpath) #recLearnedSocl, recLearnedInol, relationship
  """

  print(".. done")
  return true
end
