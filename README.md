We consider a three-opinion model with conservative learning, where individuals can hold a positive (+), negative (−), or neutral (0) opinion. The nonlinear learning behavior is controlled by the parameter $\alpha$.

We investigate the polarization probability in well-mixed populations and complex networks, including Erdős–Rényi (ER), Watts–Strogatz (WS), and Barabási–Albert (BA) networks.

(1) Monte Carlo Simulations

voter_wellmixed.m: Simulations in well-mixed populations.

voter_static_polarization_ER.m: Simulations on Erdős–Rényi (ER) networks.

voter_static_polarization_WS.m: Simulations on Watts–Strogatz (WS) networks.

voter_static_polarization_BA.m: Simulations on Barabási–Albert (BA) networks.

(2) Evolutionary Game Theory
   
fixation_curve.m: Computes fixation probability over different initial fractions of strategy A.

fixation_probability_evolutionary_game_theory.m: Calculates fixation probabilities in the Moran process.

(3) Parameter Fitting

unified_best_fitting_three_alphas.m: Parameter fitting for three selected values of $\alpha$.

unified_best_fitting_wellmixed_all_alphas.m: Parameter fitting other $\alpha$ values in well-mixed populations.

unified_best_fitting_networks.m: Parameter fitting for different network structures.

unified_best_fitting_networks_all_alphas.m: Parameter fitting other $\alpha$ values on different networks.

In all the three network models, the population size is $N=100$ and the average degree of $\bar{k}=8$. The number of undirected links is $L=400$.
The three networks are generated as follows:
(a) ER networks: Firstly, we generate a random tree to ensure connectivity. Then, additional links are randomly added between unconnected node pairs until the number of links ($L=400$) is reached.

(b) WS networks: Firstly, we generate a regular ring lattice in which each node is connected to its $\bar{k}=8$ nearest neighbors (four on each side). We assume that the rewiring probability is $p=0.3$. The new endpoint is selected randomly from the non-neighbors of the focal node. 

(c) BA networks: Firstly, we generate a complete network with $\bar{k}/2+1=5$ nodes. Then, we add additional nodes and links one by one. Each node connects to $\bar{k}/2=4$ existing nodes according to the priority connection principle. The probability of selecting an existing node is proportional to the existing node’s degree. We obtain $\bar{k}\left(\bar{k}/2+1\right)/4+\left(N-\bar{k}/2-1\right)\ast \left(\bar{k}/2\right)=390$ links. Then randomly add $L-390=10$ links to the network. 

For each combination of $\alpha$ and initial fraction of extreme opinions, we run 5000 runs (giving rise to a sample size of 5000). A new network is generated for each realization, and the initial opinions are randomly assigned. 
