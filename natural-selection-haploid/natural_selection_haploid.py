"""
Simulate a simple model for natural selection in a haploid population.
Given:
Initial:
  population size N
  two alleles with frequencies f1 and f2 
  relative fitnesses of the two alleles w1 and w2
    - w1 and w2 are expected number of offsprings per individual
  Assume N*f1*w1 and N*f2*w2 are expected values of populations in the next generation obtained by sampling from Poisson distributions of these lambdas.

Output:
  N and f1 for gens generations in list of lists
  - plot N1 and f1 versus generations
  - savfe plot and data in files
"""
import numpy as np  # we will use np.random.poisson() to sample from Poisson
import matplotlib.pyplot as plt 

def natural_selection_haploid(N, f1, w1, w2, gens=1, output_file=None):
  #initialize
  N1 = N*f1
  N2 = N-N1
  f2 = N2/N
  population_and_allele1freq = [[N,f1]]
  
  for _ in range(gens):
    # generate populations for the two allele by Poisson sampling
    [N1, N2] = np.random.poisson(lam=(N*f1*w1, N*f2*w2), size=2 )
    N = N1 + N2
    f1 = N1/(N1+N2)
    f2 = 1-f1
    # Record the data
    population_and_allele1freq.append([N, f1])
  if output_file:
    with open(output_file, "w") as filehandle:
      for listitem in population_and_allele1freq:
        filehandle.write('%s\n' % listitem)
  return population_and_allele1freq

if __name__ == "__main__":
  import os
  
  # Give parameters for simulation and saving data to file
  output_file="data.txt"
  savefig_file="natural_selection_haploid.png"
  gens = 1000
  N=10000
  f1=0.5
  w1=1.0001
  w2=0.9999

  # simulate
  scriptdir = os.path.dirname(os.path.realpath(__file__))
  population_and_allele1freq = natural_selection_haploid(N, f1, w1, w2, gens, scriptdir+"/"+output_file)

  # Plot
  # gather data to plot
  generations = range( len(population_and_allele1freq) )
  population1 = [pf[0]*pf[1] for pf in population_and_allele1freq]
  freq1 = [pf[1] for pf in population_and_allele1freq]

  # Plot N1 and f1 - share x-axis, which is generations
  fig, ax1 = plt.subplots()
  color = 'tab:red'
  ax1.set_title("Simple Model of Haploid Natural Selection")
  ax1.set_xlabel('Generations')
  ax1.set_ylabel('Population of Allele 1', color=color)
  ax1.plot(generations, population1, color=color)
  ax1.tick_params(axis='y', labelcolor=color)

  ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

  color = 'tab:blue'
  ax2.set_ylabel('Frequency of Allele 1', color=color)  # we already handled the x-label with ax1
  ax2.plot(generations, freq1, color=color)
  ax2.tick_params(axis='y', labelcolor=color)

  fig.tight_layout()  # otherwise the right y-label is slightly clipped
  plt.savefig(scriptdir+"/"+savefig_file)
  plt.show()






  




