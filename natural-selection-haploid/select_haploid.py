"""
Simulate a simple model for natural selection in a haploid population.
Given:
Initial:
  population size N
  two alleles with frequencies f1 and f2 
  relative fitnesses of the two alleles w1 and w2
    - w1 and w2 are expected number of offsprings per individual
  Assume w1 and w2 are expected values of two Poisson distributions.
    - samplig from Poisson will produce integer numbers of offsprings

Output:
  N1 and N2 for gens generations in list of lists
  - plot N1 and f1 versus generations
  - save plot and data in files
A good reference:
https://www.nature.com/scitable/knowledge/library/natural-selection-genetic-drift-and-gene-flow-15186648/
"""
import numpy as np  # we will use np.random.poisson() to sample from Poisson
import matplotlib.pyplot as plt 

from tqdm import tqdm # for progress monitoring
from random_colormap import rand_cmap

def select_haploid(N, f1, w1, w2, gens=1, output_file=None):
  """
  Simulates the natural selection process in a haploid.
  Args:
    N = initial population size
    f1 = frequency of allele 1 
    w1 = relative fitness of allele 1
    w2 = relative fitness of allele 2
    gens = number of generations to simulate
    oputput_file = name of the outfile if saving the output desired
  Returns:
    populations: List of [allele 1 population, allele 2 population] for each generation

  """
  #initialize
  N1 = int(N*f1)
  N2 = N-N1
  populations = [[N1,N2]]

  for _ in tqdm(range(gens)): # run using a progress bar
    # generate populations for the two allele by Poisson sampling
    N1 = np.sum(np.random.poisson(lam=w1, size=N1 ))
    N2 = np.sum(np.random.poisson(lam=w2, size=N2 ))
    # Record the data
    populations.append([N1, N2])

  if output_file:
    with open(output_file, "w") as filehandle:
      for pops in populations:
        filehandle.write('%s\n' % pops)

  return populations

def simulate(initial_values, num_simulations=2):
  """
  Simulates haploid simple selection for num_simulations times.
  Args:
    initial_values: dictionary with keys:
      output_file, savefig_file, generations, initial_population, initial_freq1, fitness_1, fitness2
  Returns:
    result:
      List of list(simulation number) of list([N1, N2] for each generation)
  
  """
  output_file=initial_values["output_file"]
  
  gens = initial_values["generations"]
  N=initial_values["initial_population"]
  f1=initial_values["initial_freq1"]
  w1=initial_values["fitness1"]
  w2=initial_values["fitness2"]
  scriptdir = os.path.dirname(os.path.realpath(__file__))

  output_file = scriptdir+"/"+output_file
  
  result = []
  for sim in tqdm(range(num_simulations)):
    populations = select_haploid(N, f1, w1, w2, gens, None) # don't write yet
    result.append(populations)
  if output_file:
    with open(output_file, "w") as filehandle:
      for pops in result:
        filehandle.write('%s\n' % pops)  
  return result


if __name__ == "__main__":
  import os
  
  # Give parameters for simulation and saving data to file
  num_simulations = 5

  output_file="data.txt"
  savefig_file="natural_selection_haploid.png"
  gens = 50
  N=1000
  f1=0.5
  w1=1.00
  w2=0.01
  initial_values = {
    "output_file" : output_file,
    "generations" : gens,
    "initial_population" : N,
    "initial_freq1" : f1,
    "fitness1" : w1,
    "fitness2" : w2,
  }

  result = simulate(initial_values, num_simulations)
  #print(result)

  # # simulate
  scriptdir = os.path.dirname(os.path.realpath(__file__))
  # populations = select_haploid(N, f1, w1, w2, gens, scriptdir+"/"+output_file)

  # Plot
  fig1, ax1 = plt.subplots(figsize=[10,6])
  fig2, ax2 = plt.subplots(figsize=[10,6])
  #ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
  prop_cycle = plt.rcParams['axes.prop_cycle']
  colors = prop_cycle.by_key()['color']

  def get_cmap(n, name='hsv'):#from https://stackoverflow.com/questions/14720331/how-to-generate-random-colors-in-matplotlib
      '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
      RGB color; the keyword argument name must be a standard mpl colormap name.'''
      return plt.cm.get_cmap(name, n)
  cmap = rand_cmap(len(result), type='bright', verbose=False)
  for sim_num, populations in enumerate(result):
    # gather data to plot
    generations = range( len(populations) )
    population1 = [pf[0] for pf in populations]
    freq1 = [pf[0]/(pf[0]+pf[1]) for pf in populations]
    population2 = [pf[1] for pf in populations]
    # Plot N1 and f1 - share x-axis, which is generations
    
    color = cmap(sim_num)
    line1, = ax1.plot(generations, population1, lw=3, alpha=0.8, color=color, label="N1")
    ax1.text(generations[-1], population1[-1], f"{sim_num}")
    
    #color = 'tab:green'
    line2, = ax1.plot(generations, population2, dashes=[6, 2], color=color, label="N2")
    ax1.text(generations[-1], population2[-1], f"{sim_num}")

    color = 'tab:blue'
    line3, = ax2.plot(generations, freq1, color=color, label="f1")

    #fig1.tight_layout()  # otherwise the right y-label is slightly clipped
    fig2.tight_layout()  # otherwise the right y-label is slightly clipped
  
  ax1.set_title("Simple Model of Haploid Natural Selection")
  ax1.set_xlabel('Generations')
  ax1.set_ylabel('Populations', color=color)
  ax1.tick_params(axis='y', labelcolor=color) 

  ax2.set_ylabel('Frequency of Allele 1', color=color)  
  ax2.tick_params(axis='y', labelcolor=color)

  fig1.legend((line1, line2), ("N1", "N2"))
  fig1.savefig(scriptdir+"/"+savefig_file)

  plt.show()






  




