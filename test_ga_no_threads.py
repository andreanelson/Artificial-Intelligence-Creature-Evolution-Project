# In test_ga_no_threads.py
import population
import simulation 
import genome 
import creature 
import numpy as np

pop = population.Population(pop_size=25, gene_count=5)
sim = simulation.Simulation()

for iteration in range(50):
    for cr in pop.creatures:
        cr.reset_position_history()  # Reset position history for stability calculation
        sim.run_creature(cr, 2400)
        
    # Calculate fitness for each creature
    fits = [cr.get_fitness() for cr in pop.creatures]
    links = [len(cr.get_expanded_links()) for cr in pop.creatures]
    dist = [cr.get_distance_to_center() for cr in pop.creatures if cr.positions]  # No need to pass position
    print(iteration, "fittest:", np.round(np.max(fits), 3),
          "mean:", np.round(np.mean(fits), 3), "mean links", np.round(np.mean(links)),
          "max links", np.round(np.max(links)), "dist to centre", np.round(np.min(dist), 3))

    # Create fitness map for parent selection
    fit_map = population.Population.get_fitness_map(fits)
    new_creatures = []
    for i in range(len(pop.creatures)):
        p1_ind = population.Population.select_parent(fit_map)
        p2_ind = population.Population.select_parent(fit_map)
        p1 = pop.creatures[p1_ind]
        p2 = pop.creatures[p2_ind]
        # now we have the parents!
        dna = genome.Genome.crossover(p1.dna, p2.dna)
        dna = genome.Genome.point_mutate(dna, rate=0.1, amount=0.25)
        dna = genome.Genome.shrink_mutate(dna, rate=0.25)
        dna = genome.Genome.grow_mutate(dna, rate=0.1)
        cr = creature.Creature(1)
        cr.update_dna(dna)
        new_creatures.append(cr)

    # Elitism: Keep the best creature from the previous generation
    max_fit = np.max(fits)
    for cr in pop.creatures:
        if cr.get_fitness() == max_fit:
            new_cr = creature.Creature(1)
            new_cr.update_dna(cr.dna)
            new_creatures[0] = new_cr
            filename = "elite_" + str(iteration) + ".csv"
            genome.Genome.to_csv(cr.dna, filename)
            break
    
    pop.creatures = new_creatures




# # In test_ga_no_threads.py
# import population
# import simulation 
# import genome 
# import creature 
# import numpy as np

# pop = population.Population(pop_size=22, gene_count=2)
# sim = simulation.Simulation()

# for iteration in range(100):
#     for cr in pop.creatures:
#         cr.reset_position_history()  # Reset position history for stability calculation
#         sim.run_creature(cr, 2400)
        
#     # Calculate fitness for each creature
#     fits = [cr.get_fitness() for cr in pop.creatures]
#     links = [len(cr.get_expanded_links()) for cr in pop.creatures]
#     print(iteration, "fittest:", np.round(np.max(fits), 3),
#           "mean:", np.round(np.mean(fits), 3), "mean links", np.round(np.mean(links)),
#           "max links", np.round(np.max(links)))

#     # Create fitness map for parent selection
#     fit_map = population.Population.get_fitness_map(fits)
#     new_creatures = []
#     for i in range(len(pop.creatures)):
#         p1_ind = population.Population.select_parent(fit_map)
#         p2_ind = population.Population.select_parent(fit_map)
#         p1 = pop.creatures[p1_ind]
#         p2 = pop.creatures[p2_ind]
#         # now we have the parents!
#         dna = genome.Genome.crossover(p1.dna, p2.dna)
#         dna = genome.Genome.point_mutate(dna, rate=0.1, amount=0.25)
#         dna = genome.Genome.shrink_mutate(dna, rate=0.25)
#         dna = genome.Genome.grow_mutate(dna, rate=0.1)
#         cr = creature.Creature(1)
#         cr.update_dna(dna)
#         new_creatures.append(cr)

#     # Elitism: Keep the best creature from the previous generation
#     max_fit = np.max(fits)
#     for cr in pop.creatures:
#         if cr.get_fitness() == max_fit:
#             new_cr = creature.Creature(1)
#             new_cr.update_dna(cr.dna)
#             new_creatures[0] = new_cr
#             filename = "elite_" + str(iteration) + ".csv"
#             genome.Genome.to_csv(cr.dna, filename)
#             break
    
#     pop.creatures = new_creatures








# # If you on a Windows machine with any Python version 
# # or an M1 mac with any Python version
# # or an Intel Mac with Python > 3.7
# # the multi-threaded version does not work
# # so instead, you can use this version. 


# import population
# import simulation 
# import genome 
# import creature 
# import numpy as np


# pop = population.Population(pop_size=20, 
#                             gene_count=6)
# #sim = simulation.ThreadedSim(pool_size=1)
# sim = simulation.Simulation()

# for iteration in range(10):
#     # this is a non-threaded version 
#     # where we just call run_creature instead
#     # of eval_population
#     for cr in pop.creatures:
#         sim.run_creature(cr, 2400)            
#     #sim.eval_population(pop, 2400)
#     fits = [cr.get_distance_travelled() 
#             for cr in pop.creatures]
#     links = [len(cr.get_expanded_links()) 
#             for cr in pop.creatures]
#     print(iteration, "fittest:", np.round(np.max(fits), 3), 
#             "mean:", np.round(np.mean(fits), 3), "mean links", np.round(np.mean(links)), "max links", np.round(np.max(links)))       
#     fit_map = population.Population.get_fitness_map(fits)
#     new_creatures = []
#     for i in range(len(pop.creatures)):
#         p1_ind = population.Population.select_parent(fit_map)
#         p2_ind = population.Population.select_parent(fit_map)
#         p1 = pop.creatures[p1_ind]
#         p2 = pop.creatures[p2_ind]
#         # now we have the parents!
#         dna = genome.Genome.crossover(p1.dna, p2.dna)
#         dna = genome.Genome.point_mutate(dna, rate=0.1, amount=0.25)
#         dna = genome.Genome.shrink_mutate(dna, rate=0.25)
#         dna = genome.Genome.grow_mutate(dna, rate=0.1)
#         cr = creature.Creature(1)
#         cr.update_dna(dna)
#         new_creatures.append(cr)
#     # elitism
#     max_fit = np.max(fits)
#     for cr in pop.creatures:
#         if cr.get_distance_travelled() == max_fit:
#             new_cr = creature.Creature(1)
#             new_cr.update_dna(cr.dna)
#             new_creatures[0] = new_cr
#             filename = "elite_"+str(iteration)+".csv"
#             genome.Genome.to_csv(cr.dna, filename)
#             break
    
#     pop.creatures = new_creatures
                            
