import cProfile, pstats
cProfile.run('solution = runCase("TosiFig3.2")', 'profileRes')
prof = pstats.Stats('profileRes')
prof.sort_stats('tottime').print_stats(100)
prof.sort_stats('cumtime').print_stats(100)
prof.sort_stats('calls').print_stats(100)