configfile: 'config/config.yml'

ncores = config['ncores']
ml_methods = config['ml_methods']
kfold = config['kfold']
outcome_colname = config['outcome_colname']

nseeds = config['nseeds']
start_seed = 100
seeds = range(start_seed, start_seed + nseeds)

rule targets:
    input:
        'report.md'

#include: "rules/summarize.smk", this should probably actually be a directory up instead of below, but I didn't think it through this way initially
#alternatively could create a sub-workflow snakefile that just depends on the outputs, need to figure out which is preferable workflow-wise
