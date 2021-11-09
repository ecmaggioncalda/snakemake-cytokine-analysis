# rule targets:
#     input:
#         'results/{cytokine}/{group}.{genome}.report.md'

rule preprocess_data:
    input:
        R="code/preproc.R",
        csv="data/mikropml/{cytokine}/{group}.{genome}.csv"
    output:
        rds='data/mikropml/{cytokine}/{group}.{genome}.dat_proc.Rds'
    log:
        "log/{cytokine}/{group}.{genome}.preprocess_data.txt"
    benchmark:
        "benchmarks/{cytokine}.{group}.{genome}.preprocess_data.txt"
    params:
        outcome_colname='{cytokine}'
    resources:
        ncores=ncores
    script:
        "code/preproc.R"

rule run_ml:
    input:
        R="code/ml.R",
        rds=rules.preprocess_data.output.rds
    output:
        model="results/{cytokine}/runs/{group}.{genome}.{method}_{seed}_model.Rds",
        perf=temp("results/{cytokine}/runs/{group}.{genome}.{method}_{seed}_performance.csv"),
        feat=temp("results/{cytokine}/runs/{group}.{genome}.{method}_{seed}_feature-importance.csv")
    log:
        "log/runs/run_ml.{cytokine}.{group}.{genome}.{method}_{seed}.txt"
    benchmark:
        "benchmarks/runs/run_ml.{cytokine}.{group}.{genome}.{method}_{seed}.txt"
    params:
        outcome_colname='{cytokine}',
        method="{method}",
        seed="{seed}",
        kfold=kfold
    resources:
        ncores=ncores
    script:
        "code/ml.R"

rule combine_results:
    input:
        R="code/combine_results.R",
        csv=expand("results/{{cytokine}}/runs/{{group}}.{{genome}}.{method}_{seed}_{{type}}.csv", method=ml_methods, seed=seeds)
    output:
        csv='results/{cytokine}/{group}.{genome}.{type}_results.csv'
    log:
        "log/{cytokine}/{group}.{genome}.combine_results_{type}.txt"
    benchmark:
        "benchmarks/{cytokine}/{group}.{genome}.combine_results_{type}.txt"
    script:
        "code/combine_results.R"

rule combine_hp_performance:
    input:
        R='code/combine_hp_perf.R',
        rds=expand('results/{{cytokine}}/runs/{{group}}.{{genome}}.{{method}}_{seed}_model.Rds', seed=seeds)
    output:
        rds='results/{cytokine}/{group}.{genome}.hp_performance_results_{method}.Rds'
    log:
        "log/{cytokine}/{group}.{genome}.combine_hp_perf_{method}.txt"
    benchmark:
        "benchmarks/{cytokine}/{group}.{genome}.combine_hp_perf_{method}.txt"
    script:
        "code/combine_hp_perf.R"

rule combine_benchmarks:
    input:
        R='code/combine_benchmarks.R',
        # tsv=expand(rules.run_ml.benchmark, method = ml_methods, seed = seeds)
        tsv=expand("benchmarks/runs/run_ml.{{cytokine}}.{{group}}.{{genome}}.{method}_{seed}.txt", method=ml_methods, seed=seeds)
    output:
        csv='results/{cytokine}/{group}.{genome}.benchmarks_results.csv'
    log:
        'log/{cytokine}/{group}.{genome}.combine_benchmarks.txt'
    script:
        'code/combine_benchmarks.R'

rule plot_performance:
    input:
        R="code/plot_perf.R",
        csv='results/{cytokine}/{group}.{genome}.performance_results.csv'
    output:
        plot='figures/{cytokine}/{group}.{genome}.performance.png'
    log:
        "log/{cytokine}/{group}.{genome}.plot_performance.txt"
    script:
        "code/plot_perf.R"

rule plot_hp_performance:
    input:
        R='code/plot_hp_perf.R',
        rds=rules.combine_hp_performance.output.rds,
    output:
        plot='figures/{cytokine}/{group}.{genome}.hp_performance_{method}.png'
    log:
        'log/{cytokine}/{group}.{genome}.plot_hp_perf_{method}.txt'
    script:
        'code/plot_hp_perf.R'

rule plot_benchmarks:
    input:
        R='code/plot_benchmarks.R',
        csv=rules.combine_benchmarks.output.csv
    output:
        plot='figures/{cytokine}/{group}.{genome}.benchmarks.png'
    log:
        'log/{cytokine}/{group}.{genome}.plot_benchmarks.txt'
    script:
        'code/plot_benchmarks.R'

rule render_report:
    input:
        Rmd='report.Rmd',
        R='code/render.R',
        perf_plot=rules.plot_performance.output.plot,
        #hp_plot=expand(rules.plot_hp_performance.output.plot, method = ml_methods),
        hp_plot=expand('figures/{{cytokine}}/{{group}}.{{genome}}.hp_performance_{method}.png', method=ml_methods),
        bench_plot=rules.plot_benchmarks.output.plot
    output:
        doc='results/{cytokine}/{group}.{genome}.report.md'
    log:
        "log/{cytokine}/{group}.{genome}.render_report.txt"
    params:
        nseeds=nseeds,
        ml_methods=ml_methods,
        ncores=ncores,
        kfold=kfold
    script:
        'code/render.R'

rule clean:
    input:
        rules.render_report.output,
        rules.plot_performance.output.plot,
        rules.plot_benchmarks.output.plot
    shell:
        '''
        rm -rf {input}
        '''
