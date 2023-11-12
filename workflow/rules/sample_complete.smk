def get_input_files(wildcards):
    _sruns = samples.loc[wildcards.samp]['runid']
    _layouts = runs.loc[_sruns]['layout']
    if not all(_layouts.iloc[0] == _ for _ in _layouts):
        msg = f'All runs for sample "{wildcards.samp}" do not have the same layout (PAIRED or SINGLE)'
        raise WorkflowError(msg)

    if _layouts.iloc[0] == 'SINGLE':
        return expand("results/{{dataset}}/fastq/{{samp}}/{sra_run}_1.fastq.gz", sra_run=_sruns)
    elif _layouts.iloc[0] == 'PAIRED':
        return expand("results/{{dataset}}/fastq/{{samp}}/{sra_run}_{rn}.fastq.gz", sra_run=_sruns, rn=[1,2])
    else:
        raise WorkflowError(f'Invalid layout: "{_layouts[0]}"')

rule sample_complete:
    output:
        touch("results/{dataset}/complete/{samp}/{samp}_complete.txt")
    input:
        get_input_files

