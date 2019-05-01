# consensus-variant-caller

## Quick Start
The tool(s) can be run on:
* [the command line](#running-on-the-command-line)
* [Deploit](#deploit) (recommended)

## Running on the command line

### Dependencies 
[Nextflow](https://www.nextflow.io/)
[Docker](https://www.docker.com/)

An example run of the pipeline on the command line may look like this:

```bash
nextflow run lifebit-ai/consensus-variant-caller --batch s3://lifebit-featured-datasets/pipelines/consensus-variant-caller/23-04-19-nf-testingBatch.csv --genome hg19 --annovar_protocols refGene,exac03,avsnp147,dbnsfp30a --annovar_operation gx,f,f,f --annotation s3://lifebit-featured-datasets/lifebit-featured-datasets/pipelines/consensus-variant-caller/humandb
```
**Warning: this pipeline requires at least 14GB of memory & 4 CPUs to run. It will also use a lot of storage & may take a short while to complete**

## Deploit

Deploit is a bioinformatics platform, developed by Lifebit, where you can run your analysis over the Cloud/AWS.

It is free for indivudal users to [sign-up](https://deploit.lifebit.ai/register)

To run the pipeline once logged in navigate to the pipelines page:

![deploit](https://raw.githubusercontent.com/lifebit-ai/ecw-converter/master/images/deploit.png)

### Running on Deploit

![run_vc_deploit](https://raw.githubusercontent.com/lifebit-ai/consensus-variant-caller/master/images/run_vc_deploit.gif)

The pipeline can then be run in three simple steps:
1. **Find & select the Consensus Variant Caller pipeline** this can be searched for in the `PUBLIC PIPELINES & TOOLS` section
2. **Select input data & parameters** you can click the `Try with example data & parameters` button for examples
3. **Select a project & instance** before running the analysis you must select the instance. This determines the available resources & cost. The project is used for grouping/organising multiple jobs

See [documentation](https://lifebit.gitbook.io/deploit/pipelines-documentations-and-examples-1/nextflow-pipelines/consensus-variant-caller) for running the pipeline over Deploit:
[![docs](https://raw.githubusercontent.com/lifebit-ai/consensus-variant-caller/master/images/docs.png)](https://lifebit.gitbook.io/deploit/pipelines-documentations-and-examples-1/nextflow-pipelines/consensus-variant-caller)

See [example run](https://deploit.lifebit.ai/public/jobs/5cc2e65702877100b2bb96e6) of the pipeline on Deploit:
[![job](https://raw.githubusercontent.com/lifebit-ai/consensus-variant-caller/master/images/job.png)](https://lifebit.gitbook.io/deploit/pipelines-documentations-and-examples-1/nextflow-pipelines/consensus-variant-caller)
