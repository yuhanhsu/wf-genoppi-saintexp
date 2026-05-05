# wf-genoppi-saintexp

Files for running Genoppi and SAINTexpress IP-MS analysis workflow on Terra
- R scripts (runGenoppi.r, runSAINTexpressInt.r, runSAINTexpressSpc.r) to run analysis
- Dockerfile to containerize the R scripts and install dependencies
- cloudbuild.yaml to build and push docker image to Google Artifact Registry via Cloud Build
- genoppi-saintexp.wdl to define the workflow to be run using the docker image
- .dockstore.yml to sync workflow WDL to [Dockstore](https://dockstore.org/workflows/github.com/yuhanhsu/wf-genoppi-saintexp/genoppi-saintexp:main?tab=info)

Approximate run time
- Cloud Build: ~15 minutes
- Dockstore: a few minutes
- Terra workflow: ~5 minutes

