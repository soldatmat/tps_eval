SUBMIT_JOB="qsub"

# MetaCentrum runs PBS Pro and is the Czech national grid: there is NO per-project
# core-hour / GPU-hour budget and no SLURM-style account, so (unlike karolina) there
# is no config.local.sh / SBATCH_ACCOUNT to source. The only gate is the queue, and
# PBS auto-routes by walltime / ngpus. submit_job.sh only needs $SUBMIT_JOB from here.
