# Use
- Use `tps_eval/scripts/submit_job.sh` to submit jobs on any cluster.
- Every evaluation tool has its main run script created as `tps_eval/scripts/run_<tool_name>`.
- The main run script can be run from a job script on a computational cluster created as `tps_eval/scripts/<cluster>/<jobs>/<tool_name>.sh`.

# Adding clusters
- All code specific to a computational cluster should be contained within `tps_eval/scripts/<cluster>`.
- Create a file `tps_eval/scripts/<cluster>/config.sh` with cluster-specific settings.
- Create individual job scripts as `tps_eval/scripts/<cluster>/jobs/<job_name>.sh`.

## config.sh
File `tps_eval/scripts/<cluster>/config.sh` needs to include the following:
```sh
SUBMIT_JOB="<command used to submit jobs>" # i.e. "sbatch", "qsub"
```
