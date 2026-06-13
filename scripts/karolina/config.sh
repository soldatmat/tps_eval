SUBMIT_JOB="sbatch"

# Per-install local overrides (SLURM account, etc.) — NOT committed (config.local.sh
# is gitignored). Copy config.local.sh.example -> config.local.sh and set
# SBATCH_ACCOUNT to your IT4I project. The grant id is install-specific and changes
# on allocation renewal, so it must not live in this committed file.
# (/bin/sh is bash on Karolina, so BASH_SOURCE resolves this file's dir reliably.)
__cfgdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ -f "$__cfgdir/config.local.sh" ] && . "$__cfgdir/config.local.sh"
unset __cfgdir

# Karolina's default SLURM account ('defac') has NO qcpu/qgpu association, so a real
# project account is required or jobs are rejected. Warn loudly (stderr) if unset.
if [ -z "${SBATCH_ACCOUNT:-}" ]; then
    echo "[karolina/config.sh] WARNING: SBATCH_ACCOUNT is unset — SLURM jobs will fall back to the default 'defac' account (no qcpu/qgpu association) and be rejected. Fix: cp scripts/karolina/config.local.sh.example scripts/karolina/config.local.sh and set your IT4I project, or 'export SBATCH_ACCOUNT=<project>'." >&2
fi
