#!/usr/bin/env bash
# Refresh the backup mirrors of the third-party vendor submodules.
#
# WHY: vendor/ pins third-party dependencies whose upstreams could disappear
# (repo deletion, account removal, DMCA), which would break fresh
# `git clone --recurse-submodules` / `git submodule update --init`. We keep
# independent mirrors (NOT GitHub forks — a fork shares the upstream fork
# network and can be disabled with it) on the soldatmat GitHub account as
# insurance. The `.gitmodules` URLs still point at upstream; the mirrors are
# only a fallback.
#
# WHEN TO RUN: only when you bump a vendored submodule pin to a newer upstream
# commit. The mirrors are point-in-time snapshots — nothing else makes them
# stale, so there is deliberately no scheduled sync. If an upstream dies, swap
# the URL in .gitmodules to the mirror (see tps_eval/CLAUDE.md, vendor gotcha).
#
# Only the two THIRD-PARTY submodules are mirrored here. `cif_to_pdb` and
# `pymol_scripts` are already soldatmat-owned upstreams, so they need no mirror.
#
# Requirements: git, and the `gh` CLI authenticated as soldatmat (used only to
# auto-create the mirror repo the first time). Pushes go over HTTPS.
set -euo pipefail

# upstream_url|mirror (owner/repo on github.com)
MIRRORS=(
  "https://github.com/dauparas/ProteinMPNN.git|soldatmat/ProteinMPNN"
  "https://bitbucket.org/lcbio/aggrescan3d.git|soldatmat/aggrescan3d"
)

workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT

for entry in "${MIRRORS[@]}"; do
  upstream="${entry%%|*}"
  mirror="${entry##*|}"
  name="${mirror##*/}"
  echo "=== ${name}: ${upstream} -> github.com/${mirror} ==="

  # Create the mirror repo if it does not exist yet (idempotent).
  if ! gh repo view "$mirror" >/dev/null 2>&1; then
    echo "  creating github.com/${mirror} ..."
    gh repo create "$mirror" --public \
      --description "Backup mirror of ${upstream} (vendored submodule dependency of tps_eval)"
  fi

  echo "  mirror-cloning upstream ..."
  git clone --quiet --mirror "$upstream" "${workdir}/${name}.git"

  # Push branches + tags only. Upstream GitHub mirror clones also pull
  # refs/pull/* (a read-only namespace GitHub rejects on push), so we do NOT
  # use `push --mirror`; explicit refspecs keep the backup clean.
  echo "  pushing refs/heads/* and refs/tags/* ..."
  git -C "${workdir}/${name}.git" push --quiet --force \
    "https://github.com/${mirror}.git" \
    'refs/heads/*:refs/heads/*' 'refs/tags/*:refs/tags/*'

  rm -rf "${workdir}/${name}.git"
  echo "  done."
done

echo "All vendor backup mirrors refreshed."
