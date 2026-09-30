#!/usr/bin/env bash
# Pull the container images in config/config.yaml (`containers`) to .sif files.
# Run once on a login node before `sbatch run_workflow.sbatch`: BioCloud cannot
# build images inside SLURM jobs, and turning a docker:// image into a .sif
# file is a build. Images that are already there are skipped.
#
#   conda activate snakemake_rhizopus     # for python + pyyaml
#   bash pull_containers.sh
set -euo pipefail
cd "$(dirname "$0")"

python - <<'PY' | while read -r sif image; do
import yaml

containers = yaml.safe_load(open("config/config.yaml"))["containers"]
for name, image in containers["images"].items():
    print(f"{containers['folder']}/{name}.sif", image)
PY
    if [ -s "$sif" ]; then
        echo "exists: $sif"
        continue
    fi
    mkdir -p "$(dirname "$sif")"
    echo "pulling $image -> $sif"
    apptainer pull --force "$sif.part" "$image" </dev/null
    mv "$sif.part" "$sif"
done
