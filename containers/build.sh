#!/bin/bash
# Build both images locally (Docker), and optionally convert them for Singularity/Apptainer.
# Released images are published automatically to ghcr.io by .github/workflows/containers.yml;
# local builds are only needed for development or offline use.
#   bash containers/build.sh                 # docker images poreflow-its-main:0.2.0, poreflow-its-medaka:0.2.0
#   bash containers/build.sh --sif           # also write poreflow-its-main.sif / poreflow-its-medaka.sif
set -euo pipefail
cd "$(dirname "$0")/.."
VERSION=0.2.0
for e in main medaka; do
    docker build --platform linux/amd64 -f containers/Dockerfile.$e -t poreflow-its-$e:$VERSION .
    if [[ "${1:-}" == --sif ]]; then
        singularity build poreflow-its-$e.sif docker-daemon://poreflow-its-$e:$VERSION 2>/dev/null \
            || apptainer build poreflow-its-$e.sif docker-daemon://poreflow-its-$e:$VERSION
    fi
done
