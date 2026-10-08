#!/bin/bash
# Build both images (Docker), and optionally convert them for Singularity/Apptainer.
#   bash containers/build.sh                 # docker images poreflow-its-main:0.1.0, poreflow-its-medaka:0.1.0
#   bash containers/build.sh --sif           # also write poreflow-its-main.sif / poreflow-its-medaka.sif
set -euo pipefail
cd "$(dirname "$0")/.."
VERSION=0.1.0
for e in main medaka; do
    docker build --platform linux/amd64 -f containers/Dockerfile.$e -t poreflow-its-$e:$VERSION .
    if [[ "${1:-}" == --sif ]]; then
        singularity build poreflow-its-$e.sif docker-daemon://poreflow-its-$e:$VERSION 2>/dev/null \
            || apptainer build poreflow-its-$e.sif docker-daemon://poreflow-its-$e:$VERSION
    fi
done
