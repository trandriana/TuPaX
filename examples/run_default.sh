#!/usr/bin/env bash
set -euo pipefail
python -m gray_scott.compute --nx 128 --T 2000 --steps 4000 --n-image-snapshots 6
