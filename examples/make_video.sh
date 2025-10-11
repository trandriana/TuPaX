#!/usr/bin/env bash
set -euo pipefail

gray-scott-animate output_solution/snapshots.npz --out gray_scott.mp4 --fps 20 --which both --use-fixed-limits
