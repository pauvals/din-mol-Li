#!/usr/bin/env python3
import csv
import statistics
import sys
from collections import defaultdict


def fail(message):
    raise SystemExit(f"error: {message}")


if len(sys.argv) < 2:
    fail("usage: summarize_scaling.py scaling.csv [scaling.csv ...]")

runs = defaultdict(list)
systems = set()
for path in sys.argv[1:]:
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (int(row["nodes"]), int(row["images"]), int(row["threads"]))
            systems.add((int(row["particles"]), int(row["steps"])))
            runs[key].append(
                (
                    float(row["launcher_wall_seconds"]),
                    float(row["dana_wall_seconds"]),
                    int(row["residual_overlaps"]),
                )
            )

if not runs:
    fail("no benchmark rows found")
if len(systems) != 1:
    fail(f"cannot compare different systems: {sorted(systems)}")

particles, steps = systems.pop()
baseline_key = (1, 1, 1)
baseline = None
if baseline_key in runs:
    baseline = statistics.median(value[0] for value in runs[baseline_key])

writer = csv.writer(sys.stdout, lineterminator="\n")
writer.writerow(
    [
        "nodes",
        "images",
        "threads",
        "total_cores",
        "particles",
        "steps",
        "repeats",
        "median_launcher_wall_seconds",
        "median_dana_wall_seconds",
        "speedup_vs_1x1",
        "max_residual_overlaps",
    ]
)

for nodes, images, threads in sorted(runs, key=lambda item: (item[0], item[1] * item[2], item[1])):
    values = runs[(nodes, images, threads)]
    launcher = statistics.median(value[0] for value in values)
    dana = statistics.median(value[1] for value in values)
    speedup = "" if baseline is None else f"{baseline / launcher:.6f}"
    writer.writerow(
        [
            nodes,
            images,
            threads,
            images * threads,
            particles,
            steps,
            len(values),
            f"{launcher:.6f}",
            f"{dana:.6f}",
            speedup,
            max(value[2] for value in values),
        ]
    )
