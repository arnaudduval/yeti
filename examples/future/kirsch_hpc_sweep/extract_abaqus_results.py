"""Run THIS with Abaqus's own bundled Python, not a normal python3 -- it needs
`odbAccess`, which only exists inside Abaqus:

    abaqus python extract_abaqus_results.py <job>.odb <output>.csv

Writes a plain CSV (node_id,x,y,ux,uy) for every node's final-frame
displacement in the .odb -- read back by `kirsch_sweep.py abaqus-results`,
which does NOT need odbAccess or Abaqus itself.

Written for both Abaqus's Python 2.7 (older versions) and Python 3 (newer
versions) since the Abaqus version on the target machine wasn't known when
this was written -- avoid f-strings/version-specific syntax here. If this
still fails on your Abaqus version, the two likely culprits are the odb's
step/instance naming (see the two loops below, already written to not
hard-code either) or a print-related syntax error -- report the exact
traceback back so this can be adjusted, since it could not be tested against
a real Abaqus installation.
"""
from __future__ import division, print_function

import sys

from odbAccess import openOdb


def main(odb_path, csv_path):
    odb = openOdb(path=odb_path, readOnly=True)
    try:
        step_names = list(odb.steps.keys())
        if not step_names:
            raise RuntimeError("no steps found in " + odb_path)
        step = odb.steps[step_names[-1]]

        frames = step.frames
        if len(frames) == 0:
            raise RuntimeError("no frames found in the last step of " + odb_path)
        frame = frames[-1]  # final (converged) increment

        instance_names = list(odb.rootAssembly.instances.keys())
        if not instance_names:
            raise RuntimeError("no instances found in " + odb_path)
        instance = odb.rootAssembly.instances[instance_names[0]]

        node_coords = {}
        for node in instance.nodes:
            coords = node.coordinates
            node_coords[node.label] = (float(coords[0]), float(coords[1]))

        disp_field = frame.fieldOutputs["U"]
        disp_subset = disp_field.getSubset(region=instance)

        rows = []
        for value in disp_subset.values:
            label = value.nodeLabel
            data = value.data
            ux, uy = float(data[0]), float(data[1])
            x, y = node_coords[label]
            rows.append((label, x, y, ux, uy))

        if not rows:
            raise RuntimeError("no nodal U values found -- check the .odb has "
                               "*OUTPUT, FIELD / *NODE OUTPUT / U requested")

        rows.sort(key=lambda r: r[0])

        f = open(csv_path, "w")
        try:
            f.write("node_id,x,y,ux,uy\n")
            for (label, x, y, ux, uy) in rows:
                f.write("%d,%.10g,%.10g,%.10g,%.10g\n" % (label, x, y, ux, uy))
        finally:
            f.close()

        print("wrote %s (%d nodes)" % (csv_path, len(rows)))
    finally:
        odb.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("usage: abaqus python extract_abaqus_results.py <job>.odb <output>.csv\n")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
