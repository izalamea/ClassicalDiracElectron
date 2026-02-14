#!/usr/bin/env python3
"""Plot 3D trajectory from CDE_4thRK output (columns 4,5,6 = x,y,z)."""
import argparse
import numpy as np

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    print("Need matplotlib: pip install matplotlib")
    raise SystemExit(1)

def main():
    p = argparse.ArgumentParser(description="Plot 3D trajectory from CDE_4thRK output.")
    p.add_argument("data", nargs="?", default="Trajectories/CE_trajectory_RK_0_.dat", help="Path to .dat file")
    p.add_argument("-o", "--output", help="Save figure to file instead of showing")
    p.add_argument("-t", "--title", default="Classical Dirac Electron trajectory", help="Plot title")
    args = p.parse_args()

    d = np.loadtxt(args.data)
    x, y, z = d[:, 3], d[:, 4], d[:, 5]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(x, y, z, "b-", lw=0.5)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(args.title)
    plt.tight_layout()
    if args.output:
        plt.savefig(args.output, dpi=150)
    else:
        plt.show()

if __name__ == "__main__":
    main()
