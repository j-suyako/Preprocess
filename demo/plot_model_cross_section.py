from common import Cuboid, Block
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d


if __name__ == "__main__":
    bound = Cuboid(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]))
    holes = list()
    model = Block(bound=bound, ele_size=0.3, holes=holes)
    model.initial()
    ax1 = plt.subplot2grid((1, 2), (0, 0), projection="3d")
    ax2 = plt.subplot2grid((1, 2), (0, 1), projection="3d")
    for intercept, section in model.cross_section([0.5], 0):
        for index, plane in section:
            plane.plot(ax1)
    plt.show()
