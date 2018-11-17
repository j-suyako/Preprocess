from common import *
import numpy as np


if __name__ == "__main__":
    bound = Cuboid(np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]))
    holes = list()
    bricks = list()
    output_path = ""
    models = Model(bound=bound, holes=holes)
    models.build()
    models.assign(bricks)
    models.plot()
    models.output(output_path)
