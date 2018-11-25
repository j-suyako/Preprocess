from config import LOGGER
from common import *
import numpy as np


if __name__ == "__main__":
    bound = Cuboid()
    holes = list()
    bricks = list()
    output_path = ""
    models = Model(bound=bound, holes=holes)
    models.build()
    models.assign(bricks)
    models.plot()
    models.output(output_path)
