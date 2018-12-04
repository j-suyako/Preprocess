from common import Cuboid, Model
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d
import unittest


class TestModel(unittest.TestCase):
    def test_cut_method(self):
        cuboid1 = Cuboid(origin=(0, 0, 0), side=(1, 1, 1))
        cuboid2 = Cuboid(origin=(-0.5, -0.5, -0.5), side=(1, 1, 1))
        cuboid1.cut(cuboid2)
        model = Model(bound=bound, ele_size=0.3, holes=holes)
        model.build()
        for intercept, section in model.cross_section([0.5], 0):
            for index, plane in section:
                plane.plot()
        plt.show()


if __name__ == "__main__":
    unittest.main()
