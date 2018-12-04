from common import Cuboid
import unittest


class TestCuboid(unittest.TestCase):
    def test_cut_method(self):
        cuboid1 = Cuboid(origin=(0, 0, 0), side=(1, 1, 1))
        cuboid2 = Cuboid(origin=(-0.5, -0.5, -0.5), side=(1, 1, 1))

        poly1, poly2 = cuboid1.cut(cuboid2)
        poly1.plot()
        poly2.plot()


if __name__ == "__main__":
    unittest.main()
