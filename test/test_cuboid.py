from common import Cuboid
import numpy as np
import unittest


class TestCuboid(unittest.TestCase):
    def setUp(self):
        self.cuboid = Cuboid(origin=(0, 0, 0), side=(1, 1, 1))

    def test_volume(self):
        self.assertTrue(self.cuboid.volume() == 1)

    def test_centroid(self):
        self.assertTrue(np.all(self.cuboid.centroid() == np.array([0.5, 0.5, 0.5])))

    def test_pan(self):
        cuboid = Cuboid(origin=(0, 0, 0), side=(1, 1, 1))
        vector = (0.5, 0.5, 0.5)
        self.assertTrue(np.all(cuboid.pan(vector).origin == np.array([0.5, 0.5, 0.5])))

    def test_rotate(self):
        cuboid = Cuboid(origin=(0, 0, 0), side=(1, 1, 1))
        direction = "pos_xy"
        centroid = np.array([0, 0, 0])
        self.assertTrue(np.all(cuboid.rotate(centroid, direction).origin == np.array([-1, 0, 0])))

    def test_contains(self):
        point = np.array([0.5, 0.5, 0.5])
        self.assertTrue(self.cuboid.contains(point))

    def test_is_touching(self):
        cuboid = Cuboid(origin=(0, 0, 1), side=(1, 1, 1))
        self.assertTrue(self.cuboid.is_touching(cuboid))

    def test_clarify_bound_with(self):
        Cuboid.ele_size = 0.2
        cuboid1 = Cuboid(origin=(0, 0, 0), side=(1, 1, 1))
        cuboid2 = Cuboid(origin=(0, 0, 1), side=(1, 1, 1))
        cuboid3 = Cuboid(origin=(1, 0, 0), side=(1, 1, 1))
        cuboid4 = Cuboid(origin=(0, 1, 0), side=(1, 1, 1))
        cuboid5 = Cuboid(origin=(0, 0, -1), side=(1, 1, 1))
        cuboid6 = Cuboid(origin=(-1, 0, 0), side=(1, 1, 1))
        cuboid7 = Cuboid(origin=(0, -1, 0), side=(1, 1, 1))
        cuboid1.initial()
        cuboid2.initial()
        cuboid3.initial()
        cuboid4.initial()
        cuboid5.initial()
        cuboid6.initial()
        cuboid7.initial()
        cuboid2.clarify_bound_with(cuboid1)
        cuboid3.clarify_bound_with(cuboid1)
        cuboid4.clarify_bound_with(cuboid1)
        cuboid5.clarify_bound_with(cuboid1)
        cuboid6.clarify_bound_with(cuboid1)
        cuboid7.clarify_bound_with(cuboid1)
        up_be_zero = (cuboid2.points[cuboid2.down_index] + cuboid1.points[cuboid1.up_index])[:, 2] / 2 - 1
        right_be_zero = (cuboid3.points[cuboid3.left_index] + cuboid1.points[cuboid1.right_index])[:, 0] / 2 - 1
        back_be_zero = (cuboid4.points[cuboid4.front_index] + cuboid1.points[cuboid1.back_index])[:, 1] / 2 - 1
        down_be_zero = (cuboid5.points[cuboid5.up_index] + cuboid1.points[cuboid1.down_index])[:, 2] / 2
        left_be_zero = (cuboid6.points[cuboid6.right_index] + cuboid1.points[cuboid1.left_index])[:, 0] / 2
        front_be_zero = (cuboid7.points[cuboid7.back_index] + cuboid1.points[cuboid1.front_index])[:, 1] / 2
        self.assertFalse(np.any(up_be_zero))
        self.assertFalse(np.any(right_be_zero))
        self.assertFalse(np.any(back_be_zero))
        self.assertFalse(np.any(down_be_zero))
        self.assertFalse(np.any(left_be_zero))
        self.assertFalse(np.any(front_be_zero))

    def test_rotate_then_clarify_bound_with(self):
        Cuboid.ele_size = 0.2
        cuboid1 = Cuboid(origin=(0, 0, 0), side=(1, 1, 1))
        cuboid2 = Cuboid(origin=(0, 0, 1), side=(1, 1, 1))
        cuboid1.initial()
        cuboid1.rotate(cuboid1.centroid(), direction="pos_xz")
        cuboid2.initial()
        cuboid2.clarify_bound_with(cuboid1)
        be_zero = (cuboid2.points[cuboid2.down_index] + cuboid1.points[cuboid1.up_index])[:, 2] / 2 - 1
        self.assertFalse(np.any(be_zero))


if __name__ == "__main__":
    unittest.main()
