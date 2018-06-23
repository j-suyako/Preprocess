import unittest

from Lib.GeometryLib import norm, Segment, Polyhedron
import numpy as np
from scipy.spatial import Voronoi


class TestGeometry(unittest.TestCase):

    # def __init__(self):
    #     pass

    def test_norm(self):
        a = np.array([3, 4])
        self.assertTrue((norm(a) == np.array([0.6, 0.8])).all())

    # def test_point_init(self):
    #     p = np.random.random(3)
    #     x, y, z = p
    #     with self.assertRaises(ValueError):
    #         Point(p[:2])
    #         Point(x, y)

    def test_intersect(self):
        points = np.random.random((3, 3))
        line1_start = points[1]
        line2_start = points[2]
        vector1 = points[0] - points[1]
        vector2 = points[0] - points[2]
        k1, k2 = 1 + 2 * np.random.random(2)
        line1_end = line1_start + k1 * vector1
        line2_end = line2_start + k2 * vector2
        segment1 = Segment(line1_start, line1_end)
        segment2 = Segment(line2_start, line2_end)
        self.assertTrue((np.abs(Segment.intersect(segment1, segment2) - points[0]) < 1e-5).all())

    def test_polyhedron(self):
        p = np.random.rand(30, 3)
        vor = Voronoi(p)
        index = 0
        for i, e in enumerate(vor.regions):
            if e and np.all(e >= 0):
                index = i
                break
        points = vor.vertices[vor.regions[index]]
        poly = Polyhedron(points)
        self.assertEqual(poly.min_x, min(p[:]))
        self.assertEqual(poly.min_y, 0.5)


if __name__ == '__main__':
    unittest.main()