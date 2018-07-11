import unittest

from Lib.GeometryLib import *
import numpy as np
from scipy.spatial import Voronoi


class TestGeometry(unittest.TestCase):

    # def __init__(self):
    #     pass

    # def test_norm(self):
    #     a = np.array([3, 4])
    #     self.assertTrue((norm(a) == np.array([0.6, 0.8])).all())
    def test_rotate_z(self):
        test1 = norm(np.random.random(2))
        test1 = np.append(test1, 1)
        theta_of_test1 = np.arccos(test1[0])
        test1_after_rotate = norm(np.random.random(2))
        test1_after_rotate = np.append(test1_after_rotate, 1)
        theta_of_test1_after_rotate = np.arccos(test1_after_rotate[0])
        predict = rotate_z(test1, theta_of_test1_after_rotate - theta_of_test1)
        self.assertTrue((np.abs(predict - test1_after_rotate) < 1e-5).all())

    def test_rotate_x(self):
        test1 = norm(np.random.random(2))
        test1 = np.insert(test1, 0, 1)
        theta_of_test1 = np.arccos(test1[1])
        test1_after_rotate = norm(np.random.random(2))
        test1_after_rotate = np.insert(test1_after_rotate, 0, 1)
        theta_of_test1_after_rotate = np.arccos(test1_after_rotate[1])
        predict = rotate_x(test1, theta_of_test1_after_rotate - theta_of_test1)
        self.assertTrue((np.abs(predict - test1_after_rotate) < 1e-5).all())

    def test_plane_init_less_than_3_points(self):
        points = np.random.random((2, 3))
        self.assertRaisesRegex(ValueError, "at least 3 points are needed", Plane, points)

    def test_plane_init_all_points_in_a_line(self):
        k = np.random.randint(1, 100, 20)
        base = np.random.random((1, 3))
        points = np.dot(k.reshape((-1, 1)), base)
        self.assertRaisesRegex(ValueError, "all points shouldn't be in one line", Plane, points)

    def test_plane_init_points_not_in_same_plane(self):
        points = np.random.random((10, 3))
        self.assertRaisesRegex(ValueError, "all points should be in a same plane", Plane, points)

    # def test_plane_init(self):
    #     nor_vor = norm(np.random.random(3))
    #     z_ = np.random.random()
    #     points_x_y = np.random.random((10, 2))
    #     points_z = (-z_ - np.dot(points_x_y, nor_vor[:2])) / nor_vor[-1]
    #     points = np.column_stack((points_x_y, points_z))
    #     plane = Plane(points)
    #     self.assertTrue((np.abs(plane.normal_vector - nor_vor) < 1e-5).all() or
    #                     (np.abs(plane.normal_vector + nor_vor) < 1e-5).all())
    #     self.assertTrue(np.abs(plane.z_ - z_) < 1e-5 or np.abs(plane.z_ + z_) < 1e-5)

    def test_intersect1(self):
        """检测两直线相交点

        以points[0]为交点分别连接points[1],points[2]得到两条直线，验证其相交点是否为points[0]
        """
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

    def test_intersect2(self):
        """检测直线与平面的相交点

        随机得到一个平面，随机选取平面上的一条线段，再随机选择线段上的一点作为交点inter，求出该交点与
        线段两端交点的距离比t；之后在平面两边各找到一点使其在平面上的投影为线段的两个端点，同时这
        两个点到平面的距离与t相同，验证由这两个点构成的线段与平面的交点是否为inter
        """
        nor_vor = norm(np.random.random(3))
        z_ = np.random.random()
        points_x_y = np.random.random((10, 2))
        points_z = (-z_ - np.dot(points_x_y, nor_vor[:2])) / nor_vor[-1]
        points = np.column_stack((points_x_y, points_z))
        plane = Plane(points)
        point_start = points[0] + np.random.random() * (points[1] - points[0])
        point_end = points[-1] + np.random.random() * (points[-2] - points[-1])
        t1, t2 = np.random.random(2)
        point_inter = point_start + t1 / (t1 + t2) * (point_end - point_start)
        segment_start = point_start + t1 * nor_vor
        segment_end = point_end - t2 * nor_vor
        segment = Segment(segment_start, segment_end)
        self.assertTrue((np.abs(segment.intersect(plane) - point_inter) < 1e-5).all())

    #
    # def test_polyhedron(self):
    #     p = np.random.rand(30, 3)
    #     vor = Voronoi(p)
    #     index = 0
    #     for i, e in enumerate(vor.regions):
    #         if e and np.all(e >= 0):
    #             index = i
    #             break
    #     points = vor.vertices[vor.regions[index]]
    #     poly = Polyhedron(points)
    #     self.assertEqual(poly.min_x, min(p[:]))
    #     self.assertEqual(poly.min_y, 0.5)


if __name__ == '__main__':
    unittest.main()
