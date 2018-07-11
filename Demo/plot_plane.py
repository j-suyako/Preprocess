"""
这个demo用于画平面
"""
from Lib.GeometryLib import *

if __name__ == '__main__':
    # 生成在同一平面上的点集合
    nor_vor = np.random.random(3)
    z_ = np.random.random()
    points_x_y = np.random.random((10, 2))
    points_z = (-z_ - np.dot(points_x_y, nor_vor[:2])) / nor_vor[-1]
    points = np.column_stack((points_x_y, points_z))

    # 关于这些点生成平面，同时给出这些点生成的凸包
    plane = Plane(points)
    plane.plot()
    plt.show()
