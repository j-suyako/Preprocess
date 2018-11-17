"""
这个demo用于画平面
"""
from lib.simple3D import *


def plot_finite_plane(figure=None):
    nor_vor = np.random.random(3)
    intercept = np.random.random()
    points_x_y = np.random.random((10, 2))
    points_z = (-intercept - np.dot(points_x_y, nor_vor[:2])) / nor_vor[-1]
    points = np.column_stack((points_x_y, points_z))

    plane = FinitePlane(points)
    plane.plot(figure)


def plot_infinite_plane(figure=None):
    nor_vor = np.random.random(3)
    intercept = np.random.random()

    plane = InfinitePlane(nor_vor, intercept)
    plane.plot(figure)


if __name__ == '__main__':
    plt.figure(figsize=(20, 10))
    ax1 = plt.subplot2grid((1, 2), (0, 0), projection="3d")
    ax2 = plt.subplot2grid((1, 2), (0, 1), projection="3d")
    plot_finite_plane(ax1)
    plot_infinite_plane(ax2)
    plt.show()
