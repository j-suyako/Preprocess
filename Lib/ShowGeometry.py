from Lib.GeometryLib import *
from mpl_toolkits.mplot3d import Axes3D, axes3d, art3d


def plot_segment(figure=None):
    if figure and not isinstance(figure, Axes3D):
        raise TypeError("figure argument must be a 3d layer")
    fig = Axes3D(plt.figure()) if not figure else figure
    point0 = np.array([0, 0, 0])
    point1 = np.array([1, 1, 1])
    seg = Segment(point0, point1)
    seg.plot(fig)


def plot_plane():
    fig = plt.figure()
    for i in range(16):
        ax = fig.add_subplot(4, 4, i + 1, projection='3d')
        nor_vor = np.random.random(3)
        z_ = np.random.random()
        points_x_y = np.random.random((10, 2))
        points_z = (-z_ - np.dot(points_x_y, nor_vor[:2])) / nor_vor[-1]
        points = np.column_stack((points_x_y, points_z))
        # points = np.array([[0.22678432, 0.67052894, -5.705299],
        #                    [0.73405923, 0.95565918, -7.58196795],
        #                    [0.32720054, 0.03831484, -3.71835814],
        #                    [0.74418186, 0.98012121, -7.68370589],
        #                    [0.82108931, 0.55852239, -6.3763387],
        #                    [0.44720405, 0.60459777, -5.87065239],
        #                    [0.2229608, 0.38063816, -4.70572894],
        #                    [0.86857933, 0.44226589, -6.06247031],
        #                    [0.78190689, 0.06631101, -4.62113773],
        #                    [0.76719772, 0.09653725, -4.69855088]])
        plane = Plane(points)
        plane.plot(ax)
        # print(points)


def plot_poly(figure=None):
    # points = np.random.random((8, 3))
    points = np.array([[0, 0, 0],
                       [0, 1, 0],
                       [1, 1, 0],
                       [1, 0, 0],
                       [0.5, 0.5, 0.5]])
    poly = Polyhedron(points)
    poly.plot(figure)


if __name__ == '__main__':
    # fig = plt.figure()
    # ax1 = fig.add_subplot(121, projection='3d')
    # plot_segment(ax1)
    # ax2 = fig.add_subplot(122, projection='3d')
    # plot_poly()
    plot_plane()
    plt.show()
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # X, Y, Z = axes3d.get_test_data(0.05)
    # ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
    # fig = Axes3D(fig)
    # verts = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
    # paths = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    # path = Path(verts, paths)
    # patch = patches.PathPatch(path, alpha=0.5)
    # fig.add_patch(patch)
    # art3d.pathpatch_2d_to_3d(patch, z=[0, 1, 1, 0, 0])
    # plt.show()
