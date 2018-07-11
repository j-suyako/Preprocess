"""
这个demo展示如何生成一个多面体
"""
from Lib.GeometryLib import *

if __name__ == '__main__':
    # 在空间中生成任意8个点
    points = np.random.random((10, 3))

    # 由这些点的凸包生成多面体
    poly = Polyhedron(points)
    poly.plot()
    plt.show()
