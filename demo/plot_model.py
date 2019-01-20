from common import block
from lib.simple3D import *
from application import standard_brick


if __name__ == "__main__":
    # brick0 = standard_brick(pan_vector=(10, 10, 10), direction="pos_xy", ele_size=10)
    brick0 = standard_brick(ele_size=10)
    brick0.build()
    ax = plt.subplot(111, aspect="equal", projection="3d")
    plt.axis("off")
    brick0.plot(ax)
    plt.savefig(r"../output/pictures/demo.png")
    plt.show()
