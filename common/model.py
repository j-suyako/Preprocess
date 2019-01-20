from common import Cuboid, Block
from concurrent import futures
from functools import reduce
from lib.simple3D import *
from typing import List
from util.logging import logging
import numpy as np
import time


def mirror_projection(points, face: Plane):
    if len(points.shape) == 1:
        t = -(np.dot(face.normal_vector, points) + face.intercept) / (np.dot(face.normal_vector, face.normal_vector))
        return list(points + 2 * t * face.normal_vector)
    else:
        return np.array([mirror_projection(point, face) for point in points])


def output_wrap(e):
    res = ""
    if isinstance(e, int) or isinstance(e, float):
        res += str(e)
    if isinstance(e, list) or isinstance(e, np.ndarray):
        res += str(e).lstrip("[").rstrip("]").strip(" ")
    return res + "\n"


class Model(object):
    hole_origins = np.array([[15, 15, 0], [15, 65, 0],
                             [65, 15, 0], [65, 65, 0],
                             [127.5, 15, 0], [127.5, 65, 0],
                             [190, 15, 0], [190, 65, 0]])
    hole_sides = np.array([[35, 35, 90], [35, 35, 90],
                           [47.5, 35, 90], [47.5, 35, 90],
                           [47.5, 35, 90], [47.5, 35, 90],
                           [35, 35, 90], [35, 35, 90]])

    @staticmethod
    def standard_brick(pan_vector=(0, 0, 0), direction: str=None, ele_size=5):
        holes = [Cuboid(origin=Model.hole_origins[i] + pan_vector, side=Model.hole_sides[i])
                 for i in range(Model.hole_origins.shape[0])]
        brick = Block(bound=Cuboid(origin=pan_vector, side=(240, 115, 90)), holes=holes, ele_size=ele_size)
        brick.initial()
        brick.rotate(brick.centroid(), direction)
        return brick

    def __init__(self, condition: dict):
        vecs = condition["pan"]
        rotate_directions = condition["rotate"]
        ele_size = condition["size"]
        self.layers = len(vecs)
        self.bricks = self._build_bricks(vecs, rotate_directions, ele_size)
        self.mortas, self.entity = self._build_remain([len(vec) for vec in vecs])
        self.ele_size = ele_size
        if "integer-rotate" in condition:
            self.rotate(condition["integer-rotate"])

    @staticmethod
    def _inter(this: Cuboid, that: Cuboid):
        relative_direction = norm(that.centroid() - this.centroid())
        if np.all(np.abs(relative_direction) - 1):
            raise ValueError("Models are not aligning.")
        j = np.argwhere(relative_direction).flatten()[0]
        if np.abs(that.centroid() - this.centroid())[j] < (this.side + that.side)[j] / 2:
            raise ValueError("Models are overlapped.")

        base = [this, that][0 if relative_direction[j] > 0 else 1]
        diff = np.zeros(3)
        diff[j] = base.side[j]
        origin = base.origin + diff
        side = np.array(base.side)
        side[j] = (np.abs(that.centroid() - this.centroid()) - (this.side + that.side) / 2)[j]

        inter = Block(Cuboid(origin=origin, side=side), holes=None, material=Block.morta)
        points = list()

        for p, model in enumerate([this, that]):
            direction = relative_direction * (1 if p == 0 else -1)
            intercepts = [-model.origin - model.side, model.origin]
            i = 1 if sum(direction) == -1 else 0
            touching_face = InfinitePlane(normal_vector=direction, intercept=intercepts[i][j])

            index = [["right_index", "back_index", "up_index"],
                     ["left_index", "front_index", "down_index"]]
            points.append(mirror_projection(model.points[getattr(model, index[i][j])], touching_face))

        inter.initial(np.concatenate(points))
        return inter

    @staticmethod
    def _assemble(this: Cuboid, that: Cuboid, mortas: list):
        morta = Model._inter(this, that)
        mortas.append(morta)
        return this.assemble(morta).assemble(that)

    @staticmethod
    def _build_bricks(vecs, rotate_directions, ele_size):
        bricks = list()
        for layer_vec, layer_direction in zip(vecs, rotate_directions):
            for vec, direction in zip(layer_vec, layer_direction):
                bricks.append(Model.standard_brick(pan_vector=vec, direction=direction, ele_size=ele_size))
        return bricks

    def _build_remain(self, length):
        mortas = list()
        n = 0
        layers = list()
        for l in length:
            layers.append(reduce(lambda x, y: Model._assemble(x, y, mortas), self.bricks[n:n + l]))
            n += l

        if len(layers) == 1:
            entity = layers[0].clone()
        else:
            entity = reduce(lambda x, y: Model._assemble(x, y, mortas), layers)
        return mortas, entity

    def rotate(self, direction):
        centroid = self.entity.centroid()
        for brick in self.bricks:
            brick.rotate(centroid, direction)
        for morta in self.mortas:
            morta.rotate(centroid, direction)
        self.entity.rotate(centroid, direction)
        # return self

    @logging("building")
    def build(self, pin=False):
        with futures.ProcessPoolExecutor() as pool:
            self.bricks = list(pool.map(Block.build, self.bricks))
            self.mortas = list(pool.map(Block.build, self.mortas))
        zlist, index = None, np.argwhere(self.bricks[0].holes[0].side == self.bricks[0].side)[0][0]
        if pin:
            zlist = set()
            # index = np.argmin(self.mortas[-self.layers + 1].side) if len(self.mortas) > 0 else 2
            for morta in self.mortas[-self.layers + 1:]:
                zlist.add(morta.origin[index])
                zlist.add(morta.origin[index] + morta.side[index])
        LOGGER.debug("start removing redundancy...")
        for i, brick in enumerate(self.bricks):
            self.bricks[i] = brick.remove_redundancy(zlist, index)
        LOGGER.debug("removing redundancy finished.")
        for block in list(self.bricks) + list(self.mortas):
            LOGGER.info("{} has {} units".format(block, len(block.units)))

    @logging("plotting")
    def plot(self, path=None, figure=None, **kwargs):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        figure = plt.subplot(111, aspect="equal", projection="3d") if not figure else figure

        figure.axis("off")
        for block in list(self.bricks) + list(self.mortas):
            block.plot(figure=figure, **kwargs)
        # for brick in self.bricks:
        #     brick.plot(figure=figure, **kwargs)
        # for morta in self.mortas:
        #     morta.plot(figure=figure, **kwargs)

        unit_nums = 0
        for block in list(self.bricks) + list(self.mortas):
            unit_nums += len(block.units)
        if path is not None:
            fig_name = "--".join([e1 + str(e2) for e1, e2 in zip(["length-", "weight-", "height-"], self.entity.side)]) + \
                       "--eleNum-" + str(unit_nums) + "--" + \
                       time.strftime('%Y%m%d', time.localtime(time.time())) + ".png"
            LOGGER.info("figure name is {}.".format(fig_name))
            figure.view_init(elev=17, azim=-75)
            plt.savefig(path + fig_name)

    @logging("outputing")
    def output(self, path: str):
        def _write(e):
            if isinstance(e, int) or isinstance(e, float):
                res = str(e)
            elif isinstance(e, list):
                res = str(e).lstrip("[").rstrip("]").strip(" ").replace(",", "")
            elif isinstance(e, np.ndarray):
                res = str(e).lstrip("[").rstrip("]").strip(" ")
            else:
                res = ""
            f.write(res + "\n")

        def _unit_type(unit: Polyhedron):
            unit_type = [0] * 3
            for i, (min, max) in enumerate(zip(self.entity.origin, self.entity.origin + self.entity.side)):
                if np.sum(np.abs(unit.verts[:, i] - min) < 1e-6) >= 3:
                    unit_type[i] = -1
                elif np.sum(np.abs(unit.verts[:, i] - max) < 1e-6) >= 3:
                    unit_type[i] = 1
            return unit_type

        unit_nums = 0
        for block in list(self.bricks) + list(self.mortas):
            unit_nums += len(block.units)
        doc_name = "--".join([e1 + str(e2) for e1, e2 in zip(["length-", "weight-", "height-"], self.entity.side)]) + \
                   "--eleNum-" + str(unit_nums) + "--" + \
                   time.strftime('%Y%m%d', time.localtime(time.time())) + ".txt"
        LOGGER.debug("docname is {}.".format(doc_name))
        with open(path + doc_name, 'w') as f:
            _write(unit_nums)
            _write(self.ele_size)
            _write(self.entity.origin)
            _write(self.entity.side)
            n = 0
            for block in list(self.bricks) + list(self.mortas):
                for unit in block.units:
                    f.write(str(n) + " ")
                    f.write(("1" if unit.material == Block.brick else "0") + " ")
                    _write(_unit_type(unit))
                    for line in unit.to_stream():
                        _write(line)
                    n += 1

    def __repr__(self):
        return str(self.entity)


if __name__ == "__main__":
    with open(r'../output/documents/demo.txt', 'w') as f:
        f.writelines("12")
        a = np.array([1, 2])
        f.writelines(output_wrap(a))
