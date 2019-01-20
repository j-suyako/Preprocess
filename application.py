from common import Unit, Cuboid, Block, Model
from lib.simple3D import *
from threading import Thread, Semaphore


OPTION = {"standard": {"pan": [[(0, 0, 0)]],
                       "rotate": [[None]],
                       "size": 5},  # one standard brick, no pan, no rotate

          "nine-brick": {"pan": [[(0, 0, 0), (0, 125, 0), (187.5, 62.5, 0)],
                                 [(125, 0, 100), (125, 125, 100), (-62.5, 62.5, 100)],
                                 [(0, 0, 200), (0, 125, 200), (187.5, 62.5, 200)]],
                         "rotate": [[None, None, "pos_xy"],
                                    [None, None, "pos_xy"],
                                    [None, None, "pos_xy"]],
                         "size": 5},  # nine standard brick

          "demo": {"pan": [[(0, 0, 0)],
                           [(0, 0, 100)]],
                   "rotate": [[None],
                              [None]],
                   "size": 5},  # two standard brick, align in z axis

          "vector-nine-brick": {"pan": [[(0, 0, 0), (0, 125, 0), (187.5, 62.5, 0)],
                                        [(125, 0, 100), (125, 125, 100), (-62.5, 62.5, 100)],
                                        [(0, 0, 200), (0, 125, 200), (187.5, 62.5, 200)]],
                                "rotate": [[None, None, "pos_xy"],
                                           [None, None, "pos_xy"],
                                           [None, None, "pos_xy"]],
                                "size": 5,
                                "integer-rotate": "neg_xz"}
          }
IMAGE = True
IMAGE_PATH = r'./output/pictures/'
IMAGE_SHOW = False  # it seems that plt.show() only works in main thread, so it may turn off show temporary
DOCUMENTATION = True
DOCUMENTATION_PATH = r'./output/documents/'


def plot_work(model: Model, sema: Semaphore):
    model.plot(path=IMAGE_PATH)
    sema.release()
    # if IMAGE_SHOW:
    #     LOGGER.debug("start show image in screen, that may take a while")
    #     plt.show()


def output_work(model: Model):
    model.output(path=DOCUMENTATION_PATH)


if __name__ == "__main__":
    LOGGER.info("start preprocessing...")
    model = Model(OPTION["vector-nine-brick"])
    model.build(pin=True)
    sema = Semaphore(0)
    plt.figure(figsize=(20, 15))
    if IMAGE:
        plot_task = Thread(target=plot_work, args=(model, sema))
        plot_task.start()
    if DOCUMENTATION:
        output_task = Thread(target=output_work, args=(model, ))
        output_task.start()
    LOGGER.info("preprocessing finished, plotting and outputing task is on the background, "
                "the final time is depend on the size of model.")
    if IMAGE_SHOW:
        sema.acquire()
        plt.show()
