from config import LOGGER
import time


def logging(process: str):
    def _logging(func):
        def __logging(*args, **kwargs):
            LOGGER.info("start {}...".format(process))
            start = time.time()
            res = func(*args, **kwargs)
            end = time.time()
            LOGGER.info("{} finished, elapsed time {:.2f}s".format(process, end - start))
            return res
        return __logging
    return _logging
