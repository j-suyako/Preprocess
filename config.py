import logging
import os
import sys


class Logger(object):
    def __init__(self, name: str, filename: str, level=logging.DEBUG):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(level=level)
        formatter = logging.Formatter("%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s",
                                      datefmt="%d %b %Y %H:%M:%S")
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level=level)
        console_handler.setFormatter(formatter)

        file_handler = logging.FileHandler(os.path.join(os.path.dirname(__file__), "{}.log".format(filename)), mode="w")
        file_handler.setLevel(level=level)
        file_handler.setFormatter(formatter)

        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)

    def getLogger(self):
        return self.logger


LOGGER = Logger("mylog", "application").getLogger()
