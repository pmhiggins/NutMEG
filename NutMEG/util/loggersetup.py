import logging
import logging.config


levels = {'CRITICAL' : logging.CRITICAL,
    'ERROR' : logging.ERROR,
    'WARNING' : logging.WARNING,
    'DEBUG' : logging.DEBUG,
    'INFO' : logging.INFO,
}

class loggersetup():

    def __init__(self, classname, filelevel='INFO', printlevel='WARNING'):

        logging.basicConfig(filename='EcoSysLog.log', filemode='w', level=levels[filelevel])
        self.logger = logging.getLogger(classname)

        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        # create console handler to print to stdout
        ch = logging.StreamHandler()
        ch.setLevel(levels[printlevel])
        ch.setFormatter(formatter)


        fh = logging.FileHandler('EcoSysLog.log')
        fh.setLevel(levels[filelevel])
        fh.setFormatter(formatter)

        self.logger.addHandler(ch)
        self.logger.addHandler(fh)


    @classmethod
    def get_logger(cls, classname, filelevel='INFO', printlevel='WARNING'):
        l = cls(classname, filelevel=filelevel, printlevel=printlevel)
        return l.logger

    #
    # def handler(self, classname):
    #     ch = logging.StreamHandler()
    #     ch.setLevel(levels[level])
    #     ch.setFormatter(logging.Formatter('%(asctime)s - '+classname+' - %(levelname)s - %(message)s'))
    #     self.level=level
    #     self.logger.addHandler(ch)
