import logging
log = logging.getLogger("experiment")

from pathlib import Path
# from lineage import Replicate, Treatment

class Error(RuntimeError):
    pass

class Experiment(object):
    FILENAME = 'experiment.pickle'

    def __init__(self, path, analysis_path=None, readonly=False, overwrite=False):
        # Set up paths
        if not isinstance(path, Path):
            path = Path(path)
        if not path.parent.exists():
            log.error("No parent path: {}".format(path.parent.str()))

        self.path = path

        if analysis_path is None:
            self.analysis_path = self.path
        else:
            if not isinstance(analysis_path, Path):
                analysis_path = Path(analysis_path)
            self.analysis_path = analysis_path

        if not readonly and (overwrite or not self.path.exists()):
            self._new()
        else:
            if not self.exists():
                log.error("No such path {}".format(self.filename))
                raise Error

            self._load()

    @property
    def filepath(self):
        return self.path / self.FILENAME

    @property
    def filename(self):
        return str(self.filepath)

    def exists(self):
        return (self.path / self.FILENAME).exists()

    def add_treatment(self, params):
        pass



    
