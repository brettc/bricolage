"""
"""

import logging
logging.basicConfig(
    format="%(levelname)-8s | %(message)s",
    level=logging.INFO,
    # datefmt="%H:%M:%S"
)
log = logging.getLogger("")

import pathlib, os

_COMMANDS = {}

class ProgramError(Exception):
    pass

def command(func):
    global _COMMANDS
    nm = func.__name__

    def wrapped(arguments):
        log.info("Beginning command {} ...".format(nm))
        func(arguments)
        log.info("End command {} ...".format(nm))

    _COMMANDS[nm] = wrapped
    return wrapped

def make_path(pathname):
    # The backported version of pathlib doesn't have expanduser!!
    return pathlib.Path(os.path.expanduser(pathname)).absolute()

class Arguments(object):
    """Wrapper for docopt arguments -- makes processing much simpler"""
    class Container(object):
        """Empty class allows referencing without dictionary fluff"""
        pass

    def __init__(self, arguments):
        self.arguments = arguments
        self.command = None
        self.options = Arguments.Container()
        self.arguments = Arguments.Container()

        # Turn arguments of different kinds into sensible stuff
        for arg, value in arguments.items():
            if arg.startswith('--'):
                setattr(self.options, arg[2:], value)
            elif arg.startswith('<'):
                # TODO: Move out sanitisation
                new_arg = arg.replace("<", "").replace(">", "")
                setattr(self.arguments, new_arg, value)
            else:
                if arg not in _COMMANDS:
                    log.warning("Missing command {}!".format(arg))
                if arguments[arg]:
                    self.command = arg

    def pre_run(self):
        pass

    def run(self):
        self.pre_run()
        try:
            if self.command not in _COMMANDS:
                log.error("Cannot execute command {}!".format(self.command))
                return 1
            _COMMANDS[self.command](self)
        except ProgramError:
            log.error("Error: exiting program...")
            return 1
        return 0
