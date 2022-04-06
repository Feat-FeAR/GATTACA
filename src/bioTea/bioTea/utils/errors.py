import logging
from typing import Callable
from abc import ABC, abstractmethod
import traceback

log = logging.getLogger(__name__)


# < Custom BioTea errors > #
class BioTeaError(Exception):
    """Generic class from which all BioTea Errors inherit from."""

    pass


class InvalidGeoId(BioTeaError):
    pass


class UnsupportedChip(BioTeaError):
    pass


class SanityError(BioTeaError):
    pass


class ImageNotFoundError(BioTeaError):
    pass


class InvalidPathError(BioTeaError):
    pass


# --- #


class ErrorHandler(ABC):
    """Abstract class from which error handlers can be specified."""

    error: Exception = None

    # I am not 100% sure this is the correct way to use this, but it fails
    # if "error" is not defined, so I'm happy.
    def __init_subclass__(cls, /, error, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.error = error

    @abstractmethod
    def handle(error: Exception, invoker: Callable, trace: str, inputs: dict) -> None:
        pass


class ErrorManager:
    """Class to run a fuction and handle errors that it raises.

    Since the calling function is passed to the handle, allows re-running the
    fuction again after the error is handled.
    """

    executable = None

    def __init__(self, executable: Callable) -> None:
        self.handlers = dict()
        self.executable = executable

    def add_handler(self, exception: Exception, handler: ErrorHandler) -> None:
        if exception in self.handlers.keys():
            log.warn(f"Overwriting handler for {exception}.")

        self.handlers.update({exception: handler})

    def run(self, *args, **kwargs):
        try:
            self.executable(*args, **kwargs)
        except Exception as e:
            if e in self.handlers.keys():
                target = e
            elif any(
                (
                    compatible_handlers := [
                        e for x in self.handlers.keys if issubclass(x, e)
                    ]
                )
            ):
                if sum(compatible_handlers) > 1:
                    log.warning(
                        f"Multiple equally satisfactory handlers found ({compatible_handlers}). Using the last one added."
                    )
                target = compatible_handlers[-1]

            self.handlers[target].handle(
                error=e,
                invoker=self.executable,
                trace=traceback.format_exc(),
                inputs={"args": args, "kwargs": kwargs},
            )


class RaiserHandler(ErrorHandler, error=BioTeaError):
    """Handle that just raises the error. The 'non-handle'."""

    def handle(error: Exception, invoker: Callable, trace: str, inputs: dict) -> None:
        raise error
