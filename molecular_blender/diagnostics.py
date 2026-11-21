"""Source locations for reader and chemistry errors, also suitable for Blender reports."""
import logging
import os


class InputError(ValueError):
    def __init__(self, path, section, message, *, line=None, **context):
        self.path = os.fspath(path)
        self.section = section
        self.line = line
        self.context = context
        self.detail = message
        location = self.path if line is None else f'{self.path}:{line}'
        super().__init__(f'{location}: {section}: {message}')


def report_exception(operator, exc):
    """Keep the exception chain in Blender's console and show the cause in its UI."""
    logging.getLogger(__package__).error('%s', exc, exc_info=(type(exc), exc, exc.__traceback__))
    operator.report({'ERROR'}, str(exc) or type(exc).__name__)
