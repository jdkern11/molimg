import pytest
import logging


def pytest_addoption(parser):
    """Add a command line option to disable logger."""
    parser.addoption(
        "--log-disable",
        action="append",
        default=[],
        help="disable specific loggers",
    )


def pytest_configure(config):
    """Disable the loggers."""
    for name in config.getoption("--log-disable", default=[]):
        logger = logging.getLogger(name)
        logger.propagate = False
