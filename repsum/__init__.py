# Set package info
from .version import __author__
from .version import __email__
from .version import __version__
from .version import __date__
from .version import __copyright__
from .version import __license__

# Set package level imports
#__all__ = ['defaults']
from .defaults import *
import repsum

# Generate a repertoire summarization
def summary():
    repsum.summary()

# Perform repertoire calculations and comparison
def main():
    repsum.main()
