from importlib.metadata import version, PackageNotFoundError

from .md_intrinsic_dimension import intrinsic_dimension
from .section_id import section_id
from .secondary_structure_id import secondary_structure_id

# TONI is this list correct?
__all__ = ['md_intrinsic_dimension','section_id', 'secondary_structure_id']


try:
    __version__ = version(__name__)
except PackageNotFoundError:
    # Package is not installed, fallback (useful during local development)
    __version__ = "0.0.0uninstalled"

