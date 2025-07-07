import sys
from pathlib import Path

from asdf_standard import DirectoryResourceMapping

if sys.version_info < (3, 9):
    import importlib_resources
else:
    import importlib.resources as importlib_resources

import asdf_wcs_schemas


def get_resource_mappings():
    """
    Get the resource mapping instances for the datamodel schemas
    and manifests.  This method is registered with the
    asdf.resource_mappings entry point.

    Returns
    -------
    list of collections.abc.Mapping
    """
    resources_root = importlib_resources.files(asdf_wcs_schemas) / "resources"
    if not resources_root.is_dir():
        # In an editable install, the resources directory will exist off the
        # repository root:
        resources_root = Path(__file__).parent.parent.parent / "resources"
        if not resources_root.is_dir():
            raise RuntimeError("Missing resources directory")

    return [
        DirectoryResourceMapping(
            resources_root / "schemas" / "stsci.edu" / "gwcs",
            "http://stsci.edu/schemas/gwcs/",
        ),
        DirectoryResourceMapping(
            resources_root / "manifests",
            "asdf://asdf-format.org/astronomy/gwcs/manifests/",
        ),
    ]
