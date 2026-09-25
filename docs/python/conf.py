import sys
from pathlib import Path

docs_dir = Path(__file__).resolve().parent
package_dir = docs_dir.parent.parent / "highspy"
sys.path.insert(0, str(package_dir))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "highspy"
copyright = "HiGHS developers"
author = "HiGHS developers"
release = "1.15.1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_markdown_builder",
    "myst_parser",
]

# All source pages are authored directly in Markdown (via MyST) rather than
# reStructuredText.
source_suffix = {
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "_doctrees", "Thumbs.db", ".DS_Store", "site", "docs"]

# We only ever use `autosummary` without `:toctree:` (plain summary tables
# linking back to entries already documented via `automodule`), so no stub
# pages need to be generated.
autosummary_generate = False
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

try:
    import highspy._core  # noqa: F401
except ImportError as exc:
    raise RuntimeError(
        "Could not import the compiled highspy._core extension. Build/"
        "install a real highspy before building these docs, e.g. "
        "`pip install ./highspy` from the repository root."
    ) from exc

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True
markdown_anchor_sections = True

# Give each documented class/method its own stable anchor (e.g.
# highspy.Highs.setOptionValue) so other docs (including the Julia
# Documenter site) can link directly to a specific API entry.
markdown_anchor_signatures = True
