import os

extensions = [
    "breathe",
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.graphviz",
    "sphinx.ext.mathjax",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
]

project = "Marmot"
copyright = "2026, University of Innsbruck, BOKU Vienna and other authors"

# set sphinx "read the docs" theme
html_theme = "sphinx_rtd_theme"

# set logo

html_logo = "../share/marmot_transparent.png"
html_static_path = ["_static"]
html_css_files = [
    "css/custom.css",
]

# MathJaX configuration
mathjax3_config = {
    "loader": {"load": ["[tex]/configmacros"]},
    "tex": {
        "packages": {"[+]": ["configmacros"]},
        "macros": {
            "sig": "\\boldsymbol{\\sigma}",
            "sigRate": "\\dot{\\sig}",
            "eps": "\\boldsymbol{\\varepsilon}",
            "epsE": "\\eps^{\\mathrm{el}}",
            "epsERate": "\\dot{\\eps}^{\\mathrm{el}}",
            "epsVE": "\\eps^{\\mathrm{ve}}",
            "epsVERate": "\\dot{\\eps}^{\\mathrm{ve}}",
            "epsF": "\\eps^{\\mathrm{f}}",
            "epsFRate": "\\dot{\\eps}^{\\mathrm{f}}",
            "epsDC": "\\eps^{\\mathrm{dc}}",
            "epsDCRate": "\\dot{\\eps}^{\\mathrm{dc}}",
            "epsSHR": "\\eps^{\\mathrm{shr}}",
            "epsRate": "\\dot{\\eps}",
            "Cel": "\\mathbb{C}",
            "CelInv": "\\mathbb{C}^{-1}",
            "DelNu": "\\mathbb{D}_\\nu",
            "PhiRate": "\\dot{\\Phi}",
            "JRate": "\\dot{J}",
            "kl": "\\kappa",  # local field variable
            "knl": "\\bar{\\kl}",  # nonlocal field variable
            "Nk": "\\mathbf{N}^{\\mathrm{\\knl}}",  # shape function vector for the nonlocal field variable
            "fuint": "f^{\\mathrm{u}}_\\mathrm{int,e}",  # element internal force vector associated with the displacement degrees of freedom
            "fext": "f_\\mathrm{ext,e}",  # element external load vector
            "fextb": "f_\\mathrm{ext,e}^\\mathrm{b}",  # body-force vector
            "fextp": "f_\\mathrm{ext,e}^\\mathrm{p}",  # distributed load vector
            "fk": "f^{\\mathrm{\\knl}}_\\mathrm{e}",  # element load vector associated with the nonlocal degrees of freedom
            "qu": "q^{\\mathrm{u}}",  # displacement nodal field vector
            "qk": "q^{\\mathrm{\\knl}}",  # nonlocal nodal field vector
            "ru": "r^{\\mathrm{u}}",  # displacement nodal field vector
            "rk": "r^{\\mathrm{\\knl}}",  # nonlocal nodal field vector
        },
    },
}


def getAllHeadersInFolder(folder):
    """
    Returns a list of all header files in the given folder.
    """
    headers = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith(".h"):
                headers.append(os.path.join(file))
    return headers


# Breathe configuration
breathe_default_project = "Marmot"

breathe_projects_source = {
    "Marmot": ("../", []),
    "MarmotFiniteElementCore": (
        "../modules/core/MarmotFiniteElementCore/include/Marmot",
        getAllHeadersInFolder("../modules/core/MarmotFiniteElementCore/include/Marmot/"),
    ),
    "MarmotFiniteStrainMechanicsCore": (
        "../modules/core/MarmotFiniteStrainMechanicsCore/include/Marmot",
        getAllHeadersInFolder("../modules/core/MarmotFiniteStrainMechanicsCore/include/Marmot/"),
    ),
    "MarmotMathCore": (
        "../modules/core/MarmotMathCore/include/Marmot",
        getAllHeadersInFolder("../modules/core/MarmotMathCore/include/Marmot/"),
    ),
    "MarmotMechanicsCore": (
        "../modules/core/MarmotMechanicsCore/include/Marmot",
        getAllHeadersInFolder("../modules/core/MarmotMechanicsCore/include/Marmot/"),
    ),
    "MarmotUtilitiesCore": (
        "../modules/core/MarmotUtilitiesCore/include/Marmot",
        getAllHeadersInFolder("../modules/core/MarmotUtilitiesCore/include/Marmot/"),
    ),
    "MarmotGradientMechanicsCore": (
        "../modules/core/MarmotGradientMechanicsCore/include/Marmot",
        getAllHeadersInFolder("../modules/core/MarmotGradientMechanicsCore/include/Marmot/"),
    ),
}

breathe_default_members = ("members", "private-members", "protected-members", "undoc-members")

suppress_warnings = [
    # Shared headers (e.g. Marmot.h, MarmotTypedefs.h) are picked up by every
    # per-module autodoxygenindex page and re-emit the same cpp:type directives.
    "duplicate_declaration.cpp",
    # Some Doxygen-generated declarations use nested templates that Sphinx's C++
    # domain parser cannot handle; suppress the resulting parse errors.
    "error.cpp",
    # Duplicate RST target names produced by repeated Doxygen anchors.
    "ref.duplicate",
]

# sphinxcontrib-bibtex configuration
bibtex_bibfiles = ["pages/publications.bib"]
bibtex_default_style = "unsrt"
