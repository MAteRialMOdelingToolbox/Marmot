import os

extensions = [
    "breathe",
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.graphviz",
    "sphinx.ext.mathjax",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.napoleon",
]

project = "Marmot"
copyright = "2025, University of Innsbruck and other authors"

# set sphinx "read the docs" theme
html_theme = "sphinx_rtd_theme"

# set logo
html_logo = "../share/marmot_logo.png"

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
            "epsVE": "\\eps^{\\mathrm{ve}}",
            "epsF": "\\eps^{\\mathrm{f}}",
            "epsDC": "\\eps^{\\mathrm{dc}}",
            "epsSHR": "\\eps^{\\mathrm{shr}}",
            "epsRate": "\\dot{\\eps}",
            "Cel": "\\mathbb{C}",
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
    "MarmotTopLevel": ("../include/Marmot", ["Marmot.h", "MarmotElement.h", "MarmotMaterial.h"]),
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
}

breathe_default_members = ("members", "private-members", "protected-members", "undoc-members")
