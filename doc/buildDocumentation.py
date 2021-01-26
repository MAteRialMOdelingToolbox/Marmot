import os
import argparse

def searchForDocumentedProjects( directory,
                                 motherclass = None):

    modules = []

    for f in os.scandir( directory ):
       
        doc = os.path.join( f.path,   "doc/" ) 
        header = os.path.join( f.path, "include/Marmot/" + f.name + ".h" )
        # check is documentation folder is present
        if os.path.isdir( doc ):
            if motherclass is not None:
                with open( header, "r") as g:    
                    # check if header contains given motherclass
                    if motherclass in g.read():
                        modules.append(  str( f.name ) )        
            else:
                modules.append(  str( f.name ) )
    
    return modules

def createListOfSubpages( modules ):

    markdownString = ""

    for module in modules:
        markdownString += "\n - \subpage " + module.lower().replace("marmot", "")
    
    return markdownString 


if __name__ == "__main__":

    parser = argparse.ArgumentParser( 
            description = "A python tool to build the doxygen documentation automatically" )

    parser.add_argument( "--skipDoxygen",
                         action= "store_true",
                        default = False,
                    )

    args = parser.parse_args()

    cwd = os.getcwd()

    # move to toplevel Marmot directory
    if  str( cwd ).endswith("/doc"):
        os.chdir( ".." )

    """
    Structure of the Documentation

        Content
            |- Continuum Mechanics
            |   |
            |   +- Mechanical Material Models
            |   |   |
            |   |   +- Hypo Elastic Material Models
            |   |   |
            |   |   +- Hyper Elastic Material Models
            |   |   
            |   +- Gradient Enhanced Mechanical Material Models
            |       |
            |       +- Gradient Enhanced Hypo Elastic Material Models
            |    
            |- Finite Elements
            |- Numerical Algorithms
            |   |
            |   +- Substepping Algorithms
            |   |
            |   +-Hughes Winget
            |
            |- Interfacing with Marmot
            |
            |- Others
    """

    coreModules = searchForDocumentedProjects( "modules/core" )
    
    # content page
    with open( "doc/markdown/content.md", "w+") as f:
        f.write( "\page content Content\n" )
        f.write( "The content provided by %Marmot can be found here.\n" )
        f.write( " - \subpage continuummechanics\n" )
        if "MarmotFiniteElementCore" in coreModules:
            f.write( " - \subpage finiteelementtechnology\n") 
        f.write( " - \subpage numericalalgorithms\n" )
        f.write( " - \subpage interfaces\n" )
        f.write( " - \subpage others\n" )
    
    # continuum mechanics page
    with open( "doc/markdown/continuummechanics.md", "w+" ) as f:
        f.write( "\page continuummechanics Continuum Mechanics\n" )
        f.write( "The following basic types of material models are available in %Marmot.\n" )
        f.write( " - \subpage mechanicalmaterials\n" )
        f.write( " - \subpage gradmechanicalmaterials\n" )
        f.write( " - \subpage continuummechanicsothers\n" )
    
    # mechanical materials page
    with open( "doc/markdown/mechanicalmaterials.md", "w+" ) as f:
        f.write( "\page mechanicalmaterials Mechanical Material Models\n" )
        f.write( " - \subpage hypoelastic\n" )
        f.write( " - \subpage hyperelastic\n" )
    
    # gradient enhanced mechanical materials page
    with open( "doc/markdown/gradmechanicalmaterials.md", "w+") as f:
        f.write( "\page gradmechanicalmaterials Gradient Enhanced Mechanical Material Models\n") 
        f.write( " - \subpage gradhypoelastic\n" )

    # collect and categorize materials   
    hypoelasticmaterials = searchForDocumentedProjects( "modules/materials", "MarmotMaterialHypoElastic" )
    gradhypoelasticmaterials = searchForDocumentedProjects( "modules/materials",
            "MarmotMaterialGradientEnhancedHypoElastic" )
    hyperelasticmaterials = searchForDocumentedProjects( "modules/materials", "MarmotMaterialHyperElastic" )


    with open( "doc/markdown/materials.md", "w+" ) as f:
        # hypo elastic materials
        f.write( "\page hypoelastic Hypoelastic Material Models\n" )
        f.write( "The documentation of the available hypoelastic material models in %Marmot can be found here\n")
        f.write( createListOfSubpages( hypoelasticmaterials ) )
    
        # hyper elastic materials
        f.write( "\page hyperelastic Hyperelastic Material Models\n" )
        f.write( "The documentation of the available hyperelastic material models in %Marmot can be found here\n")
        f.write( createListOfSubpages( hyperelasticmaterials ) )
        
        
        # gradient enhanced hypo elastic materials
        f.write( "\page gradhypoelastic Gradient Enhanced Hypoelastic Material Models\n" )
        f.write( "The documentation of the available gradient enhanced hypoelastic material models in %Marmot can be found here\n")
        f.write( createListOfSubpages( gradhypoelasticmaterials ) )


    # numerical algorithms page
    with open( "doc/markdown/numericalalgorithms.md", "w+" ) as f:
        f.write( "\page numericalalgorithms Numerical Algorithms\n" )
        f.write( " - \subpage substepper\n" )
        f.write( " - \subpage hugheswinget\n" )

    # execute doxygen 
    if not args.skipDoxygen:
        try:
            os.system("doxygen doc/config/dconfig")
        except:
            print( "Doxygen execution failed!" ) 
