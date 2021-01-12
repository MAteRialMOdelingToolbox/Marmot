import os
import argparse

def searchForDocumentedModules( directory ):

    modules = []

    for f in os.scandir( directory ):
       
        doc = os.path.join( f.path,   'doc/DOCUMENTATION.md' ) 
        
        # append if doc folder is available
        if os.path.isfile( doc ):
            modules.append(  str( f.name ) )
    
    modules.sort()

    return modules

def createListOfSubpages( modules ):

    markdownString = ''

    for module in modules:
        markdownString += '\n - \subpage ' + module.lower().replace('marmot', '')
    
    return markdownString 


if __name__ == "__main__":

    parser = argparse.ArgumentParser( 
            description = 'A python tool to build the doxygen documentation automatically' )

    parser.add_argument( '--skipDoxygen',
                         action= 'store_true',
                        default = False,
                    )

    args = parser.parse_args()

    cwd = os.getcwd()

    # move to toplevel Marmot directory
    if  str( cwd ).endswith('/doc'):
        os.chdir( '..' )

    
    # modules page
    with open( 'doc/modules.md', 'w+') as f:
        f.write( "\page modules Modules\n" )
        f.write( "The modules provided by %Marmot can be found here.\n" )
        f.write( " - \subpage core\n" )
        f.write( " - \subpage elements\n" )
        f.write( " - \subpage materials\n" )


    # core modules
    coreModules = searchForDocumentedModules( 'modules/core' )
    coreString =  "\page core Core\n"
    coreString += "The documentation of the available core modules in %Marmot can be found here\n"
    coreString += createListOfSubpages( coreModules )
    coreString += "\n### Installation path for all core modules\n"
    coreString += " `Marmot/modules/core/`"
    
    with open( 'doc/core.md', 'w+' ) as f:
        f.write( coreString )

    # materials page   
    materials = searchForDocumentedModules( 'modules/materials' )
    materialsString = "\page materials Materials\n"
    materialsString += "The documentation of the available material models in %Marmot can be found here\n"
    materialsString += createListOfSubpages( materials )
    materialsString += "\n### Installation path for all materials\n"
    materialsString += " `Marmot/modules/materials/`"
    
    with open( 'doc/materials.md', 'w+' ) as f:
        f.write( materialsString )

    # elements page
    elements = searchForDocumentedModules( 'modules/elements' )
    elementsString = "\page elements Elements\n"
    elementsString += "The documentation of the available elements in %Marmot can be found here\n"
    elementsString += createListOfSubpages( elements ) 
    elementsString += "\n### Installation path for all elements\n"
    elementsString += " `Marmot/modules/elements/`"

    with open( 'doc/elements.md', 'w+') as f:
        f.write( elementsString )

    # execute doxygen 
    if not args.skipDoxygen:
        os.system('doxygen doc/dconfig')
