import os
import argparse

def searchForDocumentedModules( directory ):

    modules = []

    for f in os.scandir( directory ):
       
        doc = os.path.join( f.path,   'doc' ) 
        
        # append if doc folder is available
        if os.path.isdir( doc ):
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
        f.write( "You may find the modules available in %Marmot here\n" )
        f.write( " - \subpage core\n" )
        f.write( " - \subpage elements\n" )
        f.write( " - \subpage materials\n" )


    # core modules
    coreModules = searchForDocumentedModules( 'modules/core' )
    coreString =  "\page core Core\n"
    coreString += "The documentation of the available core modules in %Marmot can be found here\n"
    coreString += createListOfSubpages( coreModules )
    
    with open( 'doc/core.md', 'w+' ) as f:
        f.write( coreString )

    # materials    
    materials = searchForDocumentedModules( 'modules/materials' )
    materialsString = "\page materials Materials\n"
    materialsString += "The documentation of the available material models in %Marmot can be found here\n"
    materialsString += createListOfSubpages( materials )
    
    with open( 'doc/materials.md', 'w+' ) as f:
        f.write( materialsString )

    # elements
    elements = searchForDocumentedModules( 'modules/elements' )
    elementsString = "\page elements Elements\n"
    elementsString += "The documentation of the available material models in %Marmot can be found here\n"
    elementsString += createListOfSubpages( elements )
    
    with open( 'doc/elements.md', 'w+') as f:
        f.write( elementsString )

    # execute doxygen 
    if not args.skipDoxygen:
        os.system('doxygen doc/dconfig')
