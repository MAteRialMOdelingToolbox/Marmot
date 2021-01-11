import os

def searchForDocumentedModules( directory ):

    modules = []

    for f in os.scandir( directory ):
       
        doc = os.path.join( f.path,   'doc' ) 
        
        # append if doc folder is available
        if os.path.isdir( doc ):
            modules.append(  f.name )

    return modules

def createListOfSubpages( modules ):

    markdownString = ''

    for module in modules:

        markdownString += '\n - \subpage ' + module.lower().replace('marmot', '')
    
    return markdownString

def writeSubpagesToFile( string, template, outputfile ):
    

    with open( template, 'r' ) as temp:

        data = temp.read().replace( 'PLACEHOLDER_FOR_SUBPAGELIST', string )
        
        with open( outputfile, 'w+') as f:
            f.write( data )
        

cwd = os.getcwd()

# move to toplevel Marmot directory
if  str( cwd ).endswith('/doc'):
    os.chdir( '..' )


# search for (documented) core modules and write core.md
coreModules = searchForDocumentedModules( 'modules/core' )
coreSubpageListString = createListOfSubpages( coreModules )
writeSubpagesToFile( coreSubpageListString, 'doc/coreTemplate.md', 'doc/core.md' )

# search for (documented) materials and write materials.md
materials = searchForDocumentedModules( 'modules/materials' )
materialsSubpageListString = createListOfSubpages( materials )
writeSubpagesToFile( materialsSubpageListString, 'doc/materialsTemplate.md', 'doc/materials.md' )

# search for (documented) elements and write elements.md
elements = searchForDocumentedModules( 'modules/elements' )
elementsSubpageListString = createListOfSubpages( elements )
writeSubpagesToFile( elementsSubpageListString, 'doc/elementsTemplate.md', 'doc/elements.md' )

# execute doxygen 
os.system('doxygen doc/dconfig')
