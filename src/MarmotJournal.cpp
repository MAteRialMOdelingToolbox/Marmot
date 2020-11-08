#include "MarmotJournal.h"

MarmotJournal& MarmotJournal::getInstance()
{
    static MarmotJournal instance; 
                                
    return instance;
}

MarmotJournal::MarmotJournal() : output( nullptr ) {}

void MarmotJournal::setMSGOutputDirection( std::ostream& newOutputStream )
{
    getInstance().output.rdbuf( newOutputStream.rdbuf() );
}

bool MarmotJournal::warningToMSG( const std::string& message )
{
    getInstance().output << message;
    return false;
}

bool MarmotJournal::notificationToMSG( const std::string& message )
{
    getInstance().output << message;
    return true;
}
