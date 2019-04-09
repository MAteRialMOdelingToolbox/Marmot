#include "bftJournal.h"

BftJournal& BftJournal::getInstance()
{
    static BftJournal instance; // Guaranteed to be destroyed.
                                // Instantiated on first use.
    return instance;
}

BftJournal::BftJournal() : output( nullptr ) {}

void BftJournal::setMSGOutputDirection( std::ostream& newOutputStream )
{
    getInstance().output.rdbuf( newOutputStream.rdbuf() );
}

bool BftJournal::warningToMSG( const std::string& message )
{
    getInstance().output << message;
    return false;
}

bool BftJournal::notificationToMSG( const std::string& message )
{
    getInstance().output << message;
    return true;
}
