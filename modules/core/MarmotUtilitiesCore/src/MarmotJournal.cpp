#include "Marmot/MarmotJournal.h"
#include <ostream>

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
  // Flushed, not merely written. Without this a warning sits in the sink's buffer until the
  // process terminates normally -- and an explicit dynamic run is routinely stopped short, by a
  // timeout, a scheduler or a user, at which point every warning it ever raised is discarded
  // unread. Warnings are rare by construction, so the flush costs nothing that matters.
  getInstance().output << message << std::endl;
  return false;
}

bool MarmotJournal::notificationToMSG( const std::string& message )
{
  getInstance().output << message << std::endl;
  return true;
}
