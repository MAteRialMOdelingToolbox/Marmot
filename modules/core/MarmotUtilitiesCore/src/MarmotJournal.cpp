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
  // A newline, but deliberately NOT a flush. Unlike a warning, a notification is raised in bulk:
  // the substeppers emit one per rejected substep, which is per quadrature point per increment,
  // and now that a consumer can actually point this stream somewhere the flush would be paid on
  // every one of them. What the reliability fix needed was for WARNINGS to survive an abnormal
  // termination; notifications are progress chatter and may sit in the buffer.
  getInstance().output << message << '\n';
  return true;
}
