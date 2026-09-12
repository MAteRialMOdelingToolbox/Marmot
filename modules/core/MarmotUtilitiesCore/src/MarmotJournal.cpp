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
  // Flushed, not merely written: otherwise a warning sits in the sink's buffer until the process
  // terminates normally, and an explicit run stopped by a timeout or a scheduler discards every
  // warning it raised. Warnings are rare, so the flush costs nothing.
  getInstance().output << message << std::endl;
  return false;
}

bool MarmotJournal::notificationToMSG( const std::string& message )
{
  // A newline, deliberately NOT a flush: notifications come in bulk -- the substeppers emit one
  // per rejected substep, per quadrature point per increment -- and only WARNINGS need to survive
  // an abnormal termination.
  getInstance().output << message << '\n';
  return true;
}
