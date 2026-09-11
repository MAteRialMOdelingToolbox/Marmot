#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include <sstream>
#include <string>
#include <vector>

using namespace Marmot::Testing;

// ---------------------------------------------------------------------------------------------
// MarmotJournal is a singleton, so setMSGOutputDirection()'s effect on warningToMSG()/
// notificationToMSG() persists across tests within this process -- each test below redirects to
// its own fresh stream first, so it is unaffected by whichever stream a previous test last set.
// ---------------------------------------------------------------------------------------------

void testWarningToMSGWritesToTheRedirectedStreamAndReturnsFalse()
{
  std::ostringstream captured;
  MarmotJournal::setMSGOutputDirection( captured );

  const bool result = MarmotJournal::warningToMSG( "a test warning message" );

  throwExceptionOnFailure( result == false,
                           "warningToMSG() must return false in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( captured.str().find( "a test warning message" ) != std::string::npos,
                           "warningToMSG() did not write its message to the redirected stream in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testNotificationToMSGWritesToTheRedirectedStreamAndReturnsTrue()
{
  std::ostringstream captured;
  MarmotJournal::setMSGOutputDirection( captured );

  const bool result = MarmotJournal::notificationToMSG( "a test notification message" );

  throwExceptionOnFailure( result == true,
                           "notificationToMSG() must return true in " + std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( captured.str().find( "a test notification message" ) != std::string::npos,
                           "notificationToMSG() did not write its message to the redirected stream in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

void testSetMSGOutputDirectionRedirectsAwayFromAPreviousStream()
{
  std::ostringstream first, second;

  MarmotJournal::setMSGOutputDirection( first );
  MarmotJournal::warningToMSG( "goes to first" );

  MarmotJournal::setMSGOutputDirection( second );
  MarmotJournal::warningToMSG( "goes to second" );

  throwExceptionOnFailure( first.str().find( "goes to second" ) == std::string::npos,
                           "setMSGOutputDirection() must redirect away from the previous stream in " +
                             std::string( __PRETTY_FUNCTION__ ) );
  throwExceptionOnFailure( second.str().find( "goes to second" ) != std::string::npos,
                           "setMSGOutputDirection() did not redirect to the new stream in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  const std::vector< std::function< void() > > tests = {
    testWarningToMSGWritesToTheRedirectedStreamAndReturnsFalse,
    testNotificationToMSGWritesToTheRedirectedStreamAndReturnsTrue,
    testSetMSGOutputDirectionRedirectsAwayFromAPreviousStream,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
