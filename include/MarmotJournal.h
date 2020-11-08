#pragma once
#include <iostream>
#include <string>

class MarmotJournal {
  private:
    static MarmotJournal& getInstance();

    std::ostream output;

    MarmotJournal();

  public:
    MarmotJournal( MarmotJournal const& ) = delete;
    void operator=( MarmotJournal const& ) = delete;

    static void setMSGOutputDirection( std::ostream& newOutputStream );

    static bool warningToMSG( const std::string& message );

    static bool notificationToMSG( const std::string& message );
};
