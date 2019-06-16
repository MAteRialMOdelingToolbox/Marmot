#pragma once
#include <iostream>
#include <string>

class BftJournal {
  private:
    static BftJournal& getInstance();

    std::ostream output;

    BftJournal();

  public:
    BftJournal( BftJournal const& ) = delete;
    void operator=( BftJournal const& ) = delete;

    static void setMSGOutputDirection( std::ostream& newOutputStream );

    static bool warningToMSG( const std::string& message );

    static bool notificationToMSG( const std::string& message );
};
