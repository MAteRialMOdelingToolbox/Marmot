#pragma once
#include <iostream>
#include <sstream>
#include <string>

class MakeString {
  public:
    std::stringstream stream;
                      operator std::string() const { return stream.str(); }

    template <class T>
    MakeString& operator<<( T const& VAR )
    {
        stream << VAR;
        return *this;
    }
};

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
