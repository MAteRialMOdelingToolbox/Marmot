#pragma once
#include <sstream>
#include <string>

// external!
extern "C" bool warningToMSG( const std::string& message );      // must return a 'false';
extern "C" bool notificationToMSG( const std::string& message ); // must return a 'true';

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
