#pragma once
#include <sstream>
#include <string>

#include "Marmot/MarmotJournal.h"

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
