#pragma once 
#include <string>

//external!
extern "C" bool warningToMSG(const std::string& message); // must return a 'false';
extern "C" bool notificationToMSG(const std::string& message); // must return a 'true';
