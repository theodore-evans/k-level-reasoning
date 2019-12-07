#include <iostream>
#include <fstream>
#include <string>
//#include <map>    // included in config.h

#include "config.h"

using namespace config;

Config::~Config() {
    nParameters_.clear();   // Clear map objects
    sParameters_.clear();   //
}

Config::Config(const char* filename) : filename_(filename)
{
    bool loadSuccess = loadParameters();
    if (!loadSuccess)
    {
        std::string errorOpenConfig("Could not open configuration file ");
        errorOpenConfig += filename;
        throw errorOpenConfig;
    }

    std::cout << filename << " loaded." << std::endl;
}

bool Config::checkNParameter(const std::string parameterName)
{
    if  ( nParameters_.find(parameterName) == nParameters_.end() )
    {
        std::string errorBadParameter("Invalid numerical parameter \'");
        errorBadParameter += parameterName;
        errorBadParameter += "\'";
        throw errorBadParameter;
        return false;
    }

    else return true;
}

bool Config::checkSParameter(const std::string parameterName)
{
    if  ( sParameters_.find(parameterName) == sParameters_.end() )
    {
        std::string errorBadParameter("Invalid std::string parameter \'");
        errorBadParameter += parameterName;
        errorBadParameter += "\'";
        throw errorBadParameter;
        return false;
    }

    else return true;
}

bool Config::loadParameters()
{
    std::ifstream fetch(filename_);

    if (!fetch) return false;

    std::string parameterType(""), parameterName(""), sParameterValue("");    // new empty key and value objects
    double nParameterValue(0);

    while ( fetch.good() )
    {
        fetch >> parameterType;     // each line has format:
        fetch >> parameterName;     // [type] [name] [value]

        if  (parameterType  == "int" || parameterType == "double")
        {
            fetch >> nParameterValue;
            nParameters_[parameterName] = nParameterValue;
        }

        else if (parameterType  == "std::string")
        {
            fetch >> sParameterValue;
            sParameters_[parameterName] = sParameterValue;
        }

        else if (parameterType == "//") {
            fetch.ignore( 256, '\n');} // ignore comments

        else
        {   // throw an error for improper type or bad formatting
            std::string errorBadType(parameterName);
            errorBadType += "has invalid type \'";
            errorBadType += parameterType;
            errorBadType += "\'";
            throw errorBadType;
        }
    }

    return true;
}

void Config::set(int &parameter, const std::string parameterName)
{
    if (checkNParameter(parameterName) == true) {
    parameter = static_cast<int>(nParameters_.find(parameterName)->second);
        // sets parameter to the value corresponding to the key parameterName in the map nParameters
    }
}

void Config::set(double &parameter, const std::string parameterName)
{
    if (checkNParameter(parameterName) == true) {
        parameter = static_cast<double>(nParameters_.find(parameterName)->second);
    }
}

void Config::set(std::string &parameter, const std::string parameterName)
{
    if (checkSParameter(parameterName) == true) {
        parameter = static_cast<std::string>(sParameters_.find(parameterName)->second);
    }
}

int Config::get(const int parameter, const std::string parameterName)
{
    if (checkNParameter(parameterName) == true) {
    return static_cast<int>(nParameters_.find(parameterName)->second);
        // sets parameter to the value corresponding to the key parameterName in the map nParameters
    }

    else return 0;
}

double Config::get(const double parameter, const std::string parameterName)
{
    if (checkNParameter(parameterName) == true) {
    return static_cast<double>(nParameters_.find(parameterName)->second);
    }

    else return 0;
}

std::string Config::get(const std::string parameter, const std::string parameterName)
{
    if (checkSParameter(parameterName) == true) {
    return static_cast<std::string>(sParameters_.find(parameterName)->second);
    }

    else return "";
}

std::ostream & config::operator<<(std::ostream &os, const Config &rhs)
{
    for (sMap_t::const_iterator sMapIt = rhs.sParameters_.begin(); sMapIt != rhs.sParameters_.end(); ++sMapIt) {
        os << sMapIt->first << ":  '" << sMapIt->second << "'" << std::endl;
    }

    for (nMap_t::const_iterator nMapIt = rhs.nParameters_.begin(); nMapIt != rhs.nParameters_.end(); ++nMapIt) {
        os << nMapIt->first << ":  " << nMapIt->second << std::endl;
    }

    return os;
}
