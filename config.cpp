#include <iostream>
#include <fstream>
#include <string>
//#include <map>    // included in config.h

#include "config.h"

using namespace std;
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
        string errorOpenConfig("Could not open configuration file ");
        errorOpenConfig += filename;
        throw errorOpenConfig;
    }

    cout << filename << " loaded." << endl;
}

bool Config::checkNParameter(const std::string parameterName)
{
    if  ( nParameters_.find(parameterName) == nParameters_.end() )
    {
        string errorBadParameter("Invalid numerical parameter \'");
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
        string errorBadParameter("Invalid string parameter \'");
        errorBadParameter += parameterName;
        errorBadParameter += "\'";
        throw errorBadParameter;
        return false;
    }

    else return true;
}

bool Config::loadParameters()
{
    ifstream fetch(filename_);

    if (!fetch) return false;

    string parameterType(""), parameterName(""), sParameterValue("");    // new empty key and value objects
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

        else if (parameterType  == "string")
        {
            fetch >> sParameterValue;
            sParameters_[parameterName] = sParameterValue;
        }

        else if (parameterType == "//") {
            fetch.ignore( 256, '\n');} // ignore comments

        else
        {   // throw an error for improper type or bad formatting
            string errorBadType(parameterName);
            errorBadType += "has invalid type \'";
            errorBadType += parameterType;
            errorBadType += "\'";
            throw errorBadType;
        }
    }

    return true;
}

void Config::set(int &parameter, const string parameterName)
{
    if (checkNParameter(parameterName) == true) {
    parameter = static_cast<int>(nParameters_.find(parameterName)->second);
        // sets parameter to the value corresponding to the key parameterName in the map nParameters
    }
}

void Config::set(double &parameter, const string parameterName)
{
    if (checkNParameter(parameterName) == true) {
        parameter = static_cast<double>(nParameters_.find(parameterName)->second);
    }
}

void Config::set(string &parameter, const string parameterName)
{
    if (checkSParameter(parameterName) == true) {
        parameter = static_cast<string>(sParameters_.find(parameterName)->second);
    }
}

int Config::get(const int parameter, const string parameterName)
{
    if (checkNParameter(parameterName) == true) {
    return static_cast<int>(nParameters_.find(parameterName)->second);
        // sets parameter to the value corresponding to the key parameterName in the map nParameters
    }

    else return 0;
}

double Config::get(const double parameter, const string parameterName)
{
    if (checkNParameter(parameterName) == true) {
    return static_cast<double>(nParameters_.find(parameterName)->second);
    }

    else return 0;
}

string Config::get(const string parameter, const string parameterName)
{
    if (checkSParameter(parameterName) == true) {
    return static_cast<string>(sParameters_.find(parameterName)->second);
    }

    else return "";
}

std::ostream & config::operator<<(std::ostream &os, const Config &rhs)
{
    for (sMap_t::const_iterator sMapIt = rhs.sParameters_.begin(); sMapIt != rhs.sParameters_.end(); ++sMapIt) {
        os << sMapIt->first << ":  '" << sMapIt->second << "'" << endl;
    }

    for (nMap_t::const_iterator nMapIt = rhs.nParameters_.begin(); nMapIt != rhs.nParameters_.end(); ++nMapIt) {
        os << nMapIt->first << ":  " << nMapIt->second << endl;
    }

    return os;
}

