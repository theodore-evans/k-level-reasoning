#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <map>
#include <string>

namespace config
{
    typedef std::map<std::string, double> nMap_t;
    typedef std::map<std::string, std::string> sMap_t;

    class Config
    {
        friend std::ostream & operator<<(std::ostream &os, const Config &config);

        private:
            const char* filename_;
            nMap_t nParameters_;    // I played around with a few different ways of storing multiple
            sMap_t sParameters_;    // variable types, but they all got too messy. This works fine.

        public:
           ~Config();
            Config() {}
            Config(const char* filename);

            bool loadParameters();
            bool checkNParameter(const std::string parameterName);
            bool checkSParameter(const std::string parameterName);

            void set(int &parameter, const std::string parameterName) ;
            void set(double &parameter, const std::string parameterName) ;
            void set(std::string &parameter, const std::string parameterName) ;

            int get(const int parameter, const std::string parameterName);
            double get(const double parameter, const std::string parameterName);
            std::string get(const std::string parameter, const std::string parameterName);

    };

    std::   ostream & operator<<(std::ostream &os, const Config &config);
};

#endif // CONFIG_H_INCLUDED

