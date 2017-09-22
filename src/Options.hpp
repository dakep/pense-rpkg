//
//  Options.hpp
//  pense
//
//  Created by David Kepplinger on 2017-05-30.
//  Copyright © 2017 David Kepplinger. All rights reserved.
//

#ifndef Options_hpp
#define Options_hpp

#include <utility>
#include <string>
#include <map>

using namespace std;

class Options
{
private:
    class OptionBase
    {
    public:
        virtual ~OptionBase()
        {}
    };

public:
    template<typename T>
    class Option : public OptionBase
    {
    public:
        Option(const T& value) : OptionBase(), value_val(value)
        {}

        const T value() const
        {
            return this->value_val;
        }
    private:
        const T value_val;
    };

    typedef map<std::string, OptionBase*> map_type;

    Options()
    {}

    template<typename T>
    static Options createSimple(const std::string& name, const T& value)
    {
        Options tmp;
        tmp.set(name, value);
        return tmp;
    }

    ~Options()
    {
        map_type::iterator it = this->optMap.begin();
        for(; it != this->optMap.end(); ++it) {
            delete it->second;
        }
    }

    template <typename T>
    void set(const std::string& name, const T& value)
    {
        Option<T>* opt = new Option<T>(value);
        this->optMap.erase(name);
        this->optMap.insert(pair<std::string, OptionBase*>(name, opt));
    }

    template <typename T>
    T get(const std::string& name, const T& def) const
    {
        map_type::const_iterator optIt = this->optMap.find(name);
        if (optIt == this->optMap.end()) {
            return def;
        }

        return (static_cast< const Option<T>* >(optIt->second))->value();
    }
private:
    map_type optMap;
};


#endif /* Options_hpp */
