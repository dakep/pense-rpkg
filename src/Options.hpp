//
//  Options.hpp
//  pense
//
//  Created by David Kepplinger on 2017-05-30.
//  Copyright © 2017 David Kepplinger. All rights reserved.
//

#ifndef Options_hpp
#define Options_hpp

#include <typeinfo>
#include <utility>
#include <map>

class Options
{
private:
    class OptionValue
    {
    private:
        class ValueBase {
        public:
            virtual ValueBase * clone() const = 0;
            virtual ~ValueBase() {}
            virtual const std::type_info& type() const = 0;
        };

        template <typename T>
        class ValueImpl : public ValueBase
        {
        public:
            ValueImpl(const T value) : value(value)
            {}

            const std::type_info& type() const
            {
                return typeid(this->value);
            }

            ValueBase* clone() const
            {
                return new ValueImpl(this->value);
            }

            T operator*() const
            {
                return this->value;
            }
        private:
            const T value;
        };

        ValueBase* value;
    public:
        OptionValue(const OptionValue& other) : value(other.value->clone())
        {}

        template <typename T>
        OptionValue(const T value) : value(new ValueImpl<T>(value))
        {}

        template <typename T>
        T get(T def) const {
            return **dynamic_cast< ValueImpl<T>* >(this->value);
        }

        const std::type_info& type() const
        {
            return this->value->type();
        }

        ~OptionValue()
        {
            delete this->value;
        }
    };

public:
    typedef std::map<std::string, OptionValue> map_type;

    Options()
    {}

    template <typename T>
    Options(const std::string& name, const T& value)
    {
        this->set(name, value);
    }

    template <typename T>
    void set(const std::string& name, const T& value)
    {
        OptionValue optval(value);
        this->optMap.erase(name);
        this->optMap.insert(std::pair<std::string, OptionValue>(name, optval));
    }

    template <typename T>
    T get(const std::string& name, const T& def) const
    {
        map_type::const_iterator optIt = this->optMap.find(name);
        if (optIt != this->optMap.end()) {
            return optIt->second.get(def);
        }

        return def;
    }
private:
    map_type optMap;
};

#endif /* Options_hpp */
