#ifndef PARSEEXCEPTION_H
#define PARSEEXCEPTION_H

#include <exception>

class IllegalArgumentException: public std::exception {
public:
    virtual const char* what() const throw() {
        return "ERROR: Wrong number of arguments supplied. Expected: 8";
    }
} illegalArgumentException;

#endif //PARSEEXCEPTION_H
