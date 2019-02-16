#include <exception>
#ifndef PARSEEXCEPTION_H
#define PARSEEXCEPTION_H

// TODO: Parse incorrect arguments into expression
class ParseException: public std::exception {
public:
    virtual const char* what() const throw() {
        return "Parsing Error occurred: Wrong number of arguments supplied.";
    }
} parseException;

#endif //PARSEEXCEPTION_H
