
#ifndef BIO_USERADMIN_H_
#define BIO_USERADMIN_H_

#define USE_SECURITY

#include "bio/defs.h"

#include <boost/utility.hpp>
#include <boost/thread.hpp>


BIO_NS_START

/**
Stores the identity of the current user on a per thread basis, allowing different 
users to be on different threads, and each being able to set and get the user
identity independently */

class userDataType
{
public:
    userDataType(const char * n,int t = 0):name(n),type(t){};
    std::string name;
    int type;
    bool validated;
    //    Used to validate users, it is returned appended to the username
    std::string cookie;
    std::string passwdHash;
};

class userDataTypePtr
{
    userDataType * udt;
public:
    userDataTypePtr(userDataType * ud) :udt(ud){};
    void set(userDataType * ud){udt = ud;};

    operator userDataType &(){return *udt;};
};

struct user_admin
{
    static userDataType & getUserData(const char * name = 0);
    static bool loadUserDataFile();


    static const char * userGet()
    {
        userDataType & u = getUserData();
        return u.name.c_str();
    }
    static void userSet(const char * name)
    {
        getUserData(name);
    }

    static bool validate(const char * validateString);
    static bool setPasswd(const char * passwd);

    static std::string getCookie();

    static int userType()
    {
        userDataType & u = getUserData();
        if (u.validated)
            return (u.type);
        else
            return 0;
    }

    static bool isAllowed()
    {
        userDataType & u = getUserData();
        return (u.validated && (u.type > 0));
    }

private:
    static boost::thread_specific_ptr<userDataTypePtr> user;

};


BIO_NS_END

#endif //BIO_USERADMIN_H_
