import time, datetime
import biopsy
import sys, os, getpass, threading, socket
import thread
import Crypto
from Crypto.Hash import MD5
import string
from random import Random


# from twisted.conch import error
# from twisted.internet import defer, protocol, reactor, posixbase

# Set to zero for testing so that all user requests are accepted.
CHECKUSER = 1

class CheckUserEngine :
    logfilename = ''

    # Common dictionary of users
    users = {}

    class UserData :
        #    Cached password and status data per user.
        def __init__(self,userType,validPassword,cookie):
            self.userType = userType
            self.validPassword = validPassword
           #    Has the password so we do not have in memory password data
            self.cookies = [];
            self.cookies.append(self.hash(cookie));
            self.retries = 0

        #    Password is hashed for storage and comparison.
        #    No static memory passwords
        def hash(self,string):
            m = MD5.new("fgsudyf67")
            m.update(string)
            return m.hexdigest()

        def cookieMatch(self,string):
            if not self.validPassword :
                return False
            hsh = self.hash(string);
            for i in range(len(self.cookies)) :
                if self.cookies[i] == hsh :
                    return True
            return False

        def setCookie(self,string):
            self.cookies.insert(0,self.hash(string))
            if len(self.cookies) > 5 :
                del self.cookies[5]

#        def setKeyCookie(self,string):
 #           self.cookie = ""
  #          self.keyCookie = string

    def genCookie(self):
        cookie = ''.join( Random().sample(string.letters+string.digits, 12) )
        return cookie

    def __init__(self) :
        #  Initialise per thread data
        self.isValidPassword = 0
        self.isValidGroup = 0
        self.initLog()


    def initLog(self):
        if os.name == 'nt' :
            self.logfilename = 'bifalog.txt'
        else :
            name = socket.gethostname()
            if name == 'wsbc.cov.warwick.ac.uk' :
                self.logfilename = 'bifalog.txt'
            else :
                pt = socket.gethostbyname(name)
                self.logfilename = 'bifalog%s.txt' % pt
        f = open(self.logfilename, 'w')
        f.write("Running")
        f.close()


    def log(self,data):
        f = open(self.logfilename, 'a')
        f.write(data)
        f.write("\n")
        f.close()

    def log2(self,data):
        f = open(self.logfilename, 'a')
#    When running on the mac, prints the current process times so that we can
#    monitor for memory leaks.
#    The windows 'ls' command is simply for testing the feature
        if os.name == 'nt' :
            s = os.popen("dir","r")
        else :
            s = os.popen("ps -l | grep \"/common\"","r")
#    Read the output from the previous command
        a = s.read()
#    and put it in the log.
        f.write(a)
        f.write("\n")
        s.close()
#    followed by the stuff we really wanted to log
        f.write(data)
        f.write("\n")
        f.close()

    def check(self,username,password):
        #    Checks if the username is valid
        now = datetime.datetime.now()
        self.log("Checking Username = %s   %s" % (username,now.strftime("%Y-%m-%d %H:%M:%S")))

        # First by looking in the local cache
        if username in self.users :
            #    We have had a previous login
            if self.users[username].cookieMatch(password) :
                self.log("Cookie matched %s " % password)
                #  and it matches... We can use existing status
                if self.users[username].userType > 0 :
                    return password,self.users[username].userType
                else :
                    return "",0

        #    Get here only if we do not have a current valid password
        #    Check it using the bifa servers built in password checking
        #    Could be a valid keyCheck hash, or a valid password
#        if biopsy.UserAdmin.validate(password) :
        if biopsy.UserAdmin.validate(password) :
            self.isValidPassword = 1
            self.userType = biopsy.UserAdmin.userType
            cookie = self.genCookie()
        else :
            cookie = biopsy.UserAdmin.getCookie()
            return cookie,-1
        #    persist the data if we have a valid user
        if username in self.users :
           self.users[username].setCookie(cookie)
           self.users[username].userType = self.userType
        else :
           user = self.UserData(self.userType,self.isValidPassword,cookie)
           self.users[username] = user
        return cookie,self.userType


if __name__ == "__main__":
    #    A test program for testing this module
    global cu
    cu = CheckUserEngine()




