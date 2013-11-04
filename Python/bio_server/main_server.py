import biopsy
import bifa_server
import base64, threading , os, socket
from sslUserCheck import CheckUserEngine

from soaplib.wsgi_soap import WSGISoapApp
from soaplib.wsgi_soap import SoapServiceBase
from soaplib.service import soapmethod
from soaplib.client import make_service_client
from soaplib.serializers.primitive import String, Integer, Array, Boolean, Float
from soaplib.serializers.binary import Attachment
from soaplib.serializers.clazz import ClassSerializer

# This does not need to be changed for local Windows testing
LOCALID = 'wsbc.warwick.ac.uk'

from tempfile import mkstemp
import os

global server
global cu
global portNo

def localIp():
    if os.name == 'nt' :
        return '127.0.0.1'
    else :
        name = socket.gethostname()
    if name == 'wsbc.cov.warwick.ac.uk' :
         return 'wsbc.warwick.ac.uk'   #return external ethernet name
    else :
         pt = socket.gethostbyname(name)
         print pt
         return pt

#class userDataClass :
#    username = String
#    password = String
#    OK = Boolean

#    Used as part of the soap interface
class Version(ClassSerializer):
    class types:
        majorVersion = Integer
        minorVersion = Integer


class BiFaWSGISoapApp(WSGISoapApp, SoapServiceBase):
    '''
    This object is a VERY simple extension of the base WSGISoapApp.
    It subclasses both WSGISoapApp, and SoapServiceBase, so that
    an object can simply subclass this single object, and it will
    be both a wsgi application and a soap service.  This is convenient
    if you want to only expose some functionality, and don't need
    complex handler mapping, and all of the functionality can be put
    in a single class.
    '''


    def onWsdl(self, environ, wsdl):
        client = make_service_client('http://%s' % (localId), BiFa())
        return client.server.wsdl('')
        '''
        This is called when a wsdl is requested
        @param the wsgi environment
        @param the wsdl string
        '''

    def __init__(self):
        self.cookie = ""
        self.state = -9

        WSGISoapApp.__init__(self)
        SoapServiceBase.__init__(self)

    def getHandler(self, environ):
        global userCheckedEvent
        global checkUserEvent
        auth = environ.get("HTTP_AUTHORIZATION")
        if auth == None :
            raise Exception("Requests must include HTTP authorization")
        if auth == '' :
            raise Exception("Requests must include HTTP authorization")
        if auth[0:6]=="Basic " :
            auth = auth[6:]
        else :
            raise Exception("Requests must include HTTP basic authorization")
        auth = base64.decodestring(auth)
        user, sep, password = auth.partition(':')

        biopsy.UserAdmin.user = user

        self.cookie,self.state = cu.check(user, password)
        if self.cookie == "" :
            print "No cookie"
            raise Exception("Invalid user")


        return self


#    Soap types that are created in the bifa_server
from bifa_server import BiFaHit,PssmInfoData

class BiFa(BiFaWSGISoapApp):

    @soapmethod(String, Float, Integer, Array(String), Boolean, String, Float, String, Boolean, Array(String), _returns=Array(BiFaHit))
    def BiFaHits(self,sequence,threshold,algorithm,phyloSequences,useConsensusSequences,matrixSpecies,phyloThreshold,matrixNameMatch,useCumulativeLikelihoods,pssmSets):
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        if pssmSets != None :
            pssms = pssmSets[0]
        else :
            pssms = ""
        hits=bifa_server.bifa_hits(sequence, threshold, algorithm, phyloSequences,
                                        useConsensusSequences, matrixSpecies, phyloThreshold, matrixNameMatch,
                                        useCumulativeLikelihoods, pssmSets)
        return hits

    @soapmethod(String, Float, String, Integer, Boolean, Array(String), Boolean, String, Float, String, Boolean, Boolean, Array(String), _returns=String)
    def BiFaAnalyser(self, sequence, threshold, title, algorithm, showLabels, phyloSequences, useConsensusSequences, matrixSpecies, phyloThreshold, matrixNameMatch, useOldAlgorithm, useCumulativeLikelihoods, pssmSets) :
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        if pssmSets != None :
            pssms = pssmSets[0]
        else :
            pssms = ""

        ps = ""
        if phyloSequences != None :
            i = 1
            for seq in phyloSequences :
                 ps += "> Seq %i\n" % i
                 ps += seq
                 ps += "\n"
                 i += 1
        str = "> RefSeq\n%s\n,%s,%f,%i,%i,%i,%s,%f,%s,%i,%i,%s\n" % (sequence, ps, threshold, algorithm, showLabels, useConsensusSequences, matrixSpecies, phyloThreshold, matrixNameMatch, useOldAlgorithm, useCumulativeLikelihoods, pssms )
        cu.log2(str)
        temp=bifa_server.bifa_tool(sequence, threshold, title, algorithm, showLabels, phyloSequences, useConsensusSequences, matrixSpecies, phyloThreshold, matrixNameMatch, useOldAlgorithm, useCumulativeLikelihoods, pssmSets)
        output_svg_file="".join([temp, ".svg"])
    	if os.path.isfile(output_svg_file):
    		f1=open(output_svg_file, 'r')
    		svg_string=f1.readlines()
    		f1.close()
    		os.remove(output_svg_file)
    		return "".join(svg_string)
    	else:
    		return "no file"

    @soapmethod(String, _returns= String)
    def returningString(self, key):
        #    If no password was provided then the cookie is used to validate the client
        #    otherwise it is a real cookie
        p1, sep, p2 = key.partition(':')
        if key == "connection_test" :
            if self.state == -1 :
                rv = "keyCheck:" + self.cookie
            else :
                rv = "established:" + self.cookie
            return rv
        elif key == "connection_info":
            if self.state == -1 :
                rv = "keyCheck:" + self.cookie
            else :
                rv = "established:" + self.cookie + ":"+self.transfacVersion() + "." + self.customPssmVersion()
            return rv
        else :
            return "unknown request"

    @soapmethod(String, String, Float, _returns=Array(Float))
    def scorePssmOnSequence(self, pssm_name, sequence, threshold):
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        return bifa_server.score_pssm_on_sequence(pssm_name, sequence, threshold)

    @soapmethod(Array(String), Array(String), Integer, _returns=Array(String))
    def scorePssmsOnSequences(self, pssmNames, sequences, algorithm):
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        return bifa_server.score_pssms_on_sequences(pssmNames, sequences, algorithm)

    @soapmethod(_returns=Float)
    def bindingPrior(self):
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        return biopsy.Environment.bindingPrior

    @soapmethod(_returns=Integer)
    def maxChainMaxNumSequences(self):
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        return biopsy.Environment.max_chain_max_num_sequences

    @soapmethod(_returns=Array(String))
    def PssmSetNames(self):
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        return bifa_server.get_pssm_set_names()

    @soapmethod(Boolean,String,String,Array(String),_returns=Array(String))
    def Pssms(self,useConsensusSequences,matrixSpecies,matrixNameMatch,pssmSets):
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        return bifa_server.pssmAccs(pssmSets,useConsensusSequences,matrixSpecies,matrixNameMatch)


    @soapmethod(String,_returns=PssmInfoData)
    def PssmInfo(self,pssmName):
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        return bifa_server.get_pssm_info(pssmName)

    @soapmethod(String,_returns=Array(String))
    def PssmFreqs(self,pssmName):
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        if (biopsy.UserAdmin.userType < 2) :
            raise Exception("Operation not allowed for current user")
        return bifa_server.get_pssm_freqs(pssmName)

    @soapmethod(String,_returns=Array(String))
    def PssmCounts(self,pssmName):
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        if (biopsy.UserAdmin.userType < 2) :
            raise Exception("Operation not allowed for current user")
        return bifa_server.get_pssm_counts(pssmName)


    @soapmethod(_returns=Version)
    def serverVersion(self):
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        v = Version()
        v.majorVersion = 5
        v.minorVersion = 1
        return v

    @soapmethod(_returns=String)
    def transfacVersion(self):
        v = "%d.%d" % (biopsy.Environment.transfac_major_version, biopsy.Environment.transfac_minor_version)
        return v

    @soapmethod(_returns=String)
    def customPssmVersion(self):
        v = biopsy.Environment.custom_PSSM_version
        return v

    @soapmethod(_returns=Integer)
    def userType(self):
        ut = biopsy.UserAdmin.userType
        return ut

    @soapmethod(String,_returns=Boolean)
    def setPassword(self,newPassword):
        if not biopsy.UserAdmin.isAllowed :
            raise Exception("Invalid user")
        return biopsy.UserAdmin.setPasswd(newPassword)


class ServerThread(threading.Thread):

    def run(self):
        from cherrypy.wsgiserver import CherryPyWSGIServer
    #    Not possible to use 127.0.0.1 as this does not link to an
    #    externally accessible interface

        global server
        server = CherryPyWSGIServer((localIp(), portNo), BiFa())

        print "Started serving:"
        server.start()

if __name__=='__main__':

#    userData = userDataClass()

#	regenerate the wsdl at every startup.  This is not
#	needed but keeps the code accessible.
#	Server address is a dummy
    client = make_service_client('http://%s' % (localIp()), BiFa())
    f1 = open("bifa.wsdl", 'w')
    f1.write(client.server.wsdl(''))
#	Need to flush before close, otherwise we only get 4096 bytes!
    f1.flush()
    f1.close()

    f2 = file("config.dat")
    for line in f2:
        line = line.strip()
        key = line[ : line.find('=') ].strip()
        value  = line[ line.find('=') + 1 : ].strip()
        if (key == "port"):
            portNo = int(value)
    f2.close()

    print "Port no %i " % portNo

    # start server on a second thread
#    it = ServerThread()
 #   it.start()

    from cherrypy.wsgiserver import CherryPyWSGIServer
    #    Not possible to use 127.0.0.1 as this does not link to an
    #    externally accessible interface

    cu = CheckUserEngine()

#    global server
    server = CherryPyWSGIServer((localIp(), portNo), BiFa())

    print "Started serving:"
    server.start()

    print "Finished"
