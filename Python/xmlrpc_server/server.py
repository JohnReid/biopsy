#
# Copyright John Reid 2010
#

"""
A XMLRPC server
"""

import ucsc_maf, logging

from SimpleXMLRPCServer import SimpleXMLRPCServer
from SimpleXMLRPCServer import SimpleXMLRPCRequestHandler

logging.basicConfig(level=logging.INFO)


# Restrict to a particular path.
class RequestHandler(SimpleXMLRPCRequestHandler):
    rpc_paths = ('/RPC2',)


# Create server
logging.info('Creating server')
server = SimpleXMLRPCServer(("localhost", 8000), requestHandler=RequestHandler)
server.register_introspection_functions()


# register functions
logging.info('Registering methods')
server.register_instance(ucsc_maf)


# Run the server's main loop
logging.info('Serving')
server.serve_forever()
