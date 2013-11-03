

#include <bio/http_server_socket.h>
#include <bio/environment.h>
USING_BIO_NS;

#include <SocketHandler.h>
#include <ListenSocket.h>

#include <iostream>
using namespace std;


void
check_http_server()
{
	cout << "******* check_http_server()" << endl;

	SocketHandler h;
	ListenSocket<HttpServerSocket> ll(h);
	if (ll.Bind(BioEnvironment::singleton().http_port))
	{
		throw "Could not bind to http server port";
	}
	h.Add(&ll);
	h.Select(1,0);
	while (h.GetCount())
	{
		h.Select(1,0);
	}
}
