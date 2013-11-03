/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/http_server_socket.h"

#include <cctype>
#include <string>
#include <iostream>
using namespace std;

#include <Parse.h>


BIO_NS_START

/** Converts a string to upper case in place. */
#define BIO_TO_UPPER(s) (std::transform((s).begin(), (s).end(), (s).begin(), (int(*)(int)) toupper), s)

/** Returns the argument converted to upper case. */
std::string bio_to_upper(const std::string & str)
{
	std::string result(str);
	return BIO_TO_UPPER(result);
}



void
HttpHeaders::SetServer(const std::string& x)
{
	m_server = x;
}

std::string
HttpHeaders::GetServer()
{
	Parse pa(m_server,":");
	std::string str = pa.getword();
	return str;
}

void
HttpHeaders::SetUserAgent(const std::string& x)
{
	m_user_agent = x;
}

std::string
HttpHeaders::GetUserAgent()
{
	return m_user_agent;
}

void
HttpHeaders::SetHttpReferer(const std::string& x)
{
	m_http_referer = x;
}

std::string
HttpHeaders::GetHttpReferer()
{
	return m_http_referer;
}



HttpServerSocket::HttpServerSocket(SocketHandler& h)
:HttpdSocket(h)
{
}


HttpServerSocket::~HttpServerSocket()
{
}

void
HttpServerSocket::Exec()
{
	// header
	AddResponseHeader("Date", GetHttpDate());
	AddResponseHeader("Server", "++ 0.001");
	AddResponseHeader("Connection", "close");
	AddResponseHeader("Content-type", "text/html");

	// status
	SetStatus("200");
	SetStatusText("OK");

	// send
	SendResponse();

	// page
	Send("<html>");
	Send("<title>Test</title>");
	Send("<body>");
	Send("<h1>status</h1>");
	//
	Send("</body>");
	Send("</html>");

	// close
	SetCloseAndDelete();
}

BIO_NS_END

