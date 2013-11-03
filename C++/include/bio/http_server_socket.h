#ifndef BIO_HTTP_SERVER_SOCKET_H_
#define BIO_HTTP_SERVER_SOCKET_H_

#include "bio/defs.h"


#include <HttpdSocket.h>
#include <SocketHandler.h>



BIO_NS_START

class ServerHandler : public SocketHandler
{
public:
	ServerHandler(const std::string& );
	~ServerHandler();

	bool Exists(const std::string& );
	int GetInt(const std::string& );
	std::string GetString(const std::string& );
	bool GetBoolean(const std::string& );
	std::string GetMimetype(const std::string& );

private:
};


struct HttpHeaders
{
	void SetServer(const std::string& x);
	std::string GetServer();
	void SetUserAgent(const std::string& x);
	std::string GetUserAgent();
	void SetHttpReferer(const std::string& x);
	std::string GetHttpReferer();

private:
	std::string m_server;
	std::string m_user_agent;
	std::string m_http_referer;
};

class HttpServerSocket : public HttpdSocket
{
public:
	HttpServerSocket(SocketHandler& );
	~HttpServerSocket();

	virtual void Exec();

};



BIO_NS_END



#endif //BIO_HTTP_SERVER_SOCKET_H_


