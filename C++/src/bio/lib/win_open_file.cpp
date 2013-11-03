/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"
#include "bio/open_file.h"
#include "bio/log.h"

#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;


#include <windows.h>
#include <shellapi.h>

#include <string>
#include <iostream>

BIO_NS_START

std::string
format_win_error(DWORD dw)
{
	LPVOID lpMsgBuf;

	if (0 == FormatMessage(
		FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		FORMAT_MESSAGE_FROM_SYSTEM,
		NULL,
		dw,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		(LPTSTR) &lpMsgBuf,
		0, NULL ))
	{
		return "Could not format error message";
	}

	const std::string result((const char *) lpMsgBuf);

	LocalFree(lpMsgBuf);

	return result;
}

std::string
get_win_last_error()
{
	return format_win_error(GetLastError());
}


void
open_file(const boost::filesystem::path & file)
{
	const std::string arg = fs::system_complete(file).normalize()._BOOST_FS_NATIVE();

	SHELLEXECUTEINFO ExecuteInfo;

	memset(&ExecuteInfo, 0, sizeof(ExecuteInfo));

	ExecuteInfo.cbSize       = sizeof(ExecuteInfo);
	ExecuteInfo.fMask        = 0;                
	ExecuteInfo.hwnd         = 0;                
	ExecuteInfo.lpVerb       = "open";                      // Operation to perform
	ExecuteInfo.lpFile       = "iexplore.exe"; //"c:\\windows\\notepad.exe";  // Application name
	ExecuteInfo.lpParameters = arg.c_str();           // Additional parameters
	ExecuteInfo.lpDirectory  = 0;                           // Default directory
	ExecuteInfo.nShow        = SW_SHOW;
	ExecuteInfo.hInstApp     = 0;

	if(ShellExecuteEx(&ExecuteInfo) == FALSE)
	{
		// Could not start application -> call 'GetLastError()'
		throw get_win_last_error();
	}

	log_stream()
		<< std::endl
		<< "If you have problems with the information bar in Internet Explorer,"
		<< std::endl
		<< "follow the instructions at:"
		<< std::endl
		<< "http://www.phdcc.com/xpsp2.htm#securityoptions"
		<< std::endl
		<< std::endl;

}

BIO_NS_END

