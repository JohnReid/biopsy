/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/application.h"
#include "bio/options.h"

namespace po = boost::program_options;

#include <xercesc/util/PlatformUtils.hpp>
XERCES_CPP_NAMESPACE_USE

#include <antlr/ANTLRException.hpp>
using namespace antlr;

#include <iostream>
#include <exception>
using namespace std;

#ifdef _WIN32
# include <windows.h>
#else //_WIN32
# include <signal.h>
#endif

BIO_NS_START

namespace detail {

static Application * app;

#ifdef _WIN32
BOOL
CtrlHandler(DWORD fdwCtrlType) 
{
	Application::CtrlSignal signal = Application::CTRL_UNKNOWN_SIGNAL;

	switch (fdwCtrlType)
	{
	case CTRL_C_EVENT:
		cout << endl << "Ctrl-C" << endl << endl;
		signal = Application::CTRL_C_SIGNAL;
		break;

	case CTRL_CLOSE_EVENT:
		cout << endl << "Ctrl-CLOSE" << endl << endl;
		signal = Application::CTRL_CLOSE_SIGNAL;
		break;

	case CTRL_BREAK_EVENT:
		cout << endl << "Ctrl-BREAK" << endl << endl;
		signal = Application::CTRL_BREAK_SIGNAL;
		break;

	case CTRL_LOGOFF_EVENT:
		cout << endl << "Ctrl-LOGOFF" << endl << endl;
		signal = Application::CTRL_LOGOFF_SIGNAL;
		break;

	case CTRL_SHUTDOWN_EVENT:
		cout << endl << "Ctrl-SHUTDOWN" << endl << endl;
		signal = Application::CTRL_SHUTDOWN_SIGNAL;
		break;

	default:
		cout << endl << "Unknown control signal" << endl << endl;
		signal = Application::CTRL_UNKNOWN_SIGNAL;
		break;
	}

	return app->ctrl_handler( signal );
} 
#else //_WIN32
void sighandler(int signum)
{
	printf( "Caught signal = %d\n", signum );
	switch( signum )
	{
	case SIGINT: app->ctrl_handler( Application::CTRL_BREAK_SIGNAL ); break;
	case SIGQUIT: app->ctrl_handler( Application::CTRL_C_SIGNAL ); break;
	case SIGABRT: app->ctrl_handler( Application::CTRL_C_SIGNAL ); break;
	case SIGKILL: app->ctrl_handler( Application::CTRL_C_SIGNAL ); break;
	case SIGTERM: app->ctrl_handler( Application::CTRL_CLOSE_SIGNAL ); break;
	case SIGSTOP: app->ctrl_handler( Application::CTRL_BREAK_SIGNAL ); break;
	default: break;
	}
}

#endif

} //namespace detail

Application::Application(const char * options_name)
: options(options_name)
{
}

bool
Application::ctrl_handler(CtrlSignal signal)
{
	return false;
}

void
Application::register_ctrl_handler()
{
	detail::app = this;
#ifdef _WIN32
	if (! SetConsoleCtrlHandler((PHANDLER_ROUTINE) detail::CtrlHandler, TRUE )) 
	{ 
		throw std::logic_error( "ERROR: Could not set control handler" ); 
	}
#else //_WIN32
	if( SIG_ERR == signal( SIGINT , detail::sighandler ) ) throw std::logic_error( "ERROR: Could not set signal handler" ); 
	//if( SIG_ERR == signal( SIGQUIT, detail::sighandler ) ) throw std::logic_error( "ERROR: Could not set signal handler" ); 
	//if( SIG_ERR == signal( SIGABRT, detail::sighandler ) ) throw std::logic_error( "ERROR: Could not set signal handler" ); 
	//if( SIG_ERR == signal( SIGKILL, detail::sighandler ) ) throw std::logic_error( "ERROR: Could not set signal handler" ); 
	//if( SIG_ERR == signal( SIGTERM, detail::sighandler ) ) throw std::logic_error( "ERROR: Could not set signal handler" ); 
	//if( SIG_ERR == signal( SIGSTOP, detail::sighandler ) ) throw std::logic_error( "ERROR: Could not set signal handler" ); 
#endif
}

void
Application::init()
{
}


int
Application::main(int argc, char * argv [])
{
	int result = 0;

	try
	{
		init();

		BioOptions::singleton().parse(argc, argv, 0, get_options(), get_positional_options());

		if (BioOptions::singleton().help)
		{
			cout << get_options().add(BioOptions::singleton().get_visible_options()) << endl;
			return 0;
		}

		if (BioOptions::singleton().version)
		{
			cout << "Version 0.1" << endl;
			return 0;
		}

		result = task();

	}
	catch (const ANTLRException & ex)
	{
		cerr
			<< to_simple_string( boost::posix_time::second_clock::local_time() )
			<< ": ANTLR Error: "
			<< ex.toString()
			<< endl;
		result = -1;
	}
	catch (const std::exception & ex)
	{
		cerr 
			<< to_simple_string( boost::posix_time::second_clock::local_time() )
			<< ": Error: " 
			<< ex.what() 
			<< endl;
		result = -1;
	}
	catch (const string & msg)
	{
		cerr 
			<< to_simple_string( boost::posix_time::second_clock::local_time() )
			<< ": Error: " 
			<< msg 
			<< endl;
		result = -2;
	}
	catch (const char * msg)
	{
		cerr 
			<< to_simple_string( boost::posix_time::second_clock::local_time() )
			<< ": Error: " 
			<< msg 
			<< endl;
		result = -3;
	}
	catch (...)
	{
		cerr 
			<< to_simple_string( boost::posix_time::second_clock::local_time() )
			<< ": Undefined error" 
			<< endl;
		result = -4;
	}

	return result;
}



po::options_description &
Application::get_options()
{
	return options;
}

po::positional_options_description &
Application::get_positional_options()
{
	return positional_options;
}




BIO_NS_END
