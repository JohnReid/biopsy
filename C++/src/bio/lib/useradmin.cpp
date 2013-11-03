/* Copyright Nigel Dyer, John Reid 2009
*/

#include "bio-pch.h"
#include "bio/useradmin.h"
#include "bio/environment.h"

#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/tokenizer.hpp>

namespace fs = boost::filesystem;
using namespace std;

BIO_NS_START

//	The useradmin class has a thread sepecfic variable 'user' which is used
//	to identify the user associated with this particular thread, the bifa
//	server being multi threaded

boost::thread_specific_ptr<userDataTypePtr> user_admin::user;

//	We maintain a map in memory of all the users and their types
typedef std::map<std::string,userDataType> userMapType;


//	Returns the userdata associated with a specific name
//	and if does not exist, male one up
boost::mutex  userDataFile_m;

userMapType & usermap()
{
	static userMapType usermap;
	return usermap;
}

bool user_admin::loadUserDataFile()
{
	userMapType::iterator i;

	for (i = usermap().begin();i != usermap().end();i++)
		(*i).second.type = 0;

	const fs::path path(
		BioEnvironment::singleton().data_dir + DIR_SEP + "info" + DIR_SEP + "config.dat"
	);

	fs::ifstream stream( path );

	if (stream == 0)
		return false;

	string line;
	//	File format

	//		user: nigeldyer,2
	//		user: fred
	//		user: sascha,1
	typedef boost::tokenizer< boost::char_separator< char > > tokenizer_t;
	static boost::char_separator< char > sep( " ,","");

	while (stream && !stream.eof())
	{
		//	Get a line of data from the file and parse it
		getline( stream, line );
		tokenizer_t tokenizer(line, sep );

		tokenizer_t::iterator tok_iter = tokenizer.begin();
		if( tok_iter != tokenizer.end() && (*tok_iter =="user:") )
		{
			tok_iter++;

			if( tok_iter != tokenizer.end())
			{
				string n = *tok_iter++;
				int userType = 1;

				if( tok_iter != tokenizer.end())
				{
					userType = atoi((*tok_iter++).c_str());
				}

				//	If the user already exists then update their
				//	type with the type from the config file
				if ((i = usermap().find(n)) != usermap().end())
					(*i).second.type = userType;
				else
				{
					//	otherwise insert the new user
					i = usermap().insert(userMapType::value_type(n,userDataType(n.c_str(),userType))).first;
				}

				if( tok_iter != tokenizer.end())
				{
					(*i).second.passwdHash = (*tok_iter).c_str();
				}
				else
					(*i).second.passwdHash = "";

			}
		}
	}
	return true;
}

userDataType & user_admin::getUserData(const char * name)
{
	if (name == 0)
	{
		//	The standalone case where there is no username
		//	Create the default user if it does not exist and return
		//	it
		userDataTypePtr * currentUser = user.get();
		if (currentUser == 0)
			return getUserData("");
		else
			return *currentUser;
	}

	userMapType::iterator i = usermap().find(name);

	//	We are about to read from the config file, so only allow one thread
	//	at a time.
	boost::lock_guard<boost::mutex> lock(userDataFile_m);

	//	If we can find the user in the in memory user data then this may
	//	be a new user that has been added to the user config file
	if (i == usermap().end() || ((*i).second.type == 0))
	{
		loadUserDataFile();
		//	Now look for our user in the file
		if ((i = usermap().find(name)) == usermap().end())
		{
			//	Still does not exist, so create one of type 0 (ie invalid)
			//	so that we don't re read the file if he tries again
			usermap().insert(userMapType::value_type(name,userDataType(name,0)));
			i = usermap().find(name);
		}

	}

	//	And finally, get the pointer to the current user
	userDataTypePtr * currentUser = user.get();
	if (currentUser == 0)
		user.reset(new userDataTypePtr(&(*i).second));
	else 
		currentUser -> set(&(*i).second);

	return (*i).second;
}

#define USE_SECURITY
#ifdef USE_SECURITY
#include <openssl/rsa.h>
#include <openssl/engine.h>

#include "time.h"
#include "sys/stat.h"
#include "sys/timeb.h"


const char * baseSF_chars = 
             "0123456789+/"
             "abcdefghijklmnopqrstuvwxyz"
             "ABCDEFGHIJKLMNOPQRSTUVWXYZ";


static inline bool is_baseSF(unsigned char c) {
  return (isalnum(c) || (c == '+') || (c == '/'));
}

#define PUBKEY "a8690C610a7T1Ch0Q03rm9qsAEI33G56CyZmYcBY4szG/byNT/H4thBjAE/2Bw/1ByQWVDgfSlGMVY2zf9lspJ7tdYr0LZa9T2EbpXYsn18rM4MaZbgAly7XzrmXH0jkZe7hPvsWiOfXyj0T8pzNvHZgYzTy+zZ56h63sFrwHN29Fk+u2Wmn0ua1004="

std::string baseSF_decode(std::string const& encoded_string) {
  int in_len = encoded_string.size();
  int i = 0;
  int j = 0;
  int in_ = 0;
  unsigned char char_array_4[4], char_array_3[3];
  std::string ret;

  while (in_len-- && ( encoded_string[in_] != '=') && is_baseSF(encoded_string[in_])) {
    char_array_4[i++] = encoded_string[in_]; in_++;
    if (i ==4) {
      for (i = 0; i <4; i++)
//        char_array_4[i] = base64_chars.find(char_array_4[i]);
          char_array_4[i] = strchr(baseSF_chars,char_array_4[i]) - baseSF_chars;

      char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
      char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
      char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

      for (i = 0; (i < 3); i++)
        ret += char_array_3[i];
      i = 0;
    }
  }

  if (i) {
    for (j = i; j <4; j++)
      char_array_4[j] = 0;

    for (j = 0; j <4; j++)
		char_array_4[j] = strchr(baseSF_chars,char_array_4[j]) - baseSF_chars;

    char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
    char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
    char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

    for (j = 0; (j < i - 1); j++) ret += char_array_3[j];
  }

  return ret;
}

std::string baseSF_encode(unsigned char const* bytes_to_encode, unsigned int in_len) {
  std::string ret;
  int i = 0;
  int j = 0;
  unsigned char char_array_3[3];
  unsigned char char_array_4[4];

  while (in_len--) {
    char_array_3[i++] = *(bytes_to_encode++);
    if (i == 3) {
      char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
      char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
      char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
      char_array_4[3] = char_array_3[2] & 0x3f;

      for(i = 0; (i <4) ; i++)
        ret += baseSF_chars[char_array_4[i]];
      i = 0;
    }
  }

  if (i)
  {
    for(j = i; j < 3; j++)
      char_array_3[j] = '\0';

    char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
    char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
    char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
    char_array_4[3] = char_array_3[2] & 0x3f;

    for (j = 0; (j < i + 1); j++)
      ret += baseSF_chars[char_array_4[j]];
//
    while((i++ < 3))
      ret += '=';

  }

  return ret;

}

string encrypt(const char * text)
{
	unsigned char * encMess = 0;
	char * mess = 0;
	RSA * R = 0;
	string ret;
	try
	{
		std::string key = baseSF_decode(PUBKEY);
		const unsigned char * pk = (const unsigned char *)key.c_str();
		R = d2i_RSAPublicKey(NULL,&pk,key.size());

		size_t s = RSA_size(R);
		encMess = (unsigned char * )malloc(sizeof(unsigned char) * s);

		mess = (char * )malloc(sizeof(char) * s + 10);
		mess[0] = 0;
		size_t p = 0;
		while(p + strlen(text) < s)
		{
			strcat(mess,text);
			p += strlen(text);
		}
		strncat(mess,text,s-p);

		int encMessLen  = RSA_public_encrypt(s,(const unsigned char *) mess,encMess,R,RSA_NO_PADDING);

		ret = baseSF_encode(encMess,encMessLen);
	}
	catch(...)
	{}
	if (R) RSA_free(R);
	free (encMess);
	free (mess);
	return ret;
}

string decrypt(const string & encrypted)
{
	unsigned char * decMess = 0;
	RSA * R = 0;
	string ret;
	try
	{
		std::string key = baseSF_decode(PUBKEY);
		const unsigned char * pk = (const unsigned char *)key.c_str();
		R = d2i_RSAPublicKey(NULL,&pk,key.size());

		int s = RSA_size(R);
		decMess = (unsigned char * )malloc(sizeof(unsigned char) * s);

		std::string enc2 = baseSF_decode(encrypted);
		pk = (const unsigned char *)enc2.c_str();

		int decMessLen  = RSA_public_decrypt(enc2.size(),pk,decMess,R,1);
		if (decMessLen >= 0)
		{
			decMess[decMessLen] = 0;
			ret = (const char *)decMess;
		}
		else
			ret = "";
	}
	catch(...){};

	if (R) RSA_free(R);
	free (decMess);
	return ret;
}



bool user_admin::validate(const char * validateString)
{
	//	Validates 'password style' information that is passed
	//	as part of a request
	userDataType & u = getUserData();
	bool valid = false;

	//	u.cookie is a cookie that we have previously sent so that it can
	//	be encrypted by the client and returned, so that we can verify it
	if (u.cookie != "")
	{
		//	After dycryption, it should be <username>:<cookie>
		string decryptCookie = decrypt(validateString);
		int colonPos = decryptCookie.find(':');
		if (colonPos)
		{
			string user = decryptCookie.substr(0,colonPos);
			string cookie = decryptCookie.substr(colonPos+1);

			if ((cookie == u.cookie) && (user == u.name))
				valid = true;
		}
	}
	if ((!valid) && (validateString) &&  (validateString[0] != 0))
	{
		string passwdHash = encrypt((u.name + validateString).c_str());

		if (u.passwdHash != "")
		{
			//	Alternatively, it may be a conventional password.
			//	The stored encrypted password is prefixed with the username
			//	so that one generated for one user cannot be used for another
			if (passwdHash == u.passwdHash)
				valid = true;
		}
		if (!valid)
		{
			//	Try reloading data
			boost::lock_guard<boost::mutex> lock(userDataFile_m);

			loadUserDataFile();
			if (passwdHash == u.passwdHash)
				valid = true;
		}

	}
	u.validated = valid;
	u.cookie = "";
	return valid;
}

bool user_admin::setPasswd(const char * passwd)
{
	userDataType & u = getUserData();

	if(u.type != 2)
		return false;

	u.passwdHash = encrypt((u.name + passwd).c_str());

    boost::lock_guard<boost::mutex> lock(userDataFile_m);

	string filename(BioEnvironment::singleton().data_dir + DIR_SEP + "info" + DIR_SEP + "config.dat");

	FILE * f = fopen(filename.c_str(),"w");

	if (f == 0)
		return false;
	//	If we can find the user in the in memory user data then this may
	//	be a new user that has been added to the user config file
	for (userMapType::iterator i = usermap().begin(); i != usermap().end();i++)
	{
		userDataType & v = (*i).second;

		if (v.type)
		{
			fprintf(f,"user: %s",v.name.c_str());
			if (v.type > 1)
			{
				fprintf(f,",%i",v.type);
				if (v.passwdHash != "")
					fprintf(f,",%s",v.passwdHash.c_str());
			}
			fprintf(f,"\n");
		}
	}
	fclose(f);

	return true;
}


std::string user_admin::getCookie()
{
	struct timeb st;
	ftime( &st );
	srand( (unsigned int) st.time ^ st.millitm );

	userDataType & u = getUserData();

	int div = RAND_MAX/64;

	char c[21];
	for (int i = 0; i < 20; i++)
		c[i] = baseSF_chars[rand()/div];

	c[20] = 0;
	u.cookie = c;
	return c;
}
#else
bool user_admin::validateCookie(const char * hashedCookie)
{
	return true;
}

std::string user_admin::getCookie()
{

	return "cookie";
}

#endif




BIO_NS_END

