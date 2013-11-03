#ifndef BIO_FILE_H_
#define BIO_FILE_H_

#include "bio/defs.h"

#include <string>



BIO_NS_START


struct filenameify
{
	char operator()(char v) const
	{
		switch (v)
		{
		case '<':
		case '>':
		case ':':
		case '"':
		case '/':
		case '\\':
		case '|':
		case '(':
		case ')':
		case '{':
		case '}':
		case '[':
		case ']':
		case ',':
		case ' ':
			return '_';

		default:
			return v;
		}
	}
};

#define BIO_FILENAMEIFY(x) (std::transform((x).begin(), (x).end(), (x).begin(), filenameify()), (x))

BIO_NS_END

#endif //BIO_FILE_H_



