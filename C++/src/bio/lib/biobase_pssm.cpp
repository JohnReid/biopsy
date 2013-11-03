/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/biobase_pssm.h"


BIO_NS_START

TableLink
BiobasePssmKeyTransformer::operator()( TableLink link ) const
{
	return link;
}


TableLink 
BiobasePssmKeyTransformer::operator()( Matrix * matrix ) const
{
	return matrix->get_link();
}



TableLink 
BiobasePssmKeyTransformer::operator()( Matrix::ptr_t matrix ) const
{
	return matrix->get_link();
}



TableLink 
BiobasePssmKeyTransformer::operator()( Matrix::map_t::value_type matrix ) const
{
	return matrix.first;
}



TableLink 
BiobasePssmKeyTransformer::operator()( Site * site ) const
{
	return site->get_link();
}



TableLink 
BiobasePssmKeyTransformer::operator()( Site::ptr_t site ) const
{
	return site->get_link();
}



TableLink 
BiobasePssmKeyTransformer::operator()( Site::map_t::value_type site ) const
{
	return site.first;
}







BIO_NS_END
