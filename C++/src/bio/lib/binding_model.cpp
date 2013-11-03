/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/binding_model.h"
#include "bio/biobase_binding_model.h"

BIO_NS_START




BindingModel::~BindingModel()
{
}



BindingModel::parameter_t::~parameter_t()
{
}

std::ostream &
operator<<( std::ostream & os, const BindingModel * model )
{
	static const std::string no_model( "<no model>" );

	return os << (model ? model->get_name() : no_model);
}


BIO_NS_END

