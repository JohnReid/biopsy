/**
@file

Copyright John Reid 2006, 2013

*/

#ifndef BIOPSY_ALIAS_LOOKUP_H_
#define BIOPSY_ALIAS_LOOKUP_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER



#include <boost/python/str.hpp>
#include <boost/python/list.hpp>
#include "biopsy/defs.h"

namespace biopsy {


boost::python::list lookup(
    boost::python::str from_db_name,
    boost::python::str from_accession,
    boost::python::str to_db_name );

void add(
    boost::python::str from_db_name,
    boost::python::str from_accession,
    boost::python::str to_db_name,
    boost::python::str to_accession );



} //namespace biopsy

#endif //BIOPSY_ALIAS_LOOKUP_H_

