#
# Copyright John Reid 2007
#

import biopsy.identifiers.build_map
reload(biopsy.identifiers.build_map)

import biopsy; T = biopsy.transfac

if '__main__' == __name__:


    map, mapper = biopsy.identifiers.build_map.default_map_and_mapper()
    t526 = T.TableLink('T00526')
    mapper(t526.as_db_ref())
    print map.links(t526.as_db_ref())

          #map = biopsy.identifiers.build_map.build_map()
