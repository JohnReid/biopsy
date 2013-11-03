/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/pathway.h"

#include <boost/assign/list_of.hpp>

#include <algorithm>
using namespace std;

BIO_NS_START

TableLinkVec interesting_pathways = boost::assign::list_of
	(TableLink(PATHWAY_DATA,  769)) //wnt
	(TableLink(PATHWAY_DATA,  723)) //VEGF 
	(TableLink(PATHWAY_DATA,  772)) //TNF 
	(TableLink(PATHWAY_DATA,  711)) //TGFbeta
	(TableLink(PATHWAY_DATA,  736)) //PDGF 
	(TableLink(PATHWAY_DATA,  750)) //Insulin 
	(TableLink(PATHWAY_DATA,  722)) //EGF 
	(TableLink(PATHWAY_DATA,  741)) //Epo 
	(TableLink(PATHWAY_DATA,  763)) //c-Kit
	(TableLink(PATHWAY_DATA,  759)) //EDAR 
	(TableLink(PATHWAY_DATA,  771)) //Fas 
	(TableLink(PATHWAY_DATA,  757)) //RANK						 ?????????RANKL??????????
	(TableLink(PATHWAY_DATA,  768)) //APP 
	(TableLink(PATHWAY_DATA, 1004)) //Aurora-A				???????????????Aurora-A cell cycle regulation??????????
	(TableLink(PATHWAY_DATA,  770)) //Beta-catenin 
	(TableLink(PATHWAY_DATA,  767)) //Notch  
	(TableLink(PATHWAY_DATA,  849)) //p38  
	(TableLink(PATHWAY_DATA,  719)) //AhR-signaling / Xenobiotic		???????????????AhR?????????????????
	(TableLink(PATHWAY_DATA,  715)) //HIF-1alpha signaling / Hypoxia-induced pathway  
	(TableLink(PATHWAY_DATA,  869)) //p73  
	(TableLink(PATHWAY_DATA,  712)) //PPAR  
	//(TableLink(PATHWAY_DATA, ??)) //Gs-coupled receptor			??????????
	//(TableLink(PATHWAY_DATA, ??)) //MAPK (Ras/Raf/MEK/ERK)     ????????????????????????????
	;


bool
IsInterestingPathway::operator()(const TableLink & link) const
{
	bool result = interesting_pathways.end() != find(interesting_pathways.begin(), interesting_pathways.end(), link);
	return result;
}

BIO_NS_END
