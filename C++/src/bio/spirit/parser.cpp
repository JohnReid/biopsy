/**
@file

Copyright John Reid 2006

*/

#include "bio-pch.h"

//#define BOOST_SPIRIT_DEBUG

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/factor.h"

#include <boost/config.hpp>

#if BOOST_VERSION >= 13800
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_closure.hpp>
#include <boost/spirit/include/classic_push_back_actor.hpp>
#include <boost/spirit/include/phoenix1.hpp>
#else
#include <boost/spirit/core.hpp>
#include <boost/spirit/attribute/closure.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>
#include <boost/spirit/phoenix.hpp>
#endif //BOOST_VERSION >= 013800

#include <iostream>

BIO_NS_START

namespace spirit {

const std::string test_factor_table = 
	"VV  TRANSFAC FACTOR TABLE, Release 11.3 - licensed - 2007-09-10, (C) Biobase GmbH\n"
	"XX\n"
	"//\n"
	"AC  T00001\n"
	"XX\n"
	"ID  T00001\n"
	"XX\n"
	"DT  16.09.1996 (created); ewi.\n"
	"CO  Copyright (C), Biobase GmbH.\n"
	"XX\n"
	"FA  AAF\n"
	"XX\n"
	"OS  human, Homo sapiens\n"
	"OC  eukaryota; animalia; metazoa; chordata; vertebrata; tetrapoda; mammalia; eutheria; primates\n"
	"XX\n"
	"SF  similar to GAF;\n"
	"XX\n"
	"FF  induced by interferon-alpha (15-30'), inhibited by 2-AP;\n"
	"XX\n"
	"BS  R02116; AAF$CONS; Quality: 6.\n"
	"BS  R03064; HS$GBP_02; Quality: 6; GBP, G000264; human, Homo sapiens.\n"
	"XX\n"
	"DR  TRANSPATH: MO000026034.\n"
	"XX\n"
	"RN  [1]; RE0000446.\n"
	"RX  PUBMED: 1901265.\n"
	"RA  Decker T., Lew D. J., Mirkowitch J., Darnell J. E.\n"
	"RT  Cytoplasmic activation of GAF, an IFN-gamma-regulated DNA-binding factor\n"
	"RL  EMBO J. 10:927-932 (1991).\n"
	"RN  [2]; RE0001471.\n"
	"RX  PUBMED: 1833631.\n"
	"RA  Decker T., Lew D. J., Darnell J. E.\n"
	"RT  Two distinct alpha-interferon-dependent signal transduction pathways may contribute to activation of transcription of the guanylate-binding protein gene\n"
	"RL  Mol. Cell. Biol. 11:5147-5153 (1991).\n"
	"XX\n"
	"//\n"
	"AC  T00002\n"
	"XX\n"
	"ID  T00002\n"
	"XX\n"
	"DT  28.10.1992 (created); ewi.\n"
	"DT  02.10.2002 (updated); hom.\n"
	"CO  Copyright (C), Biobase GmbH.\n"
	"XX\n"
	"FA  ACE2\n"
	"XX\n"
	"SY  ACE2; YLR131C.\n"
	"XX\n"
	"OS  yeast, Saccharomyces cerevisiae\n"
	"OC  Eukaryota; Fungi; Ascomycota; Hemiascomycetes; Saccharomycetales; Saccharomycetaceae; Saccharomyces.\n"
	"XX\n"
	"GE  G004153; ACE2.\n"
	"XX\n"
	"HO  SWI5 (S. cerevisiae).\n"
	"XX\n"
	"CL  C0001; CH; 2.3.3.0.1.\n"
	"XX\n"
	"SZ  770 AA; 86.6 kDa (gene) (calc.).\n"
	"XX\n"
	"SQ  MDNVVDPWYINPSGFAKDTQDEEYVQHHDNVNPTIPPPDNYILNNENDDGLDNLLGMDYY\n"
	"SQ  NIDDLLTQELRDLDIPLVPSPKTGDGSSDKKNIDRTWNLGDENNKVSHYSKKSMSSHKRG\n"
	"SQ  LSGTAIFGFLGHNKTLSISSLQQSILNMSKDPQPMELINELGNHNTVKNNNDDFDHIREN\n"
	"SQ  DGENSYLSQVLLKQQEELRIALEKQKEVNEKLEKQLRDNQIQQEKLRKVLEEQEEVAQKL\n"
	"SQ  VSGATNSNSKPGSPVILKTPAMQNGRMKDNAIIVTTNSANGGYQFPPPTLISPRMSNTSI\n"
	"SQ  NGSPSRKYHRQRYPNKSPESNGLNLFSSNSGYLRDSELLSFSPQNYNLNLDGLTYNDHNN\n"
	"SQ  TSDKNNNDKKNSTGDNIFRLFEKTSPGGLSISPRINGNSLRSPFLVGTDKSRDDRYAAGT\n"
	"SQ  FTPRTQLSPIHKKRESVVSTVSTISQLQDDTEPIHMRNTQNPTLRNANALASSSVLPPIP\n"
	"SQ  GSSNNTPIKNSLPQKHVFQHTPVKAPPKNGSNLAPLLNAPDLTDHQLEIKTPIRNNSHCE\n"
	"SQ  VESYPQVPPVTHDIHKSPTLHSTSPLPDEIIPRTTPMKITKKPTTLPPGTIDQYVKELPD\n"
	"SQ  KLFECLYPNCNKVFKRRYNIRSHIQTHLQDRPYSCDFPGCTKAFVRNHDLIRHKISHNAK\n"
	"SQ  KYICPCGKRFNREDALMVHRSRMICTGGKKLEHSINKKLTSPKKSLLDSPHDTSPVKETI\n"
	"SQ  ARDKDGSVLMKMEEQLRDDMRKHGLLDPPPSTAAHEQNSNRTLSNETDAL\n"
	"XX\n"
	"SC  SwissProt #P21192\n"
	"XX\n"
	"FT        1    284   PF00478; IMP dehydrogenase / GMP reductase domain.\n"
	"FT      603    627   PF00096; zf-C2H2.\n"
	"FT      603    627   SM00355; c2h2final6.\n"
	"FT      603    632   PS50157; ZINC_FINGER_C2H2_2.\n"
	"FT      605    627   zinc finger 1.\n"
	"FT      633    657   PF00096; zf-C2H2.\n"
	"FT      633    657   SM00355; c2h2final6.\n"
	"FT      633    662   PS50157; ZINC_FINGER_C2H2_2.\n"
	"FT      635    657   zinc finger 2.\n"
	"FT      664    683   zinc finger 3.\n"
	"XX\n"
	"SF  3 zinc finger motifs;\n"
	"XX\n"
	"FF  cell cycle-regulated activator of metallothionein and chitinase genes;\n"
	"FF  localisation is regulated via phosphorylation state [4];\n"
	"FF  daughter cell-specific localization of Ace2p and coordination of Ace2p-dependent transcription with the mitotic exit network is mediated by the Mob2p and Cbk1p protein kinases [5];\n"
	"XX\n"
	"IN  T03200; ASH1; yeast, Saccharomyces cerevisiae.\n"
	"XX\n"
	"BS  R09745; Y$CTS1_01; Quality: 2; CTS1, G002207; yeast, Saccharomyces cerevisiae.\n"
	"BS  R09744; Y$SIC1_02; Quality: 2; SIC1, G002206; yeast, Saccharomyces cerevisiae.\n"
	"XX\n"
	"DR  EMBL: M55619; SCACE2.\n"
	"DR  SWISSPROT: P21192; ACE2_YEAST.\n"
	"DR  PIR: S12943; TWBYA2.\n"
	"XX\n"
	"RN  [1]; RE0000710.\n"
	"RX  PUBMED: 1730413.\n"
	"RA  Dohrmann P. R., Butler G., Tamai K., Dorland S., Greene J. R., Thiele D. J., Stillman B.\n"
	"RT  Parallel pathways of gene regulation: Homologous regulators SWI5 and ACE2 differentially control transcription of HO and chitinase\n"
	"RL  Genes Dev. 6:93-104 (1992).\n"
	"RN  [2]; RE0001513.\n"
	"RX  PUBMED: 1986241.\n"
	"RA  Butler G., Thiele D. J.\n"
	"RT  ACE2, an activator of yeast metallothionein expression which is homologous to SWI5\n"
	"RL  Mol. Cell. Biol. 11:476-485 (1991).\n"
	"RN  [3]; RE0014327.\n"
	"RX  PUBMED: 10409653.\n"
	"RA  McBride H. J., Yu Y., Stillman D. J.\n"
	"RT  Distinct regions of the Swi5 and Ace2 transcription factors are required for specific gene activation\n"
	"RL  J. Biol. Chem. 274:21029-21036 (1999).\n"
	"RN  [4]; RE0015821.\n"
	"RX  PUBMED: 10517323.\n"
	"RA  O'Conallain C., Doolin M. T., Taggart C., Thornton F., Butler G.\n"
	"RT  Regulated nuclear localisation of the yeast transcription factor Ace2p controls expression of chitinase (CTS1) in Saccharomyces cerevisiae\n"
	"RL  Mol. Gen. Genet. 262:275-282 (1999).\n"
	"RN  [5]; RE0018039.\n"
	"RX  PUBMED: 12196508.\n"
	"RA  Weiss E. L., Kurischko C., Zhang C., Shokat K., Drubin D. G., Luca F. C.\n"
	"RT  The Saccharomyces cerevisiae Mob2p-Cbk1p kinase complex promotes polarized growth and acts with the mitotic exit network to facilitate daughter cell-specific localization of Ace2p transcription factor.\n"
	"RL  J. Cell Biol. 158:885-900 (2002).\n"                    
	"XX\n"
	"//\n"
	;

using namespace boost::spirit;
#if BOOST_VERSION >= 103500
 using namespace BOOST_SPIRIT_CLASSIC_NS;
#endif


//
// TransData
//
struct transdata_closure : closure< transdata_closure, TransData > {
    member1 val;
};
struct transdata : public grammar< transdata, transdata_closure::context_t >
{
    template < typename ScannerT >
    struct definition
    {
        definition( transdata const& self )
        {
			_rule =
				str_p("C")[self.val=COMPEL_DATA]
				| str_p("ev")[self.val=EVIDENCE_DATA]
				| str_p("T")[self.val=FACTOR_DATA]
				| str_p("FR")[self.val=FRAGMENT_DATA]
				| str_p("G")[self.val=GENE_DATA]
				| str_p("MO")[self.val=MOLECULE_DATA]
				| str_p("M")[self.val=MATRIX_DATA]
				| str_p("CH")[self.val=PATHWAY_DATA]
				| str_p("XN")[self.val=REACTION_DATA]
				| str_p("RE")[self.val=REFERENCE_DATA]
				| str_p("S")[self.val=S_DATA]
				| str_p("R")[self.val=SITE_DATA]
				| str_p("")[self.val=CELL_DATA]
				;
		}
		rule< ScannerT > _rule;
		rule<ScannerT> const& start() const { return _rule; }
	};
};



//
// TableLink
//
struct tablelink_closure : closure< tablelink_closure, TableLink > {
    member1 val;
};
struct tablelink : public grammar< tablelink, tablelink_closure::context_t >
{
    template < typename ScannerT >
    struct definition
    {
        definition( tablelink const& self )
        {
			_rule = 
				_transdata[ phoenix::bind(&TableLink::table_id)(self.val) = phoenix::arg1 ] 
				>> uint_p[ phoenix::bind(&TableLink::entry_idx)(self.val) = phoenix::arg1 ]
				;
		}
		rule< ScannerT > _rule;
		transdata _transdata;
		rule<ScannerT> const& start() const { return _rule; }
	};
};




//
// FactorLink
//
/*  Commented out because code cannot handle multiple species
struct factorlink_closure : closure< factorlink_closure, FactorLink > {
    member1 val;
};
struct factorlink : public grammar< factorlink, factorlink_closure::context_t >
{
    template < typename ScannerT >
    struct definition
    {
        definition( factorlink const& self )
        {
			//this is what we are trying to parse: "T00526; MyoD; Species: mouse, Mus musculus."
			using namespace phoenix;
			_rule = 
				_tablelink[ bind(&FactorLink::link)(self.val) = arg1 ] 
				>> "; "
				>> (+(anychar_p - ';'))[ bind(&FactorLink::name)(self.val) = construct_<std::string>(arg1, arg2) ]
				>> "; "
				>> "Species: "
				>> (+(anychar_p - '.'))[ bind(&FactorLink::species[0])(self.val) = construct_<std::string>(arg1, arg2) ]
				>> '.'
				;
		}
		rule< ScannerT > _rule;
		tablelink _tablelink;
		rule<ScannerT> const& start() const { return _rule; }
	};
};
*/



} //namespace spirit
BIO_NS_END

int main()
{
	using namespace BIO_NS;
	using namespace BIO_NS::spirit;

	TransData parsed_transdata;
	//parse("MA", transdata()[assign_a(parsed_transdata)]);
	
	TableLink parsed_tablelink;
	//parse("M001", tablelink()[assign_a(parsed_tablelink)]);
	
	FactorLink parsed_factorlink;
	//parse("T00526; MyoD; Species: mouse, Mus musculus.", factorlink()[assign_a(parsed_factorlink)]);

	return 0;
}

