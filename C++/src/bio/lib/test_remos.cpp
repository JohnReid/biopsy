/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/remo.h"


BIO_NS_START

RemoMap test_remos;

void
build_test_remos() {

	static bool already_built = false;
	if (! already_built) {
		RemoPtr remo;
		
		remo.reset(new Remo);
		remo->map[MOUSE_SPECIES] = 
			"GGAGCGGGTCGCAGGGTGGAGGTGCCCACCACTCTTGGATGGGAGGGCTTCACGTCACTCCGGGTCCTCCCGGCCGGTCCTTCCATATTAGGGCTTCCTGCTTCCCATATATGGCCATGTACGTCACGGCGGAGGCGGGCCCGTGCTGTTCCAGACCCTTGAAATAGAGGCCGATTCGGGGAGTC";
		remo->map[HUMAN_SPECIES] = 
			"CCCTAGGGTGCAGGATGGAGGTGCCGGGCGCTGTCGGATGGGGGGCTTCACGTCACTCCGGGTCCTCCCGGCCGGTCCTGCCATATTAGGGCTTCCTGCTTCCCATATATGGCCATGTACGTCACGACGGAGGCGGACCCGTGCCGTTCCAGACCCTTCAAATAGAGGCGGATCCGGGGAGTCGC";
		remo->map[DOG_SPECIES] = 
			"CGGGGTGCGGGACCGCGGTGCTCCCCGATCTCGGGTAGGGGGCTTCACGTCACTCCGGGTCCTCCCGGCGGGTCTTGCCATATTAGGGCTTCCTGCTTCCCATATATGGCCATGTACGTCACGACGGAGGCGGACCCGTGGCGCTCCAGACCCTTCAAATAGAGGCGGATCCGGGGAGTCGCGAG";
		remo->map[COW_SPECIES] = 
			"GGGATGGCAGTGCCTACAGTTCTTGGAGAGGGGGCTTCACGTCACTCCGGGTCCTCCAGGCCGGTCTTGCCATATTAGGGCTTCCTGCCTCCCATATATGGCCATGTACGTCACGACGGAGGCGGGTCCATGCCGTTCCAGACCCTTCAAATAGAGGCGGATCCGGGGAGTCGCG";
		remo->map[CHICKEN_SPECIES] = 
			"CCCCACCGCCCGGCTTTCCCCGGTGCCAGGCGCTGACGCCAGGCGGGGTCCTCCGGCTCGCCTCCAAATAAGGGCTTCCTGCTTCCCATATATGGCCATGTACGTCAGGGCCGGGGCGGGCTGCGGCCGGTGGGATCCGATATATAGAGGCGGCTCCGTAGCGGGCGAGG";
		test_remos["egr1_1"] = remo;

		remo.reset(new Remo);
		remo->map[MOUSE_SPECIES] = 
			"GGGGGTGTGCGCCGACCCGGAAACGCCATATAAGGAGCAGGAAGGATCCCCCGCCGGAACAGACCTTATTTGGGCAGCGCCTTATATGGAGTGGCCCAATATGGCCCTGCCGCTTCCGGCTCTGGGAGGAGGGGCGAGCGGGGGTTGGGG";
		remo->map[HUMAN_SPECIES] = 
			"GGAGCAACCAGCTGCGACCCGGAAATGCCATATAAGGAGCAGGAAGGATCCCCCGCCGGAACAACCCTTATTTGGGCAGCACCTTATTTGGAGTGGCCCGATATGGCCCGGCCGCTTCCGGCTCTGGGAGGAGGGAAGAAGGCGGAGGGAGGGGC";
		remo->map[DOG_SPECIES] = 
			"GGAGGAGGGAGGGAGCCAGGGAGCAGCCAGCAACGACCCGGAAACGCCATATAAGGAGCAGGAAGGATCCCCCGCCGGAACAACCCTTATTTGGGCAGCGCCTTGTTTGGAACGGCCCGATATGGCCCGGCCGCTTCCGGCTCCGGGAGGAGGGGAGAAGGCGGAGGGAAGAGGCGGGGGCAGTG";
		remo->map[COW_SPECIES] = 
			"GGAGGGAGAGAGCCAGGGAGCAGAAAGCAGAGACCCGGAAACGCCATATAAGGAGCAGGAAGGATCCCCCGCCGGAACAAACCTTATTTGGGCAGCGCCTTATTTGGAGCGCCCCGATAAGGCCTGGCTGCTTCCGGCTCCGGGAGGAGGGAAGAAGGCGGAGGG";
		remo->map[FUGU_SPECIES] = 
			"AGCGGCGATCAGGAAGTTCCAAATAAGGGCAGGGAAATGGAAACCCACCGGAATAGGTCCTTATGTAGTCACGACCTTATAAGGCACGGCGGAGCCCGGCTTTCC";
		test_remos["egr1_2"] = remo;

		remo.reset(new Remo);
		remo->map[MOUSE_SPECIES] = 
			"GGTGGAAGACTTGGACACATAATACTTCAGACTGCTGATTTGCTATTTTCTTAATTCTGGCTCTTCTTTTTCCATATGAATCAAACAGTCTTTTGTCCTTCAGGGTTTGATTTGATCGCTACTGTTCACATTTCACTTTTGTACTTTCCCTCTTCCTATACCATTACTCATTGAGTGCCCTGTGAAATTACAATCGTACATTTTCAACTCAGCAACCCATTTGCAGTACAAAAAATAGGGTCTAAATAATGGCTGAATTAGCCCTACTGGACAGTTTCAGATGTAACACTCTGTAATAATTATATTGCAGGCTGGATTAGGATGCTATTATCATAATCTGGACGTTTACAATTATCTGTAATTTGCAAAGATGCGCCAGGTCTTGATTACAGCAGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTATCACGCTAACCATCACTAAACAGTGACAGTAATAACAGCTAATTTTGCTGGCAATATAAGAGGTGCTGGGGTGTGCAAACAATTTCACACCTGGATGTGCTCACTCAACCAAGAATATAGAGACAGAGCCTCTGCCCTGAGACTCAGAGAAACGCTCTCCTGTGCCT";
		remo->map[RAT_SPECIES] = 
			"TTGAAGAGTTGGACACATAATACTTCAGACCACTGATTTGCTATTTTCTTAATTCTGGTTCTTCTTTTTCCATATGAATCAAACAGTCTTTTGTCCTTCAGGGTTTGATTTGATCGCTACTGTTCACATTTCACTTTTGTATTTTCCTCTTCCTATACTATTACTCATTGAGTGCCCCGTGAAATTACAATCGTACATTTTCAACTCAGCAACCCATTTGCAGTACAAAAAATAGGGTCTAAATAATGGCTGAATTAGCCCTACTGGACAGTTTCAGATGTAACACTCTGTAATAATTATATTGCAGGCTGGATTAGGATGCTATTATCATAATCTGGACGTTTACAATTATCTGTAATTTGCAAAGATGCGCCAGGTCTTGATTACAGCAGCTTTTTTTTTTTTTTGTATCACGCTAACCATCACTAAACAGTGACAGTAATAACAGCTAATTTTGCTGGCAATATAAGAGGTGCTGGGGTGTGCAAACAATTTCACACCTGGATGTGCTCACTCAACCAAGAATATAGAGACAGAGCTTCTGCCCCGAGACTCAGAAAAATACTCTCCTGTGC";
		remo->map[HUMAN_SPECIES] = 
			"CACATCAGCCCTGCACTATTGAAGATTTGGACATGAAATGCTTCAGACCACTGATTTGCTATTTTCTTAATTCTCGCTCTTCTTTTTCCATATGAATCAAACAGTCTTTTGTCCTTCAGGGTTTGATTTGATCGCTACTGTTCACATTTCACTTTTGTATTTTCCTCTTCCTATACCATTACTCATTGAGTGCCCTGTGAAATTACAATCGTACATTTTCAGCTCAGCAACCCATTTGCAGTGCAAAAATAGGGTCTAAATAATGGCTGAATTAGCCCTACTGGACAGTTTCAGATGTAACACTCTGTAATAATTATATTGCAGGCTGGATTAGGATGCTATTATCATAATCTGGACGTTTACAATTATCTGTAATTTGCAAAGATGCGCCAGGTCTTGATTACAGCAGTTTTTTTTTTTTTTTTTTTTTTTTTTGTACCACGCTAACCATCACTAAACAGTGACAGTAATAACAGCTAATTTTGCTGGCAATATAAGAGGTGCTGGGGTGTGCAAACAATTTCACACCTGGATGTGCTCACTCAACCAAGAATATAGAGAAAGAGCTTCTGCCCTGAGACTCAGAAAAATATTCTCCTGTGCTTTGGTTCAGTATAGATTTCTAAACCCTGATCATTGCTTAAGAGATATTCACTGAGGGTAAGTTTTTATTTCTTGCA";
		remo->map[DOG_SPECIES] = 
			"TACAAACACATCAGCCCTACACATATTGAAGATTTGTACACGTAATGCTTCAGACCACTGATTTGCTATTTTCTTAATTCCCGCTCTTCTTTTTCCATATGAATCAAACAGTCTTTTGTCCTTCAGGGTTTGATTTGATCGCTACTGTTCACATTTCACTTTTGTATTTTCCTCTTCCTATACCATTACTCATTGAGTGCCCTGTGAAATTACAATCGTACATTTTCAACTCAGCAACCCATTTGCAGTGCAAAAATAGGGTCTAAATAATGGCTGAATTAGCCCTACTGGACAGTTTCAGATGTAACACTCTGTAATAATTATATTGCAGGCTGGATTAGGATGCTATTATCATAATCTGGACGTTTACAATTATCTGTAATTTGCAAAGATGCGCCAGGTCTTGATTACAGCAGTTTTTTTTTTTTTTTTGTATCACGCTAACCATCACTAAACAGTGACAGTAATAACAGCTAATTTTGCTGGCAATATAAGAGGTGCTGGGGTGTGCAAACAATTTCACACCTGGATGTGCTCACTCAACCAAGAATATAGAGAAAGAGCTTCTGCCCTGAGACTCAGACAAATATTCTCCTCTGCTTCTGTTCAGTATAGATTTCTAGACCCTGATCATTGCTTAAGAGATATTCACTGGTGGTAAGTTTTTATTTCTTCCATACTTCAAGGCACAGTTCATGTTTCTCTTCTCACACTTTCCAA";
		remo->map[XENOPUS_SPECIES] = 
			"AGAATCAATCCTTTGTCCTTTATCGCTCCTCTACCCTTCCCACTTTTGTATTTCCCTCAGCCTAAACTGTCGCTCAATGAGGGCCCAGCGAAATGACAATTGTGCATTTTCCACTCAGCAACCAGTCTGCCGTACAAAATAGGGGCTAAATAATAGTAGAATTAGCCCTACTGGCCGGTTTCAAATGTAACACTCTGTAATAATTATATTGCAGGCTGGATTAGGATGCTATTATCATAATCTGGACGTTTACAATTATCTGTAATTTGCAAAGATGCTCCAGGTCTTGATTACAGCAGTTTTTTTTATCAAGCTAACCATCACTAAACAGTGACAGTAATAACAGCTAATTTTGCTGGCAATATAAGAGGTGCCGGGGTGTGCAAACAATTTCACACCTGGATATGCTCACTCAACCAAGAATGCAAAGAGAGAGTCTGCCCCGACACTCGGGAAGGTATTCTGCCTTTCTCCTGCCCGGCACCCTCCATGCTCCTGCCCGGCA";
		remo->map[TETRAODON_SPECIES] = 
			"CGCCCTGGCTGTGCCATTACATCGAGAGCCCCGTTTATTTCCAACCGAATATTTTCCACTCAGCCAGTTACCCCGTTTTGCAAATCGAGTCTAAATAATGACAGTTAGTTTTACTGGGCAGCTTCAACTGTAACACTCCACAATAATTATTCTCCAGGATTAAGGCCGCCATTATCATAATCTNNNNNNNNNNNN";
		test_remos["dlx5_1"] = remo;

		remo.reset(new Remo);
		remo->map[MOUSE_SPECIES] = 
			"ATACACTCACAGTGGTTTGGCATATATTTGGTGAAATTTTTTAAGGAAAAATTAGTGTTGGTTTCGATATATGGTAGCTTTTTCTCTAACATAATTTGAATAATTCAGCAAAGCCCTACTACCAGCTGTACTTCTGCAGCCTCTTCCATTCTTTCCAGCATTATAATTTTGGTTAATTTTCAATTTTAGGTCCTACGTCTCTGCAATTTGTGTATGAATAACAGAATAATTTCCCTCTTTTGTTTCGCCTTTCCTGTTCCTGAATCTAAATAAAGATGGCTTTTTAGTATTAAAAGTGGAAGAAAATTACAGGTAATTATCTTTGACGGTAAAAACGCTGTAATCAGCGGGCTACATGAAAAATTACTCTAATTATGGCTGCATTTAAGAGAATGGAAAAAAACCTTCTTGTGGATAAAAACCTTAAATTGTCCCCAATGTCTGCTTCAAATTGGATGGCACTGCAGCTGGAGGCTTTGTTCAGAATTGATCCTGGGGAGCTACGAACCCAAAGTTTCACAGTAGGAAGGGGGAAAAAAGAAAAGAAAACATTTTTCCTAATGTAACAATGCGAATGCTAGAAAATGACAAGACTGATCGGTTTTAAACCATTCTGAAGACTGACTGAGCGTGGAAGTTGCTCAACAAAAAAGGGAACGGGGATATTGGTAAGTAGTCTTTGCTTTGTGTCTAATTATAGTGACTGTCCTTTTGCACGACTGTATGTCGCTGTCATTATATGGAGCACTCTGAGTATCTCTATTGACTTCTGATAAATGGCTCAGTGGTTTCAAGGTTCATAATTTGAAGGATGCCCAGGTCTGTAGCCATTAATTCTCACAGTATGGGGCTTCCTGCTTGTATGACGAAAGGACTTTTCATTTTCATTATTTCTTTGCTTTCCACGAGATCGGAAGGATGTGTTTCTTGCAGGATAGAACATTTTCCTACCTTGCTAGCTTGCAACGGGTTTTTAACATTAAGTGGTCAAAATTGTGCACTTTGCTTTCTGGGTTTGTTATTCT";
		remo->map[RAT_SPECIES] = 
			"ATACACTCACAGTGGTTTGGCATATATTTGGTGAATTTATTTTTAAGGAAAAAAATTAGTGTTGGTTTCGATATATGGTAGCTTTCTCTCTAACATAATTTGAATAATTCAGCAAAGCCCCACTACCAGCTGCTCTTCTGCAGCCCTTTTCCATTCTTTTCAGCATTATAATTTTGGTTAATTTTCAATTTTAGGTCCTACGTCTCTGCAATTTGTGTATGAATAACAGAATAATTTCCCTCTTTTGTTTCGCCTTTCCTGTTCCTGAATCTAAATAAAGATGGCTTTTTAGTATTAAAAGTGGAAGAAAATTACAGGTAATTATCTTTGACGGTAAAAACGCTGTAATCAGCGGGCTACATGAAAAATTACTCTAATTATGGCTGCATTTAAGAGAATGGAAAAAAACCTTCTTGTGGATAAAAACCTTAAATTGTCCCCAATGTCTGCTTCAAATTGGATGGCACTGCAGCTGGAGGCTTTGTTCAGAATTGATCCTGGGGAGCTACAAACCCAAAGTTTCACAGTAGGAAGGGGGAAAAAAGAAAAGAAAACATTTTTCCTAATGTAACAATGCGAATGGTAGAAAATGACAAGACTGATCGGTTTTAAACCATTCTGAAGACTGACTGAGTGTGGAAGTTGCTCAACAAAAAAAGGAACGGGGATATTGGTAAGTAGTCTTTGCTTTGTGTCTAATTATAGTGACTGTCCTTTTGCACGACTGTATGTCGCTGTCATTATATGGAGCACTCTGAGTATCTCTATTGACTTCTGATAAATGGCTCAGTGGTTTCAAGGTTCATAATTTGAAGGATGCCCAGGTCTGTAGCCATTAATTCTCACAGTATGGGGCTTCCTGCTTGTATGACGAAAGGACTTTTCATTTTCATTATTTCTTTGCTTTCCACCAGATCGGAAGGATGTGTTTCTTGCAGGCTAGAGCATTTTCCTACCTTGCTAGCTTGCAACGGGTTCTCAACATTAAGCAGTCAAAATTGTGCACTTCGCTTTCTGGGTTTGCATTCTT";
		remo->map[HUMAN_SPECIES] = 
			"ATGCATTTGCCGGGTTTCAGAATGTTATGCACTCACAGTGGTTTGGCATGCATCTGGTGAATTTTTTTTAACGAAAAATTAGTGTTGGTTTCGATGTATGGTAGCATTCTCCCTAACGTAATTTGAATAATTCAGCAAAGCCCCACTACCAGCTGTACTTCTGCAGCCTCTTCCATTCTTTTCAGCATTATAATTTTGGTTAATTTTCAATTTTAGGTCCTACGTCTCTGCAATTTGTGTATGAATAACAGAATAATTTCCCTCTTTTGTTTCGCCTTTCCTGTTCCTGAATCTAAATAAAGATGGCTTTTTAGTATTAAAAGTGGAAGAAAATTACAGGTAATTATCTTTGACGGTAAAAACGCTGTAATCAGCGGGCTACATGAAAAATTACTCTAATTATGGCTGCATTTAAGAGAATGGAAAAAAACCTTCTTGTGGATAAAAACCTTAAATTGTCCCCAATGTCTGCTTCAAATTGGATGGCACTGCAGCTGGAGGCTTTGTTCAGAATTGATCCTGGGGAGCTACGAACCCAAAGTTTCACAGTAGGAAGGGGGAAAAAAGAAAAGAAAACATTTTTCCTAATGTAACAATGCGAATGCTAGAAAATGACAAGACTGATCGGTTTTAAACCATTCTGAAGACTGACTGAGCGTGGAAGTTGCTCAACAAAAAAAGGAACGGGTATATTGGTAAGTAGTCTTTGCTTTGTGTCTAATTATAGTGACTGTCCTTTTGCACGACTGTATGTCGCTGTCATTATATGGAGCACTCTGAGTATCTCTATTGACTTCTGATAAATGGCTCAGTGGTTTCAAGGTTCATAATTTGAAGGATGCCCAGGTCTGTAGCCATTAATTCTCGCAGTATGGGGCTTCCTGCTTGTATGACGAAAGGACTTTTCATTTTCATTATTTCTGTGCTTTCCCCAAGATCGGAAGGATGTATTTCATCCAGCATAGTGCATTTTCCTACCTAGCTAGCTGGCAATGGGTTTTAAACATTAA";
		remo->map[DOG_SPECIES] = 
			"TGGGTTTCAGAATTCCATGCACTCACAGTGGTTTGGCATACGTCTGGTAGATTTTTTTTTTTTTTAAGAAAAATGAGTGTTGGTTTCGATGTATGGTAGCTTTTTCTCTAACGTAATTTGAATAATTCAGCAAAGCCCTACTACCAGCTGTACTTCTGCAGCCTCTTCCATTCTTTTCAGCATTATAATTTTGGTTAATTTTCAATTTTAGGTCCTACGTCTCTGCAATTTGTGTATGAATAACAGAATAATTTCCCTCTTTTGTTTCGCCTTTCCTGTTCCTGAATCTAAATAAAGATGGCTTTTTAGTATTAAAAGTGGAAGAAAATTACAGGTAATTATCTTTGACGGTAAAAACGCTGTAATCAGCGGGCTACATGAAAAATTACTCTAATTATGGCTGCATTTAAGAGAATGGAAAAAAACCTTCTTGTGGATAAAAACCTTAAATTGTCCCCAATGTCTGCTTCAAATTGGATGGCACTGCAGCTGGAGGCTTTGTTCAGAATTGATCCTGGGGAGCTACGAACCCAAAGTTTCACAGTAGGAAGGGGGAAAAAAGAAAAGAAAACATTTTTCCTAATGTAACAATGCGAATGCTAGAAAATGACAAGACTGATCGGTTTTAAACCATTCTGAGACTGAGCGTGGAAGTTGCTCAACAAAAAAAGGAACGGGTATATTGGTAAGTAGTCTTTGCTTTGTGTCTAATTATAGTGACTGTCCTTTTGCACGACTGTATGTCGCTGTCATTATATGGAGCACTCTGAGTATCTCTATTGACTTCTGATAAATGGCTCAGTGGTTTCAAGGTTCATAATTTGAAGGATGCCCAGGTCTGTAGCCATTAATTCTCGCAGTATGGGGCTTCCTGCTTGTATGACGAAAGGGCTTTTCATTTTCATTATTTCTGTGCTTTCCACAAGATCGGAAGGATGTATTTCATCCAGGATCATGCATTTTCCTACCTAGCTAGCTTGCAAAGGGTGTAAAACATTAAGTAGTTATAATTGTGCATTTACTTTGTTTT";
		remo->map[XENOPUS_SPECIES] = 
			"CCCCGCACAAAGCTCCAGGTTTAGCATAATGTGAATAATTCAGCGATGCTCTGAGCAGCTGTAATTCTGAAATCTTTCTCCCTTCTTTTTACTATTATAATTCTGGTTAATTTTCAATTTTAGGTCCTACATCTCTGCAATTTGTGTATGAATAACAGAATAATTTCCCTTTTTTGTTTCGCTTTTCCTGTTCCTGAATCTAAATAAAGATGGCTTTTTAGTATTAAAAGTGGAAGAAAATTACAGGTAATTATCTTTGACGGTAAAAACGCTGTAATCAGCGGGCTACATGAAAAATTACTCTAATTATGGCTGCATTTAAGAGAATGGGGGAAAAAAAACCTTCTTCTGGATAAAAACGTTAAATTGTCCCCAATGTGTGCTTGTAACTGGATGGGACGGCAGCTGCAGGCTTTGTTCAGCAATGATCCTGGGGAGATGTGAACCCAAGATGCACAGTGGGAAGGGGGAGAGATAAGGCGGCACTGTTCAATAACATAACAATACGAGACATAGAAATGACAAGACTGACCGGTTTAAGCCATTCTGAAGACTGGAACTTGCTCAGCGTCAAGAACTAAAAAAAATAAATATAAAAAACAAAAACACACACTCACACACACAAAACACTGGTGTATTGGTAAGTAATCTTTGCTTTGTGCCTAATTATAGTGACTGTCCTTTTGCACGACTGTGTGTAGCTGTCATTATATCGAGCACTCTCTGACTGTCTATTGACTCCTCATAAGTGGCTCAGTGCTTTCAAGGTTCAGGATTTGAAAGATGCCCAGGTCGCTCCTATTAATTCTCCCAGCAAGGGGCATCTAACGAACATGGAAATGCAGCTTCTTTATTTTAGCTTTGACTAGATTTGC";
		remo->map[TETRAODON_SPECIES] = 
			"ACTTGACAGCTGCGCGGAGCTTTTCAGCGTCGTAATTTCAGATAATATCCCGAGCGCTCACTCCTCTCCGGCAATTTGTATATGAATAACCGAATAATTTCCCCCTTTTGTTCCATCTGTGCTACCTCAAATCCAAATAAAGATGCCTTTTAGTATTAAAAGTGGTAGAAAATTACAGGTAATTATCTTTGACGGTAAAAACGCTGTAATCAGCGGGCTACATCAAAAATTACCCTAATTATGCCTGCATTTATGAGAATGGGGGGAAAAAAAGCTTCACTTGGATAAAAACCCACAAATTGTCCCAAATATCTCCTTCATATTGAACGCCGCTGCAGCTGGCTGCTTTGTTCCACATCGATCCCGCGGAGTTTGGCATTTAGCGGCCCGCATTAGGACGAAGAGACGGGAGAAAACATTTTGTCAATAC";
		test_remos["dlx5_2"] = remo;

		already_built = true;
	}
};

BIO_NS_END