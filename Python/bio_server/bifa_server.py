import biopsy.transfac, biopsy.pssm_paths
import os.path, re, webbrowser
import biopsy.analyse_remos.seqs as seqs
import biopsy.analyse_remos.paths as paths
from tempfile import mkstemp

from soaplib.serializers.primitive import String, Integer, Array, Boolean, Float
from soaplib.serializers.clazz import ClassSerializer

import math

OTT_ALGORITHM=0
BAYESIAN_ALGORITHM=1

def pssmAccs(pssm_sets,use_consensus_sequences,matrix_species,matrix_name_match):
	if pssm_sets != None :
		pssm_accs = biopsy.SequenceVec()
		#	Make a straight copy of the longer transfac sequence first, and then append
		#	any other sequence sets as required
		for pssm_set in pssm_sets :
			if pssm_set == "transfac" :
				pssm_filter = biopsy.transfac.PssmFilter(use_consensus_sequences,matrix_species,matrix_name_match)
				pssm_accs = biopsy.get_transfac_pssm_accessions( pssm_filter )
		for pssm_set in pssm_sets :
			if pssm_set == "all-custom" :
				biopsy.append_sequences(pssm_accs,biopsy.get_all_custom_pssms())
			elif pssm_set != "transfac" :
				biopsy.append_sequences(pssm_accs,biopsy.get_custom_pssms(pssm_set))
	else :
		pssm_filter = biopsy.transfac.PssmFilter(use_consensus_sequences,matrix_species,matrix_name_match)
		pssm_accs = biopsy.get_transfac_pssm_accessions( pssm_filter )
	return pssm_accs


def bifa_tool(sequence,threshold,title,algorithm,show_labels,
			  phylo_sequences,use_consensus_sequences,matrix_species,phylo_threshold,
			  matrix_name_match,use_old_algorithm,use_cumulative_likelihoods,pssm_sets):

	biopsy.PssmParameters.use_p_value = (algorithm == OTT_ALGORITHM)
	
	if (algorithm == OTT_ALGORITHM) :
		biopsy.PssmParameters.use_cumulative_dists = 1
	else :
		biopsy.PssmParameters.use_cumulative_dists = use_cumulative_likelihoods

	if None == title :
		title = ""

	pssm_accs = pssmAccs(pssm_sets,use_consensus_sequences,matrix_species,matrix_name_match)

#	hit_map = { }
#	include_max_chain=0
	open_svg=False
#	tag="all"


	fd, fname = mkstemp()
	os.close(fd)
	temp=fname.split(os.sep)
	remo=temp[2]

#	Build svg args common to both old and new algorithms
	notes = 'Species: %s\nMatrix filter: %s\nThreshold: %f\nPhylo: %f' % (matrix_species,matrix_name_match,threshold,phylo_threshold)
	build_svg_args = biopsy.BuildSvgArgs(
									file = '%s.svg' % ( remo ),
									title = title,
									#  max_threshold = 0..  C++ and Python default.
									min_threshold = threshold,
									max_num_factors = 16,   	# Python default = 10, Original dll = 16
									show_labels=show_labels,
									notes = notes,
									open_file = open_svg
							)

	if not use_old_algorithm:

		sequences = biopsy.SequenceVec()
		sequences.append(sequence)
		if phylo_sequences!=None:
			for each in phylo_sequences:
				sequences.append(each)

		hits, max_chain, unadjusted_hits = biopsy.score_pssms_on_phylo_sequences \
			(pssm_accs, sequences,threshold = threshold, phylo_threshold = phylo_threshold)
#	if include_max_chain:
#		hit_map[ remo ] = (hits, max_chain)
#		else:
#		hit_map[ remo ] = hits

#		print len(hits)

		if len( hits ):
			biopsy.build_analysis_svg(hits,
									  max_chain,
									  sequence,
									  args = build_svg_args
									  )
#		print 'Done remo %s with phylo_threshold = %f' % ( remo, phylo_threshold )
	else :
		#	Old algorithm
		if use_cumulative_likelihoods :
			raise Exception("Cumulative likelihoods not supported with old algorithm")

		sequences = biopsy.SequenceVec()
		if phylo_sequences!=None:
			for each in phylo_sequences:
				sequences.append(each)

		pssm_match_args = biopsy.PssmMatchArgs(threshold,
							 use_consensus_sequences,
							 matrix_species,
							 matrix_name_match,
							 (algorithm == BAYESIAN_ALGORITHM))

		biopsy.pssm_match(sequence,
						  sequences,
						  pssm_match_args,
						  build_svg_args,
						  phylo_threshold)


	return remo


class BiFaHit(ClassSerializer):
	class types:
		Identifier = String
		Name = String
		Position = Integer
		Length = Integer
		Positive_strand = Boolean
		Value = Float


def bifa_hits(sequence,threshold,algorithm,
			  phylo_sequences,use_consensus_sequences,matrix_species,phylo_threshold,
			  matrix_name_match,use_cumulative_likelihoods,pssm_sets):

	biopsy.PssmParameters.use_p_value = (algorithm == OTT_ALGORITHM)

	if (algorithm == OTT_ALGORITHM) :
		biopsy.PssmParameters.use_cumulative_dists = 1
	else :
		biopsy.PssmParameters.use_cumulative_dists = use_cumulative_likelihoods

	pssm_accs = pssmAccs(pssm_sets,use_consensus_sequences,matrix_species,matrix_name_match)

	sequences = biopsy.SequenceVec()
	sequences.append(sequence)
	if phylo_sequences!=None:
		for each in phylo_sequences:
			sequences.append(each)

	hits, max_chain, unadjusted_hits = biopsy.score_pssms_on_phylo_sequences \
			(pssm_accs, sequences,threshold = threshold, phylo_threshold = phylo_threshold)

	output = []
	for hit in hits :
		h = BiFaHit()
		h.Identifier = hit.binder
		h.Position = hit.location.position
		h.Length = hit.location.length
		h.Positive_strand = hit.location.positive_strand
		h.Value = hit.p_binding

		if (h.Identifier.find('-') > 0) :
			pssm = biopsy.get_custom_pssm(hit.binder)
			h.name = pssm.name
		elif h.Identifier[0] == 'M' :
			pssm = biopsy.transfac.Matrix(hit.binder)
			h.Name = pssm.name
		elif h.Identifier[0] == 'R' :
			site = biopsy.transfac.Site(hit.binder)
			h.Name = site.name
		else :
			h.Name = 'Unknown'
		output.append(h)

	return output

def score_pssm_on_sequence(pssm_name,sequence,threshold):
#	biopsy.score_pssm_on_sequence.__doc__
#	biopsy.PssmParameters.use_p_value = 1
#	biopsy.PssmParameters.use_cumulative_dists = 1

	hits = biopsy.HitVec()

	# print "method score_pssm_on_sequence called"
	biopsy.score_pssm_on_sequence(pssm_name=pssm_name, sequence=sequence, threshold=threshold, result=hits)

	output=[]
	counter=0
	for hit in hits :
		if hit.location.position == int(counter/2) :
			output.append(hit.p_binding)
			counter += 1
		else :
			raise Exception ("Error marshaling hit data")
	return output

def score_pssms_on_sequences(pssm_names,sequences,algorithm):
#	biopsy.score_pssm_on_sequence.__doc__

#	If using p-value then it must also be cumulative.  Non cumulative does not work
	biopsy.PssmParameters.use_p_value = (algorithm == OTT_ALGORITHM)
	biopsy.PssmParameters.use_cumulative_dists = (algorithm == OTT_ALGORITHM)

	pssmNo = len(pssm_names)
	seqNo = len(sequences)

	output = []		# print "method score_pssm_on_sequence called"
	if pssmNo > 1 :
		if seqNo > 1 :
			raise Exception("Can't have multiple sequences and pssms")

		for i in range(pssmNo):
			hitArray = ""
			hits = biopsy.HitVec()
			try :
				biopsy.score_pssm_on_sequence(pssm_name=pssm_names[i], sequence=sequences[0], threshold = 0, result=hits)
				j=0
				for hit in hits :
					if hit.location.position == int(j/2) :
						hitArray += '%g,' % hit.p_binding
						j = j+1
					else :
						raise Exception("")
				output.append(hitArray)
			except Exception :
				output.append("")
	else :
		for seq in sequences :
			hitArray = ""
			hits = biopsy.HitVec()
			try :
				biopsy.score_pssm_on_sequence(pssm_name=pssm_names[0], sequence=seq, threshold = 0, result=hits)
				j=0
				for hit in hits :
					if hit.location.position == int(j/2) :
						hitArray += '%g,' % hit.p_binding
						j += 1
					else :
						raise Exception("")
				output.append(hitArray)
			except Exception :
				output.append("")
	return output



def get_pssm_set_names():
	return biopsy.get_custom_pssm_set_names()


class PssmInfoData(ClassSerializer):
	class types:
		name = String
		pathway = String
		length = Integer
		pseudo_count = Float
		sites = Integer
		factors = Array(String)


def get_pssm_info(pssm_name):
	v = PssmInfoData()
	pssm = biopsy.get_pssm(pssm_name)
	v.pseudo_count = pssm.pseudo_count
	v.sites = pssm.sites

	if (pssm_name.find('-') > 0) :
		pssm = biopsy.get_custom_pssm(pssm_name)
		v.name = pssm.name
		v.length = len(pssm.counts)
		v.pathway = "Unknown"
		v.factors = []
		return v
	if pssm_name[0]== 'R' :
		pssm = biopsy.transfac.Site(pssm_name)
	elif pssm_name[0]== 'M' :
		pssm = biopsy.transfac.Matrix(pssm_name)
	else :
		v.name = 'Unknown'
		return v
	#Common code for both sites and Matrixes
	p = pssm.pathway
	if p.known :
		v.pathway = pssm.pathway.name
	else :
		v.pathway = "Unknown"
	v.name = pssm.name
	v.length = pssm.size
	factors = []
	for i in range(len(pssm.factors)):
		s = str(pssm.factors[i].link) + "; " + pssm.factors[i].name + ";"
		if (len(pssm.factors[i].species) > 0) :
		    s = s + " Species: "
		for j in range(len(pssm.factors[i].species)):
			if (j > 1):
				s = s + ", "
			s = s + pssm.factors[i].species[j]
		factors.append(s)
	v.factors = factors

	return v

def get_pssm_freqs(pssm_name):
	pssm = biopsy.get_pssm(pssm_name)
	dists = pssm.dists
	freqs = []
	for n in dists :
		f = ""
		f += "%f," % n.get_freq(0)
		f += "%f," % n.get_freq(1)
		f += "%f," % n.get_freq(2)
		f += "%f," % n.get_freq(3)
		freqs.append(f)
	return freqs

def get_pssm_counts(pssm_name):
	pssm = biopsy.get_pssm(pssm_name)
	p = pssm.counts
	counts = []
	for n in p :
		f = ""
		f += "%f," % n.get(0)
		f += "%f," % n.get(1)
		f += "%f," % n.get(2)
		f += "%f," % n.get(3)
		counts.append(f)
	return counts

