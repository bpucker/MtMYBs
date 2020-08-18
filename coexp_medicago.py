### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python coexp_ath.py
					--in <FULL_PATH_TO_CANDIDATE_FILE(one ID per line)>
					--exp <FULL_PATH_TO_EXPRESSION_TABLE>
					--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
					
					optional:
					--mapping <GENE_NAME_MAPPING_TABLE>
					--ann <FULL_PATH_TO_ANNOTATION_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

from operator import itemgetter
import numpy as np
import re, math, sys, os
from scipy import stats

# ---- end of imports --- #

def load_expression_values( filename ):
	"""! @brief load all expression values """
	
	expression_data = {}
	with open( filename, "r" ) as f:
		tissues = f.readline().strip().split('\t')[1:]
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			expression = {}
			for idx, each in enumerate( parts[1:] ):
				expression.update( { tissues[  idx ] : float( parts[ idx+1 ] ) } )
			line = f.readline()
			expression_data.update( { parts[0]: expression } )
	return expression_data


def compare_candidates_against_all( candidate, gene_expression ):
	"""! @brief compare candidate gene expression against all genes to find co-expressed genes """
	
	tissues = sorted( gene_expression[ gene_expression.keys()[0] ].keys() )
	coexpressed_genes = []
	for i, gene2 in enumerate( gene_expression.keys() ):
		if candidate != gene2:
			values = []
			total_expression = 0
			for tissue in tissues:
				try:
					x = gene_expression[ candidate ][ tissue ]
					y = gene_expression[ gene2 ][ tissue ]
					total_expression += y
					if not math.isnan( x ) and not math.isnan( y ) :
						values.append( [ x, y ] )
				except KeyError:
					pass	#print tissue
			r, p = stats.spearmanr( values )
			if not math.isnan( r ) and total_expression > 30:
				if r > 0.3 and p < 0.05:
					coexpressed_genes.append( { 'id': gene2, 'correlation': r, 'p_value': p } )
		
	return coexpressed_genes


def load_annotation( annotation_file ):
	"""! @brief load annotation mapping table """
	
	annotation_mapping_table = {}
	
	with open( annotation_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			annotation_mapping_table.update( { parts[0]: parts[1] } )
			line = f.readline()
	return annotation_mapping_table


def load_gene_name_mapping_table( gene_name_mapping_table_file ):
	"""! @brief load name mapping table """
	
	mapping_table = {}
	with open( gene_name_mapping_table_file, "r" ) as f:
		line =f.readline()
		while line:
			parts = line.strip().split('\t')
			mapping_table.update( { parts[0]: parts[1] } )
			line = f.readline()
	return mapping_table


def main( arguments ):
	"""! @brief run everything """
	
	bv_expression_file = arguments[ arguments.index( '--exp' )+1 ]
	output_prefix = arguments[ arguments.index( '--out' )+1 ]
	candidate_gene_file = arguments[ arguments.index( '--in' )+1 ]
	
	
	if "--mapping" in arguments:
		gene_name_mapping_table_file = arguments[ arguments.index( '--mapping' )+1 ]
		gene_name_mapping_table = load_gene_name_mapping_table( gene_name_mapping_table_file )
	else:
		gene_name_mapping_table = {}
	
	if output_prefix[-1] != "/":
		output_prefix += "/"
	
	if not os.path.exists( output_prefix ):
		os.makedirs( output_prefix )
	
	if '--ann' in arguments:
		annotation_file = arguments[ arguments.index( '--ann' )+1 ]
		annotation_mapping_table = load_annotation(  annotation_file)
	else:
		annotation_mapping_table = {}
	
	gene_expression = load_expression_values( bv_expression_file )
	
	high_impact_candidates = [ ]
	
	with open( candidate_gene_file, "r" ) as f:
		line = f.readline()
		while line:
			high_impact_candidates.append( line.strip() )
			line = f.readline()
	
	number_of_genes = float( len( gene_expression.keys() ) )
	for candidate in high_impact_candidates:
		coexpressed_genes = sorted( compare_candidates_against_all( candidate, gene_expression ), key=itemgetter( 'correlation' ) )[::-1]
		coexpression_output_file = output_prefix + candidate + ".txt"
		with open( coexpression_output_file, "w", 0 ) as out:
			out.write( 'CandidateGene\tGeneID\tSpearmanCorrelation\tadjusted_p-value\tFunctionalAnnotation\n' )
			for entry in coexpressed_genes:
				try:
					ann = annotation_mapping_table[ gene_name_mapping_table[ entry['id'] ] ]
				except KeyError:
					try:
						ann = annotation_mapping_table[ entry['id'] ]
					except KeyError:
						"N/A"
				
				try:
					gene1 = gene_name_mapping_table[ candidate ]
				except KeyError:
					gene1 = candidate
				
				try:
					gene2 = gene_name_mapping_table[ entry['id'] ]
				except KeyError:
					gene2 = entry['id']
				
				out.write( "\t".join( map( str, [ 	gene1,
																	gene2,
																	entry['correlation'],
																	entry['p_value'] * number_of_genes,
																	ann
																] ) ) + '\n' )


if '--exp' in sys.argv and '--ann' in sys.argv and '--out' in sys.argv and '--in' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
