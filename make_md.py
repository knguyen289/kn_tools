#Copy paste into terminal:
#python make_md.py -fnames basic_tools/basic_tools.py,ternary_tools/ternary_tools.py,rna_tools/rna_tools.py -o README.md

import argparse
import os

def getOptions():
	"""Function to pull in arguments"""
	parser = argparse.ArgumentParser()

	parser.add_argument('-fnames',dest='fnames',type=str,help='List of filenames')
	parser.add_argument('-o',dest='output',type=str,help='Name of output file')

	args = parser.parse_args()
	return args

def main(args):
	fnames = args.fnames
	fnames = fnames.rstrip()
	fnames = fnames.split(',')
	doc = '# kn_tools\nKim Nguyen\'s tools for Wang Lab, includes tools for plotting and data management\n'
	for fname in fnames:
		try:
			f = open(fname,'r')
		except:
			print fname, 'could not be opened.'
			break

		doc += '## ' + fname + '\n'

		in_desc = False
		func = ''
		desc = ''
		param = False
		ret = False
		params = []
		rets = []

		for line in f:
			line = line.lstrip()
			line = line.rstrip()

			if len(line) > 0:
				if line.find('\"\"\"') != -1:
					in_desc = not in_desc
				if line.find('def ') == 0:
					if func != '':
						doc += '### ' + func + '()\n'
						doc += '#### Description:\n' + '* ' + desc + '\n\n'

					if len(params) > 0:
						doc += '#### Parameters:\n' 
						for p in params:
							paren1 = p.index('(')
							paren2 = p.index(')')
							p_name = p[:paren1-1]
							p_type = p[paren1:paren2+1]
							p_rest = p[paren2+1:]
							doc += '* **' + p_name + '** *' + p_type  + '*' + p_rest + '\n\n'
					
					if len(rets) > 0:
						doc += '#### Returns:\n'
						for r in rets:
							paren1 = r.index('(')
							paren2 = r.index(')')
							r_name = r[:paren1-1]
							r_type = r[paren1:paren2+1]
							r_rest = r[paren2+1:]
							doc += '* **' + r_name + '** *' + r_type  + '*' + r_rest + '\n'

					doc += '\n'
					func = ''
					desc = ''
					param = False
					ret = False
					params = []
					rets = []
					func = line.split(' ')[1].split('(')[0]
				elif in_desc:
					if func != '' and desc == '':
						desc = line[3:]
					elif line.find('Parameters:') != -1:
						param = True
						ret = False
						pass
					elif line.find('Returns:') != -1:
						ret = True
						param = False
						pass
					elif func != '' and desc != '' and param:
						params.append(line)
					elif func != '' and desc != '' and not param and ret:
						rets.append(line)
		if func != '':
			doc += '### ' + func + '()\n'
			doc += '#### Description:\n' + '* ' + desc + '\n\n'
			
		if len(params) > 0:
			doc += '#### Parameters:\n' 
			for p in params:
				paren1 = p.index('(')
				paren2 = p.index(')')
				p_name = p[:paren1-1]
				p_type = p[paren1:paren2+1]
				p_rest = p[paren2+1:]
				doc += '* **' + p_name + '** *' + p_type  + '*' + p_rest + '\n\n'
		
		if len(rets) > 0:
			doc += '#### Returns:\n'
			for r in rets:
				paren1 = r.index('(')
				paren2 = r.index(')')
				r_name = r[:paren1-1]
				r_type = r[paren1:paren2+1]
				r_rest = r[paren2+1:]
				doc += '* **' + r_name + '** *' + r_type  + '*' + r_rest + '\n'

	try:
		f = open(args.output,'w')
	except:
		os.makedirs(args.output)
		f = open(args.output,'w')

	f.truncate()
	f.write(doc)

if __name__ == '__main__':
	args = getOptions()
	main(args)
		


