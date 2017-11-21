import os
import sys
from pdb2zernike import zernike_moments
from search_pdb import set_default_db_path
from multiprocessing import Pool
from libtbx import easy_pickle
from scitbx.math import zernike_align_fft as fft_align
from scitbx import math
from sastbx.interface import get_input
import iotbx.phil
import time

global nmax 
global nlm_array_ref
global coefs 
global nlm_total
global codes
global pdbfile



base_path = os.path.split(sys.path[0])[0]

master_params = iotbx.phil.parse("""\
	retrieval{
	pdbfile = None
	.type=path
	.help="the protein model"

	dbpath = None
	.type = path
	.help = "the path of the protein database"

	nmax = 10
	.type=int
	.help="maximum order of zernike polynomial: FIXED for the existing database"

	db_prefix='mydb'
	.type=path
	.help="the prefix of database filename"

	prefix = "retrivel"
	.type=path
	.help="the output prefix"


	}
	"""
	)

banner = "-------------------Searching the protein DATABASE for similar shapes-------------------"

def help( out=None ):
  if out is None:
    out= sys.stdout
  print >> out, "\nUsage: \n"
  print >> out, "   sastbx.retrieval pdbfile=pdbfile [rmax=rmax nmax=nmax scan=True*/False buildmap=True*/False pdb=pdbfile path=database_path]\n"
  print >> out, "   The pdb file is the only required input file  (in theory)\n"
  print >> out, "   Optional control parameters:"
  print >> out, "     nmax     : maximum order of the zernike polynomial expansion (<=20 for precomputed database; 10 is the default)"
  print >> out, "     dbpath     : path to the database (this MUST be correct to execute the searching)"
  print >> out, "     db_user_prefix: database filename"
  print >> out, "     prefix   : the output prefix\n\n"


def set_default_db_path():
  global stdfile
  import libtbx.env_config
  env = libtbx.env_config.unpickle()
  sastbx_path = env.dist_path("sastbx")
  path = sastbx_path+'/smalldatabase/'
  print  "\nATTENTION: Database path was set to : >>%s<<"%path

  return path

def calcc(modelId):
	nlm_array_mov=math.nlm_array(nmax)
	nlm=nlm_array_ref.nlm()
  	nlm_array_mov.load_coefs(nlm, coefs[modelId][0:nlm_total])
	align_obj=fft_align.align(nlm_array_ref, nlm_array_mov, nmax=nmax, refine=True)
	cc=align_obj.get_cc()
	print "c.c. between ", os.path.split(pdbfile)[-1], "and ", codes[modelId],"is ", cc
	return cc




def run(args):
	targetfile = os.path.join(os.path.split(sys.path[0])[0],"retrieval.txt")
	with open(targetfile,"w") as f:
		f.truncate()

	time1 = time.time()
	global nmax 
	global nlm_array_ref
	global coefs 
	global nlm_total
	global codes
	global pdbfile
	params = get_input(args, master_params, "retrieval", banner, help)
	if( params is None):
		exit()
	pdbfile = params.retrieval.pdbfile
	dbpath = params.retrieval.dbpath
	nmax = params.retrieval.nmax
	dbprefix = params.retrieval.db_prefix
	prefix = params.retrieval.prefix

	print "=============process the protein model=============="
	with open(targetfile,"a") as f:
		print >> f, "=============process the protein model=============="


	zernike_moments(pdbfile,nmax=nmax)
	
	queryCoefFile = pdbfile.split(".")[0]+".nlm.pickle"
	#queryCoefFile=pdbfile.replace("pdb", "nlm.pickle")
	queryCoef=easy_pickle.load(queryCoefFile)
	with open(targetfile,"a") as f:
		print >> f,"=============load database==============="
	print "=============load database==============="
	if (dbpath is None):
		dbpath = set_default_db_path()
		codes = easy_pickle.load(os.path.join(dbpath,dbprefix+".codes"))
		coefs = easy_pickle.load(os.path.join(dbpath,dbprefix+".nlm"))

	else:
		codes = easy_pickle.load(os.path.join(dbpath,dbprefix+".codes"))
		coefs = easy_pickle.load(os.path.join(dbpath,dbprefix+".nlm"))
	
	with open(targetfile,"a") as f:
		print >>f, "=============database============="
		print >>f, os.path.join(dbpath,dbprefix+".codes")	
		print >>f, os.path.join(dbpath,dbprefix+".nlm")
		print >>f, "=================================="


	print "=============database============="
	print os.path.join(dbpath,dbprefix+".codes")
	print os.path.join(dbpath,dbprefix+".nlm")
	print "=================================="


	nmodels = len(coefs)
	nlm_array_ref=math.nlm_array(nmax)
	nlm=nlm_array_ref.nlm()
	nlm_total=nlm_array_ref.coefs().size()
	nlm_array_ref.load_coefs(nlm, queryCoef[0:nlm_total])




	p = Pool(8)
	cclist = p.map(calcc,range(nmodels))
	

	distlist=[1-cc for cc in cclist]
	rankedlist=sorted(range(nmodels), key=lambda k:distlist[k])	
	rankedcodes=[codes[rank] for rank in rankedlist]
	sortedcclist = sorted(cclist, reverse = True)
	
	with open(targetfile,"a") as f:
		print >> f, "=========Tope 10 models matching the input protein model============"
		
	print "=========Tope 10 models matching the input protein model============"
	
	with open(targetfile,"a") as f:
		for i in range(10):
			print "top ",(i+1), " ", rankedcodes[i], "c.c.", sortedcclist[i]
			print >>f , "top ",(i+1), " ", rankedcodes[i], "c.c.", sortedcclist[i]


	time2 = time.time()
	print "time used:", time2-time1
	with open(targetfile,"a") as f:
		print >>f ,"time used: ", time2-time1
	
	# pdbfile = args[0]
	# #dbpath = args[1]
	# #nmax = args[2]
	# dbpath = "/Users/songna/Desktop/smalldatabase/"
	# nmax = 20
	# pdb_rank(pdbfile,dbpath,nmax)

if __name__ =="__main__":
	args = sys.argv[1:]
	time1 = time.time()
	run(args)
	time2 = time.time()
	print "time used:", time2-time1




