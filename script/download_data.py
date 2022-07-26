import argparse
import os
import sys
import logging
import wget

try:
	from urllib.request import urlretrieve
except ImportError:
	from urllib import urlretrieve

def download_data(db_path):
	"""
	Download a file from a url
	"""

	DirLocal = db_path + "/DATABASE"

	if not os.path.isdir(DirLocal):
		try:
			os.makedirs(DirLocal)
		except:
			logging.info("ERROR: Unable to create folder for database install: " + DirLocal)

	DBdict={
#		'UniqnessMap':{'Server':'http://xinzhouneuroscience.org/wp-content/uploads/2019/05/Uniqness_map.tar.gz', 'Local':db_path + "/DATABASE/Uniqness_map.tar.gz"},
		'RefGenome':{'Server':'https://www.dropbox.com/s/cepdgqth94s3hut/GSE132024-GPL16791_series_matrix.txt.gz?dl=0', 'Local':db_path + "/DATABASE/refdata-GRCh38-2.1.0.tar.gz"}		
	}

	for key in DBdict.keys():
		file_download=DBdict[key]['Local']
		file_url=DBdict[key]['Server']
		if not os.path.isfile(file_download):
			try:
				print("\nDownloading " + file_url + "\n")
				logging.info("\nDownloading " + file_url + "\n")
				#file, headers = urlretrieve(file_url, file_download)
				wget.download(file_url, file_download)
			except:
				logging.info("\nERROR: Unable to download " + file_url + "\n")
		else:
			logging.info("\nFile {} already present!\n".format(file_download))


download_data("/tmp/local/cschaoyang/PIPELINE/LRTK/TEST/")
