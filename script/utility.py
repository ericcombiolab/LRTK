import argparse
import datetime
import os
import subprocess
import sys


def logging(message):
	print(f"[{datetime.datetime.now()}] {message}", flush=True)


def run_cmd(command, name, verbose=True):
	if verbose:
		logging(f'{name} starts: {" ".join(command)}')
	ret = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
	if ret.returncode:
		logging(f'{name} fails: {" ".join(command)}')
		sys.exit(1)
	elif verbose:
		logging(f'{name} ends: {" ".join(command)}')


