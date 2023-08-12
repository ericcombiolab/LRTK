import argparse
import datetime
import os
import subprocess
import sys


def logging(message):
	print(f"[{datetime.datetime.now()}] {message}", flush=True)


def run_cmd(command, name, verbose=True):
	if verbose:
		logging(f'{name} starts')
	ret = subprocess.call(command, shell=True)
	if ret != 0:
		logging(f'{name} fails: {"".join(command)}')
		sys.exit(1)
	elif verbose:
		logging(f'{name} ends')


