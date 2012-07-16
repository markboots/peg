from __future__ import absolute_import

import subprocess
from celery import Celery

celery = Celery()

@celery.task
def runCalculation(argumentString):
	return subprocess.call(["./pegSerial"] + argumentString.split())
