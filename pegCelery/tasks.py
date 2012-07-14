from __future__ import absolute_import

from pegCelery.celery import celery
from subprocess import call


@celery.task
def runCalculation(argumentString):
	return call("./pegSerial", argumentString.split())
