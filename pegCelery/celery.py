from __future__ import absolute_import

from celery import Celery

celery = Celery(broker='amqp://',
                backend='amqp://',
                include=['pegCelery.tasks'])

# Configuration required to work with PHP-Celery
celery.conf.update(
    CELERY_RESULT_SERIALIZER = "json",
    CELERY_TASK_RESULT_EXPIRES = None
)
