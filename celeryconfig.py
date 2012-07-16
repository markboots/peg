## Broker settings.
BROKER_URL = "amqp://guest:guest@localhost:5672//"

# List of modules to import when celery starts.
CELERY_IMPORTS = ("pegCelery", )

## Using amqp also to send task state and results.
CELERY_RESULT_BACKEND = "amqp"

# Settings required for PHP AMQP extension
CELERY_RESULT_SERIALIZER = "json"
CELERY_TASK_RESULT_EXPIRES = None

