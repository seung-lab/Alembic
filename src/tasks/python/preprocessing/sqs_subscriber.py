import boto.sqs
from boto.sqs.message import Message
import preprocessing
import time
import json
import os
import sys

def get_connection():
	secrets_subpath = ".cloudvolume/secrets/aws-secret.json"
	secret_path = os.path.join(os.path.expanduser('~'), secrets_subpath)
	with open(secret_path) as f:    
	    aws_secrets = json.load(f)
	return boto.sqs.connect_to_region(
			"us-east-1",
			aws_access_key_id=aws_secrets['AWS_ACCESS_KEY_ID'],
			aws_secret_access_key=aws_secrets['AWS_SECRET_ACCESS_KEY'])

def get_queue(queue_name):
	conn = get_connection()
	return conn.get_queue(queue_name)

def receive_messages(queue_name):
	"""Receives messages from a pull subscription."""
	q = get_queue(queue_name)
	while True:
		messages = q.get_messages(1)
		if len(messages) > 0:
			m = messages[0]
			message = m.get_body().encode('latin-1')
			preprocessing.main(message)
			q.delete_message(m)
		else:
			time.sleep(5)

if __name__ == "__main__":
    receive_messages(sys.argv[1])