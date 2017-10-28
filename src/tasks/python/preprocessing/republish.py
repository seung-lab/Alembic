import boto.sqs
from boto.sqs.message import Message
import sys
import os
import json

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

def main(queue_name, fn):
	q = get_queue(queue_name)
	with open(fn, 'r') as f:
		for ln in f.readlines():
			# Data must be a bytestring
			data = ln.encode('utf-8')
			m = Message()
			m.set_body(data)
			q.write(m)

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])