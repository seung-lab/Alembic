import boto.sqs
from boto.sqs.message import Message
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

def main(queue_name, fn):
	q = get_queue(queue_name)
	rs = q.get_messages(10)
	while len(rs) > 0:
		with open(fn, 'a') as f:
			for r in rs:
				m = r.get_body().encode('latin-1')
				f.write(m + '\n')
				q.delete_message(r)
		rs = q.get_messages(10)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])