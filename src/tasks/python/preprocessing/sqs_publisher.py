import boto.sqs
from boto.sqs.message import Message
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

def create_preprocessing_task(origin, chunk_dims):
    args = {}
    args['inlayer'] = 'gs://neuroglancer/pinky100_v0/image'
    args['outlayer'] = 'gs://neuroglancer/pinky100_v0/image_corrected'
    args['contrastlayer'] = 'gs://neuroglancer/pinky100_v0/contrast'
    args['masklayer'] = 'gs://neuroglancer/pinky100_v0/edge_mask'
    args['orderlayer'] = 'gs://neuroglancer/pinky100_v0/z_order_corrected'
    args['x_slice'] = [origin[0], origin[0]+chunk_dims[0]]
    args['y_slice'] = [origin[1], origin[1]+chunk_dims[1]]
    args['z_slice'] = [origin[2], origin[2]+chunk_dims[2]]
    return json.dumps(args)

def publish_messages(queue_name):
	"""Publishes multiple messages to SQS queue."""
	q = get_queue(queue_name)
	chunk_dims = [1024,1024,64]
	x_start = 35840
	x_stop = 35840+100352 - chunk_dims[0]
	y_start = 26624
	y_stop = 26624+62464 - chunk_dims[1]
	z_start = 1
	z_stop = 1+2240 - chunk_dims[2]

	n = 0
	for x in range(x_start, x_stop, chunk_dims[0]):
	    for y in range(y_start, y_stop, chunk_dims[1]):
	        for z in range(z_start, z_stop, chunk_dims[2]):
		        data = create_preprocessing_task([x,y,z], chunk_dims)
		        # Data must be a bytestring
		        data = data.encode('utf-8')
		        m = Message()
		        m.set_body(data)
		        q.write(m)
		        n += 1

	print('Published {0} messages.'.format(str(n)))

def main(queue_name):
	publish_messages(queue_name)

if __name__ == "__main__":
	main(sys.argv[1])