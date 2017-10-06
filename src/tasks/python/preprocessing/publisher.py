from google.cloud import pubsub_v1
import json

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

def publish_messages(project, topic_name):
	"""Publishes multiple messages to a Pub/Sub topic."""
	publisher = pubsub_v1.PublisherClient()
	topic_path = publisher.topic_path(project, topic_name)
	chunk_dims = [1024,1024,64]
	# x_start = 35840
	# x_stop = 35840+100352 - chunk_dims[0]
	# y_start = 26624
	# y_stop = 26624+62464 - chunk_dims[1]
	# z_start = 1
	# z_stop = 64 #1+2240 - chunk_dims[2]

	# test
	x_start = 35840 
	x_stop = 35840 + 48*64 #78272 # 24 chunks
	y_start = 75200 - 48*64
	y_stop = 75200
	z_start = 1
	z_stop = 64

	n = 1
	for x in range(x_start, x_stop, chunk_dims[0]):
	    for y in range(y_start, y_stop, chunk_dims[1]):
	        for z in range(z_start, z_stop, chunk_dims[2]):
		        data = create_preprocessing_task([x,y,z], chunk_dims)
		        # Data must be a bytestring
		        data = data.encode('utf-8')
		        publisher.publish(topic_path, data=data)

	print('Published messages.')

def main():
	project = 'neuromancer-seung-import'
	topic_name = 'alembic_preprocessing'
	publish_messages(project, topic_name)

if __name__ == "__main__":
	main()