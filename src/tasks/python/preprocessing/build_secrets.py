import json
import os

google_keys = ['type', 'project_id', 'private_key_id', 'private_key', 
                'client_email', 'client_id', 'auth_uri', 'token_uri',
                'auth_provider_x509_cert_url', 'client_x509_cert_url']
aws_keys = ['AWS_ACCESS_KEY_ID', 'AWS_SECRET_ACCESS_KEY']

def make_json(keys):
	return {k: os.environ[k] for k in keys}

def main():
	with open('/root/.cloudvolume/secrets/google-secret.json', 'w') as f:
		json.dump(make_json(google_keys), f)
	with open('/root/.cloudvolume/secrets/aws-secret.json', 'w') as f:
		json.dump(make_json(aws_keys), f)

if __name__ == '__main__':
	main()