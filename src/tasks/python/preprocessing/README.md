# Preprocess Chunk Task  
Use [cloud-volume](https://github.com/seung-lab/cloud-volume) to run the 
following preprocessing tasks on Precomputed chunks:
* Contrast stretching
* Edge masking
* Section removal & transposing

## Scheduling  
Task scheduling based on AWS SQS.  
Edit `sqs_publisher.py` with task iterator.
Run `python sqs_publisher.py <SQS_QUEUE_NAME>`.
Most recent Docker image created @ `macrintr/alembic-preprocessing`.  
Create cluster in gcloud. Suggested parameters:   
* 1 vCPU
* >5GB RAM
* ~32GB boot disk
Attach secrets to the cluster:
`kubectl create secret generic secrets --from-file=<HOME>/.cloudvolume/secrets/aws-secret.json --from-file=<HOME>/.cloudvolume/secrets/google-secret.json`
Run `job.yaml`  
Scale deployment