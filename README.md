[![Build Status](https://travis-ci.org/seung-lab/Julimaps.svg?branch=master)](https://travis-ci.org/seung-lab/Julimaps)

# Alembic
ALignment of Electron Microscopy By Image Correlograms
A set of tools for elastic image registration in Julia.

## Alignment
Register coarsely-aligned sections to each other, and k-nearest neighboring sections, by blockmatching. Filter matches automatically and manually. Cover each section image in a triangle mesh, and globally deform all the meshes elastically to accommodate the correspondences. Render with a piecewise affine transform.

## Requirements  
* julia v0.6
* [cloud-volume](https://github.com/seung-lab/cloud-volume)

## Getting started  
Setup the Julia wrapper for [CloudVolume](https://github.com/seung-lab/CloudVolume.jl).

Here's a simple startup script:
```
using Alembic
load_params(joinpath(homedir), ".julia/v0.6/Alembic/src/params/pinky40_test.json"); # specifies CloudVolume dirs and params
ms = make_stack(); # make meshset given the params specified in the JSON file
match!(ms); # blockmatch between the meshes in the meshset, per the params
elastic_solve!(ms); # relax the spring system
render(ms); # render the images and save to CloudVolume
```

## Task Scheduling  
You can schedule and run blockmatching and render tasks in parallel across a cluster. We recommend you use one of the following Docker containers & Kubernetes.
* [Alembic with Julia-MKL precompiled on a Google Compute Engine instance](https://hub.docker.com/u/macrintr/alembic:julia-mkl-gce)
* [Generic Alembic](https://hub.docker.com/u/sergiypopo/alembic:generic)

If you'd like to create a container for precompiled Julia+MKL on a different instance type, do the following:
* On the instance type, install Docker.
* Clone this repository: `git clone https://github.com/seung-lab/Alembic.git`
* `cd Alembic`
* `docker build -t [TAG NAME] .`


### Publishing/Subsribing to Tasks
We've included blockmatch and render tasks in `src/tasks/tasks.jl`, along with task publisher/subscriber scripts based on AWS SQS. We assume that you've setup CloudVolume to include the file `~/.cloudvolume/secrets/aws-secret.json`. You can publish tasks in the following manner:

```
julia ~/.julia/v0.6/Alembic/src/tasks/sqs_publisher.jl [PARAMS JSON PATH] [QUEUE NAME]
```

You can view your scheduled tasks in the [AWS SQS console](https://console.aws.amazon.com/sqs/).

You can setup a queue subscriber in a similar fashion:

```
julia ~/.julia/v0.6/Alembic/src/tasks/sqs_subscriber.jl [QUEUE NAME] [NO. OF PROCS TO USE]
```

Pay attention to the `Default Visibility Timeout` attribute of your SQS queue. You'll want to set it to be longer than the time your task is expected to take. That way, you won't have workers working on the same task simultaneously.

### Setting up Kubernetes
Once you have tasks scheduled in the queue, you can set up a cluster that can be managed with Kubernetes. Here's how to do it using Google Container Engine.

First, make sure you have [kubectl](https://kubernetes.io/docs/tasks/tools/install-kubectl/) installed on your local workstation.

Now, create a cluster in Google Cloud (use the console or CLI). We recommend using the following parameters:
```
Zone: us-east1-b
Cluster Version: 1.7.8-gke.0 (default)
Machine type: 8 vCPU, 30 GB memory
Boot disk size in GB (per node): 32
```

Creating cluster through gcloud CLI:
```
gcloud container clusters create sergiy-test-alembic --zone=us-east1-b --machine-type=n1-standard-8 --disk-size=32
```

Connecting to cluster through glcoud CLI:
```
gcloud container clusters get-credentials your-cluster-name 
```




It's important that you make sure the zone is set to match your storage zone. Size your machine type & boot disk according to the size of your data. Boot disk is not used beyond container installation, but it can become a constraining resource for your entire project, so set it accordingly.

Once your cluster has been created, within the console, select the `Connect` button to access the `kubectl` commands to connect your local workstation to your cluster. Now you can mount your secrets to your cluster nodes:

```
kubectl create secret generic secrets --from-file=$HOME/.cloudvolume/secrets/aws-secret.json --from-file=$HOME/.cloudvolume/secrets/google-secret.json
```

Now you can intialize the container+commands to be run on each node (we recommend a kube deployment, rather than a job). See `deployment.yaml` for an example deployment definition. We can run that deployment with the following command:

```
kubectl create -f $HOME/.julia/v0.6/Alembic/src/tasks/deployment.yaml
```

And now you're running on a cluster! If you want to adjust the number of pods that are supporting the deployment, you can use this following command:

```
kubectl scale deployment alembic-tasks --replicas=[NO. OF PODS]
```

You can investigate the different Kubernete resources with the following commands:

```
kubectl describe deployment
kubectl describe pods
kubectl get deployment
kubectl get pods
kubectl logs [POD NAME]  # 
```

For more commands, see the [kubectl Cheat Sheet](https://kubernetes.io/docs/user-guide/kubectl-cheatsheet/).

If you ever need to replace the deployment (e.g. you updated the docker container or you want to change the command you pass to the container), just delete the deployment and create a new one.

```
kubectl delete deployment alembic-tasks
kubectl create -f [HOMEDIR FULL PATH]/.julia/v0.6/Alembic/src/tasks/deployment.yaml
```
