# Use an official Python runtime as a parent image
FROM python:2.7

# Attach author
MAINTAINER macrintr

# Download git & other packages
RUN apt-get update
RUN apt-get install -y libblas-dev liblapack-dev liblapacke-dev gfortran git

# Copy git repo
RUN git clone https://github.com/seung-lab/Alembic.git

# Set the working directory to ../preprocessing
WORKDIR /Alembic/src/tasks/python/preprocessing

# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt

# Setup cloud-volume dirs & vars
RUN mkdir -p ~/.cloudvolume/secrets/

# Create symlink to /run/secrets/.. for cloud-volume
RUN ln -s /run/secrets ~/.cloudvolume/secrets/

# Run run.py when the container launches
ENTRYPOINT ["python", "run.py"]