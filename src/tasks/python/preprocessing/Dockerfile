# Use an official Python runtime as a parent image
FROM python:2.7-alpine

# Attach author
MAINTAINER macrintr

# Download git
RUN apk update && apk upgrade && \
    apk add --no-cache bash git py-pip

# Copy git repo
RUN git clone https://github.com/seung-lab/Alembic.git

# Set the working directory to ../preprocessing
WORKDIR /Alembic/src/tasks/python/preprocessing/

# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt

# Setup cloud-volume dirs & vars
RUN mkdir -p ~/.cloudvolume/secrets/

# Create symlink to /run/secrets/.. for cloud-volume
RUN ln -s /run/secrets ~/.cloudvolume/secrets/

# Run run.py when the container launches
CMD ["python", "run.py"]