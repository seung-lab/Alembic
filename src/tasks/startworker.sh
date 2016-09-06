#! /bin/bash
cd /home/ubuntu/.julia/v0.4/SimpleTasks/src/
git pull
sed -i -e 's/pmap/map/g' ./services/datasource.jl
cd /home/ubuntu/.julia/v0.4/ImageRegistration/src/
git pull
cd /home/ubuntu/Alembic/src/
git pull
sudo -i -u ubuntu -H sh -c "stdbuf -oL -eL julia -p 35 /home/ubuntu/Alembic/src/tasks/runawsdaemon.jl | tee -a /home/ubuntu/daemon.out &"
