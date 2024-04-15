FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

# Install R
RUN apt-get update -y && \
    apt-get install -y software-properties-common

# Install tabix for bgzip command
RUN apt-get install -y tabix

# Install python packages
RUN pip install --upgrade pip
COPY requirements.txt /root/requirements.txt
RUN python3 -m pip install -r requirements.txt

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
