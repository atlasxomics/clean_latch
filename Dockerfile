FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:dd8f-main

# Install R
RUN apt-get update -y && \
    apt-get install -y software-properties-common && \
    add-apt-repository "deb http://cloud.r-project.org/bin/linux/debian buster-cran40/" && \
    apt-get install -y \
        r-base \
        r-base-dev \
        r-cran-devtools 

# Install tabix for bgzip command
RUN apt-get install -y tabix

RUN R -e "install.packages(c('dplyr'))"

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN echo "hi"
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root

