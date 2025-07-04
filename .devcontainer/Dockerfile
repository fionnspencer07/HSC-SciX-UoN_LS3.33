FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install R and dependencies
RUN apt-get update && apt-get install -y \
    dirmngr gnupg apt-transport-https ca-certificates software-properties-common \
    libssl-dev libcurl4-openssl-dev libxml2-dev libjpeg-dev libpng-dev \
    sudo wget locales gdebi-core \
    && locale-gen en_US.UTF-8

# Add R repository and install R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 'E298A3A825C0D65DFD57CBB651716619E084DAB9' \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/' \
    && apt-get update \
    && apt-get install -y r-base

# Install RStudio Server
RUN wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2023.12.1-402-amd64.deb \
    && gdebi -n rstudio-server-2023.12.1-402-amd64.deb \
    && rm rstudio-server-2023.12.1-402-amd64.deb

# Set password for the default user (codespace)
RUN echo 'codespace:rstudio' | chpasswd

# Expose RStudio Server port
EXPOSE 8787
