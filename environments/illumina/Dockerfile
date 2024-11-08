FROM condaforge/miniforge3:latest AS condabuild
LABEL authors="Matt Bull" \
      description="Docker image containing all requirements for an Illumina ncov2019 pipeline"

COPY environments/extras.yml /extras.yml
COPY environments/illumina/environment.yml /environment.yml

RUN conda env create -f /environment.yml

FROM debian:buster-slim
RUN apt-get update && \
apt-get install -y git procps && \
apt-get clean -y 
COPY --from=condabuild /opt/conda/envs/artic-ncov2019-illumina /opt/conda/envs/artic-ncov2019-illumina
ENV PATH=/opt/conda/envs/artic-ncov2019-illumina/bin:$PATH
