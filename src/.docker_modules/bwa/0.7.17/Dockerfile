FROM lbmc/sambamba:0.6.9
MAINTAINER Laurent Modolo

ENV BWA_VERSION=0.7.17
ENV SAMBLASTER_VERSION=0.1.24

ENV PACKAGES curl=7.58* \
    ca-certificates=20180409 \
    build-essential=12.4* \
    zlib1g-dev=1:1.2.11*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN curl -k -L https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2 -o bwa-v${BWA_VERSION}.tar.bz2 && \
tar xjf bwa-v${BWA_VERSION}.tar.bz2 && \
cd bwa-${BWA_VERSION}/ && \
make && \
cp bwa /usr/bin && \
cd .. && \
rm -R bwa-${BWA_VERSION}/

RUN curl -k -L https://github.com/GregoryFaust/samblaster/releases/download/v.${SAMBLASTER_VERSION}/samblaster-v.${SAMBLASTER_VERSION}.tar.gz -o samblaster-v.${SAMBLASTER_VERSION}.tar.gz && \
tar xvf samblaster-v.${SAMBLASTER_VERSION}.tar.gz && \
cd samblaster-v.${SAMBLASTER_VERSION}/ && \
make && \
cp samblaster /usr/bin && \
cd .. && \
rm -R samblaster-v.${SAMBLASTER_VERSION}/
