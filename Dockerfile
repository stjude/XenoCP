FROM ubuntu:18.04 as builder

RUN apt-get update \
    && apt-get --yes install \
        build-essential \
        openjdk-11-jdk-headless \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN cd /tmp \
    && wget https://github.com/lh3/bwa/releases/download/v0.7.13/bwa-0.7.13.tar.bz2 \
    && echo "559b3c63266e5d5351f7665268263dbb9592f3c1c4569e7a4a75a15f17f0aedc *bwa-0.7.13.tar.bz2" | sha256sum --check \
    && tar xf bwa-0.7.13.tar.bz2 \
    && cd bwa-0.7.13 \
    && make -j$(nproc) \
    && mv bwa /usr/local/bin

# bz2 and lzma support is for CRAM files. curses is for `samtools tview`.
RUN cd /tmp \
    && wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && echo "083f688d7070082411c72c27372104ed472ed7a620591d06f928e653ebc23482 *samtools-1.9.tar.bz2" | sha256sum --check \
    && tar xf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && ./configure \
        --prefix /usr/local \
        --disable-bz2 \
        --disable-lzma \
        --without-curses \
    && make -j$(nproc) \
    && make install

RUN cd /tmp \
    && wget https://services.gradle.org/distributions/gradle-5.3.1-bin.zip \
    && echo "1c59a17a054e9c82f0dd881871c9646e943ec4c71dd52ebc6137d17f82337436 *gradle-5.3.1-bin.zip" | sha256sum --check \
    && unzip gradle-5.3.1-bin.zip \
    && mv gradle-5.3.1 /opt/gradle

RUN mkdir -p /opt/picard/lib \
    && cd /opt/picard/lib \
    && wget -O picard-2.6.0.jar https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar \
    && echo "671d9e86e6bf0c28ee007aea55d07e2456ae3a57016491b50aab0fd2fd0e493b *picard-2.6.0.jar" | sha256sum --check

ENV PATH /opt/gradle/bin:${PATH}

COPY common-java-genome /tmp/xenocp/common-java-genome
COPY common-java-sam /tmp/xenocp/common-java-sam
COPY util-java /tmp/xenocp/util-java
COPY tools-sam /tmp/xenocp/tools-sam
COPY xenocp /tmp/xenocp/xenocp
COPY build.gradle /tmp/xenocp/build.gradle
COPY settings.gradle /tmp/xenocp/settings.gradle

RUN cd /tmp/xenocp \
    && gradle installDist \
    && cp -r tools-sam/build/install/tools-sam /opt \
    && cp -r xenocp/build/install/xenocp /opt

FROM ubuntu:18.04

RUN apt-get update \
    && apt-get --yes install \
        gawk \
        nodejs \
        openjdk-11-jre-headless \
        python3 \
        python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install cwlref-runner html5lib

RUN ln -s /usr/bin/gawk /bin/awk

COPY --from=builder /usr/local/bin/bwa /usr/local/bin/bwa
COPY --from=builder /usr/local/bin/samtools /usr/local/bin/samtools
COPY --from=builder /opt/picard /opt/picard
COPY --from=builder /opt/tools-sam /opt/tools-sam
COPY --from=builder /opt/xenocp /opt/xenocp

COPY bin/java-settmp.sh /usr/local/bin/java-settmp.sh
COPY bin/java.sh /usr/local/bin/java.sh
COPY bin/qclib.sh /usr/local/bin/qclib.sh
COPY bin/view_awk_picard.sh /usr/local/bin/view_awk_picard.sh

COPY cwl /opt/xenocp/cwl

COPY mapping-standard/src/main/scripts/merge_markdup_index.sh /usr/local/bin/merge_markdup_index.sh
COPY mapping-standard/src/main/scripts/qc_bam.sh /usr/local/bin/qc_bam.sh

COPY tools-sam/src/main/scripts/sam_to_single.awk /usr/local/bin/sam_to_single.awk
COPY tools-sam/src/main/scripts/sort_flagstat.sh /usr/local/bin/sort_flagstat.sh

COPY xenocp/src/main/scripts/bwa_alignse_onlymapped.sh /usr/local/bin/bwa_alignse_onlymapped.sh

ENTRYPOINT ["cwl-runner", "--outdir", "results", "/opt/xenocp/cwl/xenocp.cwl"]
