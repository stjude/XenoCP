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
    && wget https://github.com/biod/sambamba/releases/download/v0.7.0/sambamba-0.7.0-linux-static.gz \
    && echo "5a739ea53ef296639825831e110e359eab6ff421e108c7e3f4df0d67859e3024 *sambamba-0.7.0-linux-static.gz" | sha256sum --check \
    && gunzip sambamba-0.7.0-linux-static.gz \
    && mv sambamba-0.7.0-linux-static /usr/local/bin/sambamba \
    && chmod 0755 /usr/local/bin/sambamba

RUN cd /tmp \
    && wget https://services.gradle.org/distributions/gradle-5.4-bin.zip \
    && echo "c8c17574245ecee9ed7fe4f6b593b696d1692d1adbfef425bef9b333e3a0e8de *gradle-5.4-bin.zip" | sha256sum --check \
    && unzip gradle-5.4-bin.zip \
    && mv gradle-5.4 /opt/gradle

RUN mkdir -p /opt/picard/lib \
    && cd /opt/picard/lib \
    && wget -O picard-2.6.0.jar https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar \
    && echo "671d9e86e6bf0c28ee007aea55d07e2456ae3a57016491b50aab0fd2fd0e493b *picard-2.6.0.jar" | sha256sum --check

ENV PATH /opt/gradle/bin:${PATH}

COPY bin /tmp/xenocp/bin
COPY src /tmp/xenocp/src
COPY dependencies /tmp/xenocp/dependencies
COPY build.gradle /tmp/xenocp/build.gradle
COPY settings.gradle /tmp/xenocp/settings.gradle

RUN cd /tmp/xenocp \
    && gradle installDist \
    && cp -r build/install/xenocp /opt

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
COPY --from=builder /usr/local/bin/sambamba /usr/local/bin/sambamba
COPY --from=builder /opt/picard /opt/picard
COPY --from=builder /opt/xenocp /opt/xenocp
COPY --from=builder /opt/xenocp/bin/* /usr/local/bin/

COPY cwl /opt/xenocp/cwl

ENTRYPOINT ["cwl-runner", "--outdir", "results", "/opt/xenocp/cwl/xenocp.cwl"]
