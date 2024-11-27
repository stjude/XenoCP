FROM ubuntu:20.04 as builder

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get --yes install \
        build-essential \
        openjdk-11-jdk-headless \
        unzip \
        wget \
        python3 \
        python3-pip \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --ignore-installed \
        --prefix /usr/local \
        cwlref-runner \
        html5lib \
        urllib3==1.26.15

RUN cd /tmp \
    && wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
    && echo "de1b4d4e745c0b7fc3e107b5155a51ac063011d33a5d82696331ecf4bed8d0fd *bwa-0.7.17.tar.bz2" | sha256sum --check \
    && tar xf bwa-0.7.17.tar.bz2 \
    && cd bwa-0.7.17 \
    && make -j$(nproc) \
    && mv bwa /usr/local/bin

RUN cd /tmp \
    && wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10a.tar.gz \
    && echo "af0df8fdc0e7a539b3ec6665dce9ac55c33598dfbc74d24df9dae7a309b0426a *2.7.10a.tar.gz" | sha256sum --check \
    && tar xf 2.7.10a.tar.gz \
    && mv STAR-2.7.10a/bin/Linux_x86_64_static/STAR /usr/local/bin

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
COPY gradle /tmp/xenocp/gradle
COPY gradlew /tmp/xenocp/gradlew
COPY build.gradle /tmp/xenocp/build.gradle
COPY settings.gradle /tmp/xenocp/settings.gradle

RUN cd /tmp/xenocp \
    && ./gradlew installDist \
    && cp -r build/install/xenocp /opt

FROM ubuntu:20.04

RUN apt-get update \
    && apt-get --yes install --no-install-recommends \
        gawk \
        nodejs \
        openjdk-11-jre-headless \
        python3 \
        python3-distutils \
        python-is-python3 \
        file \
    && rm -rf /var/lib/apt/lists/*

COPY --chmod=755 --from=builder /usr/local/bin/cwl* /usr/local/bin/
COPY --chmod=755 --from=builder /usr/local/lib /usr/local/lib/
COPY --chmod=755 --from=builder /usr/local/bin/bwa /usr/local/bin/bwa
COPY --chmod=755 --from=builder /usr/local/bin/STAR /usr/local/bin/STAR
COPY --chmod=755 --from=builder /usr/local/bin/samtools /usr/local/bin/samtools
COPY --chmod=755 --from=builder /usr/local/bin/sambamba /usr/local/bin/sambamba
COPY --chmod=755 --from=builder /opt/picard /opt/picard
COPY --chmod=755 --from=builder /opt/xenocp /opt/xenocp
COPY --chmod=755 --from=builder /opt/xenocp/bin/* /usr/local/bin/

COPY --chmod=755  cwl /opt/xenocp/cwl