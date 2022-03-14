##############################################################
# Experimental Docker file for the Artemis software suite.
#
# Please see usage documentation at Docker Hub:
#
# https://hub.docker.com/r/sangerpathogens/artemis
#
##############################################################
FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

LABEL maintainer="path-help@sanger.ac.uk"

# Build and environment configuration
ARG ARTEMIS_BUILD_DIR=/artemis-build
ARG INSTALL_DIR=/opt/artemis
ARG BLAST_INSTALL_DIR=/usr/local/ncbi/blast/
ARG BROWSER_INSTALL_DIR=/usr/local/firefox
ARG ARTEMIS_WORKDIR=/artemis
ENV JAVA_OPTS="-Dsun.jnu.encoding=UTF-8 -Dfile.encoding=UTF-8"
ENV LANG="en_US.UTF-8"
ENV LANGUAGE="en_US.UTF-8"
ENV LC_ALL="en_US.UTF-8"
ENV LC_COLLATE="en_US.UTF-8"
ENV LC_CTYPE="en_US.UTF-8"

# Mount points
RUN mkdir -p $ARTEMIS_WORKDIR && chmod 777 $ARTEMIS_WORKDIR

# Install main dependencies
RUN apt-get update --quiet --assume-yes
RUN apt-get upgrade --quiet --assume-yes
RUN apt-get install --quiet --assume-yes locales openjdk-11-jdk libxtst6 libxrender1 libxext6 libexpat1 fonts-dejavu-core fontconfig-config libfontconfig1 libfreetype6 libpng16-16 curl wget maven ant
RUN update-alternatives --config java
RUN update-alternatives --config javac

RUN echo "Java version installed: `java -version`"

# Set locale
RUN echo en_US.UTF-8 UTF-8 >> /etc/locale.gen
RUN locale-gen en_US.UTF-8

# Build and install Artemis
# No tests are run currently as some of them pop up windows and hence require X
RUN mkdir -p $INSTALL_DIR
RUN mkdir -p $ARTEMIS_BUILD_DIR
COPY . $ARTEMIS_BUILD_DIR
RUN cd $ARTEMIS_BUILD_DIR && mvn validate && mvn -Djarsigner.skip=true -Dskip.tests=true package -P release
RUN tar xvf $ARTEMIS_BUILD_DIR/target/release-artifacts/unix-release/artemis-unix-release-*.tar.gz --strip-components=1 -C ${INSTALL_DIR}
RUN rm -rf $ARTEMIS_BUILD_DIR && rm -rf ~/.m2

# Install Blast+
RUN mkdir -p $BLAST_INSTALL_DIR
RUN curl --fail -L ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz | tar xzf - --strip-components=1 -C $BLAST_INSTALL_DIR

# Install Firefox for use with pfam etc (use firejail for increased security)
RUN apt-get install --quiet --assume-yes libgtk-3-0 libdbus-glib-1-2 libxt6 firejail
RUN mkdir -p $BROWSER_INSTALL_DIR
RUN wget -O firefox.tar.bz2 "https://download.mozilla.org/?product=firefox-latest-ssl&os=linux64&lang=en-GB" && \
    tar xvf /firefox.tar.bz2 && \
    rm -f /firefox.tar.bz2 && \
    mv /firefox /usr/local/

RUN echo "firejail --private --dns=8.8.8.8 --dns=8.8.4.4 $BROWSER_INSTALL_DIR/firefox -no-remote "'"$@"' > /usr/bin/firefox && chmod 555 /usr/bin/firefox

# Cleanup
RUN apt-get autoremove
RUN apt-get clean

# Environment
ENV PATH=$INSTALL_DIR:$INSTALL_DIR/dist:$INSTALL_DIR/etc:$BLAST_INSTALL_DIR/bin:$PATH
ENV CLASSPATH=$INSTALL_DIR:$INSTALL_DIR/dist:$INSTALL_DIR/etc:${CLASSPATH}
ENV HOME=$ARTEMIS_WORKDIR
WORKDIR $ARTEMIS_WORKDIR

# Define default command.
CMD   echo 'Usage:  docker run -d -e DISPLAY="<your display name>:0" -v <your user home directory>:/home/artuser -v <your data folder>:/artemis' && \
      echo '          --user $(id -u):$(id -g) -e ARTEMIS_JVM_FLAGS="-Duser.home=/home/artuser -Djava.io.tmpdir=/tmp" --tmpfs /tmp' && \
      echo '          --rm artemis <program name [art|act|bamview|dnaplotter]> [program arguments]' && \
      echo && \
      echo 'For help, please go to http://sanger-pathogens.github.io/Artemis/'

