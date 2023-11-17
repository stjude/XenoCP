#!/usr/bin/env bash

# If the classpath is already set, then delegate directly to java
if [ "$CLASSPATH" != "" ]; then exec java "$@"; fi

# Otherwise, build an appropriate classpath
# This section assumes you are running inside the container
for arg in "$@"; do
    case $arg in
    org.stjude.compbio.*)
        CLASSPATH=/opt/xenocp/lib/*:$CLASSPATH
        break
        ;;
    picard.cmdline.*)
        CLASSPATH=/opt/picard/lib/picard-2.6.0.jar:$CLASSPATH
        break
        ;;
    esac
done

java -classpath "$CLASSPATH" "$@"
