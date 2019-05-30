#!/usr/bin/env bash

for arg in "$@"; do
    case $arg in
    org.stjude.compbio.sam.*)
        CLASSPATH=/opt/tools-sam/lib/*:$CLASSPATH
        break
        ;;
    org.stjude.compbio.xenocp.*)
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
