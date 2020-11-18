#!/bin/sh

hostname=`hostname --short`
date=`date --iso-8601=ns`

# Debug sleep to check sync / async
# sleep 7

echo "$hostname: $date"