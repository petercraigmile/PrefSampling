#!/bin/bash

tar xvf Original/ghcnm.tavg.latest.qcf.tar.gz


cat ghcnm.v4.0.1.20210902/ghcnm.tavg.v4.0.1.20210902.qcf.inv | grep "^USC" |  sed "s/*//g" > Derived/US_stations.txt


cat ghcnm.v4.0.1.20210902/ghcnm.tavg.v4.0.1.20210902.qcf.dat | grep "^USC" |  sed "s/*//g" > Derived/US_monthly_temps.txt
