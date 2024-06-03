#!/usr/bin/env bash

datadir='/home/bodrug-a/bin/beacon2-ri-api/deploy/data/'

# AIC EXOME DATA
if [[ "$1" == *"aic"* ]] ; then
  echo "Loading... ", $1
  echo ${datadir}${1}
  for term in analyses biosamples cohorts datasets individuals runs genomicVariations; 
    do
    sudo docker cp ${datadir}${1}/${term}.json rimongo:tmp/${term}.json;
    sudo docker exec rimongo mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/${term}.json --collection $term ; 
    done
  echo "Loaded "
# 1000 GENOMES DATA, EXOMES
elif [[ "$1" == *"onekgenomes"* ]]; then
  echo "Loading...", $1
  echo ${datadir}${1}
  for term in analyses biosamples cohorts datasets individuals runs genomicVariations; 
    do
    sudo docker cp ${datadir}${1}/${term}.json rimongo:tmp/${term}.json;
    sudo docker exec rimongo mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/${term}.json --collection $term ; 
    done
  echo "Loaded "
# 1kGP hg38
elif [[ "$1" == *"1kGP-hg38"* ]]; then
  echo "Loading...", $1
  echo ${datadir}${1}
  for term in analyses biosamples cohorts datasets individuals runs genomicVariations; 
    do
    #cp /home/bodrug-a/Devlopment/aic_metadata/bff_1kGP_GRCh38/${term}.json ${datadir}${1}/${term}.json
    sudo docker cp ${datadir}${1}/${term}.json rimongo:tmp/${term}.json;
    sudo docker exec rimongo mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/${term}.json --collection $term ; 
    done
  echo "Loaded "
# CINECA DATA 
# Cineca data aw directly in the data folder, i moved it
elif [[ "$1" == "cineca" ]] ; then
  echo "Loading CINECA data..."
  echo ${datadir}${1}
  for term in analyses biosamples cohorts datasets individuals runs genomicVariations; 
    do
    sudo docker cp ${datadir}${1}/${term}.json rimongo:tmp/${term}.json;
    sudo docker exec rimongo mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/${term}.json --collection $term ; 
    done
  echo "Loaded"

fi

sudo docker exec beacon python beacon/reindex.py
sudo docker exec beacon python beacon/db/extract_filtering_terms.py
