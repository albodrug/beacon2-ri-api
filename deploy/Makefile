build:
	docker-compose up -d --build beacon permissions db
	sleep 10
	docker exec rimongo /bin/bash -c 'mongo beacon -u root -p example --authenticationDatabase admin --eval "db.filtering_terms.deleteMany({})"'
	docker exec rimongo /bin/bash -c 'mongo beacon -u root -p example --authenticationDatabase admin --eval "db.analyses.deleteMany({})"'
	docker exec rimongo /bin/bash -c 'mongo beacon -u root -p example --authenticationDatabase admin --eval "db.biosamples.deleteMany({})"'
	docker exec rimongo /bin/bash -c 'mongo beacon -u root -p example --authenticationDatabase admin --eval "db.cohorts.deleteMany({})"'
	docker exec rimongo /bin/bash -c 'mongo beacon -u root -p example --authenticationDatabase admin --eval "db.datasets.deleteMany({})"'
	docker exec rimongo /bin/bash -c 'mongo beacon -u root -p example --authenticationDatabase admin --eval "db.genomicVariations.deleteMany({})"'
	docker exec rimongo /bin/bash -c 'mongo beacon -u root -p example --authenticationDatabase admin --eval "db.individuals.deleteMany({})"'
	docker exec rimongo /bin/bash -c 'mongo beacon -u root -p example --authenticationDatabase admin --eval "db.runs.deleteMany({})"'
	docker cp data/analyses.json rimongo:tmp/analyses.json
	docker cp data/biosamples.json rimongo:tmp/biosamples.json
	docker cp data/cohorts.json rimongo:tmp/cohorts.json
	docker cp data/datasets.json rimongo:tmp/datasets.json
	docker cp data/genomicVariations.json rimongo:tmp/genomicVariations.json
	docker cp data/individuals.json rimongo:tmp/individuals.json
	docker cp data/runs.json rimongo:tmp/runs.json
	docker exec rimongo mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/datasets.json --collection datasets
	docker exec rimongo mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/individuals.json --collection individuals
	docker exec rimongo mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/cohorts.json --collection cohorts
	docker exec rimongo mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/analyses.json --collection analyses
	docker exec rimongo mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/biosamples.json --collection biosamples
	docker exec rimongo mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/runs.json --collection runs
	docker exec rimongo mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" --file /tmp/genomicVariations.json --collection genomicVariations
	docker cp /data/vault/bio-scratch/arnau/beacon/beacon2-ri-tools-v2_test_anot/tcga-coad/analyses.json mongoega:tmp/analyses.json
	docker cp /data/vault/bio-scratch/arnau/beacon/beacon2-ri-tools-v2_test_anot/tcga-coad/biosamples.json mongoega:tmp/biosamples.json
	docker cp /data/vault/bio-scratch/arnau/beacon/beacon2-ri-tools-v2_test_anot/tcga-coad/cohorts.json mongoega:tmp/cohorts.json
	docker cp /data/vault/bio-scratch/arnau/beacon/beacon2-ri-tools-v2_test_anot/tcga-coad/datasets.json mongoega:tmp/datasets.json
	docker cp /data/vault/bio-scratch/arnau/beacon/beacon2-ri-tools-v2_test_anot/tcga-coad/genomicVariations.json mongoega:tmp/genomicVariations.json
	docker cp /data/vault/bio-scratch/arnau/beacon/beacon2-ri-tools-v2_test_anot/tcga-coad/individuals.json mongoega:tmp/individuals.json
	docker cp /data/vault/bio-scratch/arnau/beacon/beacon2-ri-tools-v2_test_anot/tcga-coad/runs.json mongoega:tmp/runs.json
	docker cp /data/vault/bio-scratch/arnau/beacon/beacon2-ri-tools-v2_test_anot/tcga-coad/filtering_terms.json mongoega:tmp/filtering_terms.json
	docker exec mongoega mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27020/beacon?authSource=admin" --file /tmp/datasets.json --collection datasets
	docker exec mongoega mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27020/beacon?authSource=admin" --file /tmp/individuals.json --collection individuals
	docker exec mongoega mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27020/beacon?authSource=admin" --file /tmp/cohorts.json --collection cohorts
	docker exec mongoega mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27020/beacon?authSource=admin" --file /tmp/analyses.json --collection analyses
	docker exec mongoega mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27020/beacon?authSource=admin" --file /tmp/biosamples.json --collection biosamples
	docker exec mongoega mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27020/beacon?authSource=admin" --file /tmp/runs.json --collection runs
	docker exec mongoega mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27020/beacon?authSource=admin" --file /tmp/genomicVariations.json --collection genomicVariations
	docker exec mongoega mongoimport --jsonArray --uri "mongodb://root:example@127.0.0.1:27020/beacon?authSource=admin" --file /tmp/filtering_terms.json --collection filtering_terms
	docker exec beacon python beacon/reindex.py
	docker exec beacon python beacon/db/extract_filtering_terms.py