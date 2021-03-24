.PHONY: help check-queue
.DEFAULT_GOAL := help
SHELL := /bin/bash

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

ifneq (,$(wildcard ./.env))
    include .env
    export
endif

help:  ## get help
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

check-queue: ## check queue size
	python3 -c "import redis; redis_conn = redis.from_url('redis://localhost:6379'); print(redis_conn.llen('queue'))"

delete-queue: ## check queue size
	python3 -c "import redis; redis_conn = redis.from_url('redis://localhost:6379'); print(redis_conn.delete('queue'))"

pop-queue: ## check queue size
	python3 -c "import redis; redis_conn = redis.from_url('redis://localhost:6379'); print(redis_conn.lpop('queue'))"	