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

delete-job: ## delete job from database
ifndef job_id
	$(error job_id is not set)
endif
	psql -d pypkaserver -c "delete from job where job_id = $(job_id);"

restart-services: ## restarts the services
	sudo service pypka-api restart
	sudo service pypka-fastapi restart
	


