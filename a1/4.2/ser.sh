#!/bin/bash
#@ wall_clock_limit = 00:20:00
#@ job_name = nhp-serial
#@ job_type = MPICH
#@ class = fattest
#@ output = ser_$(jobid).out
#@ error = ser_$(jobid).out
#@ node = 1
#@ total_tasks = 1
#@ node_usage = not_shared
#@ energy_policy_tag = lulesh
#@ minimize_time_to_solution = yes
#@ island_count = 1
#@ notification = error
#@ notify_user = phu.nguyen@tum.de
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

./lulesh2.0

