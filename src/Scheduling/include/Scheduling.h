/**
 * @file Scheduling.h
 * @brief Memory-aware scheduling for DAGs
 *
 * Performs scheduling of DAG operations to minimize peak memory
 * usage or the bytes transferred between device and host. Uses a CP-SAT
 * solver (brute-force solver for small graphs for verification).
 */

#ifndef DAG_PARTITIONING_SCHEDULING_H
#define DAG_PARTITIONING_SCHEDULING_H

#include "BruteForceSolver.h"
#include "CPSATMemorySolver.h"
#include "CPSATTransferMemorySolver.h"
#include "HyperGraph.h"
#include "Objective.h"
#include "PartitionFootprints.h"
#include "RPOScheduler.h"
#include "Scheduler.h"
#include "Utils.h"

#endif // DAG_PARTITIONING_SCHEDULING_H
